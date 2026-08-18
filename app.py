"""
========================================================================================
PetroRock: Reservoir Quality & SCAL Evaluator
----------------------------------------------------------------------------------------
AI DEVELOPMENT DOCUMENTATION & METHODOLOGY:
1. AI Tools Utilized:
   - Google Gemini 3.7 Flash & Google AI Studio Agentic Engineering Workflow.
   - Code validation for petrophysical formula correctness and numerical stability.

2. Key Prompts Executed During Architecture & Design:
   - Prompt 1: "Design a high-precision reservoir engineering calculation engine for
     Reservoir Quality Index (RQI), Pore Volume-to-Grain Volume Ratio (phi_z), Flow Zone
     Indicator (FZI), and Winland R35 pore throat radius classification with strict
     unit handling."
   - Prompt 2: "Formulate Leverett J-Function capillary pressure Pc(Sw) curves coupled
     with Archie water saturation Sw modeling, incorporating interfacial tension (sigma),
     contact angle (theta), and rock permeability-porosity scaling."
   - Prompt 3: "Implement robust edge-case validation, division-by-zero safeguards,
     interactive dual-axis/faceted Plotly visualizations, and synthetic multi-plug core
     analysis with instant CSV export."

3. Manual Edge-Case Verifications & Numerical Safeguards:
   - Porosity Range Validation: Enforced strict bounds (0 < phi < 1.0 fraction or 0.1% to 99.9%).
     Protected phi_z = phi / (1 - phi) from singularity at phi >= 1.0.
   - Permeability Safeguards: Protected log10(k) and sqrt(k/phi) against negative and zero
     permeability values (k > 0 enforced, minimum epsilon = 1e-6).
   - Archie Sw Clamping: Clamped water saturation between 0.0 and 1.0 (0% - 100%) to prevent
     non-physical saturation values when Rt or phi are extreme.
   - Unit Consistency: Explicit conversion between percentage porosity and fractional
     porosity, field units conversion factor (0.21645 for psi, dynes/cm, mD, and fraction).
========================================================================================
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION & THEMING
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="PetroRock: Reservoir Quality & SCAL Evaluator",
    page_icon="🪨",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for polished petrophysical engineering aesthetics
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E293B;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.05rem;
        color: #64748B;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 16px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .rock-badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 16px;
        font-weight: 600;
        font-size: 0.85rem;
    }
</style>
""", unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# CORE PETROPHYSICAL CALCULATION ENGINE
# -----------------------------------------------------------------------------
def calculate_rqi(perm_md: float, phi_frac: float) -> float:
    """Calculates Reservoir Quality Index (RQI) in micrometers."""
    if phi_frac <= 0 or perm_md <= 0:
        return 0.0
    return 0.0314 * np.sqrt(perm_md / phi_frac)

def calculate_phi_z(phi_frac: float) -> float:
    """Calculates Pore Volume to Grain Volume Ratio (phi_z)."""
    if phi_frac <= 0 or phi_frac >= 1.0:
        return 0.0
    return phi_frac / (1.0 - phi_frac)

def calculate_fzi(rqi: float, phi_z: float) -> float:
    """Calculates Flow Zone Indicator (FZI) in micrometers."""
    if phi_z <= 0:
        return 0.0
    return rqi / phi_z

def calculate_winland_r35(perm_md: float, phi_pct: float) -> float:
    """
    Calculates Winland R35 pore throat radius at 35% mercury saturation (microns).
    Equation: log10(R35) = 0.732 + 0.588 * log10(k) - 0.864 * log10(phi_pct)
    """
    if perm_md <= 0 or phi_pct <= 0:
        return 0.0
    log_r35 = 0.732 + 0.588 * np.log10(perm_md) - 0.864 * np.log10(phi_pct)
    return float(10 ** log_r35)

def classify_winland(r35_microns: float) -> tuple[str, str, str]:
    """
    Classifies pore system into Winland Rock Types.
    Returns: (Category Name, Description, Hex Color)
    """
    if r35_microns > 10.0:
        return "Megaporous", "Excellent reservoir quality; coarse pore network with high flow capacity (>10 µm)", "#059669"
    elif 2.5 <= r35_microns <= 10.0:
        return "Macroporous", "Good to high reservoir quality; moderate to high flow capacity (2.5 - 10 µm)", "#2563EB"
    elif 0.5 <= r35_microns < 2.5:
        return "Mesoporous", "Moderate reservoir quality; fine pore network with moderate flow capacity (0.5 - 2.5 µm)", "#D97706"
    else:
        return "Microporous", "Low to tight reservoir quality; micro-pore throats dominated by capillary trap (<0.5 µm)", "#DC2626"

def calculate_archie_sw(a: float, rw: float, rt: float, phi_frac: float, m: float, n: float) -> float:
    """
    Calculates Water Saturation (Sw) using Archie's classic equation.
    Sw = ((a * Rw) / (Rt * phi^m))^(1/n)
    Clamped strictly between 0.0 and 1.0.
    """
    if phi_frac <= 0 or rt <= 0 or rw <= 0 or a <= 0 or n <= 0:
        return 1.0
    val = (a * rw) / (rt * (phi_frac ** m))
    sw = val ** (1.0 / n)
    return float(np.clip(sw, 0.0, 1.0))

def generate_leverett_curve(perm_md: float, phi_frac: float, sigma_dynes: float, theta_deg: float, points: int = 100):
    """
    Generates Leverett J(Sw) and corresponding Capillary Pressure Pc(Sw) curves.
    Standard Leverett J-function empirical representation: J(Sw) = 0.21645 * (Pc / (sigma * cos(theta))) * sqrt(k / phi)
    Pc (psi) = [J(Sw) * sigma * cos(theta)] / [0.21645 * sqrt(k / phi)]
    """
    # Sw from irreducible water saturation (e.g. 0.08) to 1.0
    sw_array = np.linspace(0.08, 1.0, points)
    
    # Power-law synthetic Leverett J(Sw) normalized model: J(Sw) = 0.5 * ( (Sw - 0.05)/(0.95) )^(-0.7)
    sw_norm = np.clip((sw_array - 0.05) / 0.95, 0.01, 1.0)
    j_sw = 0.45 * (sw_norm ** (-0.75))
    
    # Capillary pressure conversion (Field Units: Pc in psi, sigma in dynes/cm, k in mD, phi fractional)
    cos_theta = np.cos(np.radians(theta_deg))
    if cos_theta <= 0:
        cos_theta = 0.01 # Prevent negative/zero capillary pressure scaling
        
    sqrt_k_phi = np.sqrt(max(perm_md, 1e-4) / max(phi_frac, 1e-4))
    # Conversion factor C = 0.21645 (psi to dynes/cm and mD units)
    pc_psi = (j_sw * sigma_dynes * cos_theta) / (0.21645 * sqrt_k_phi)
    
    return pd.DataFrame({
        "Water_Saturation_Sw": sw_array,
        "Leverett_J_Function": j_sw,
        "Capillary_Pressure_Pc_psi": pc_psi
    })

# -----------------------------------------------------------------------------
# SIDEBAR CONTROLS & INPUT VALIDATION
# -----------------------------------------------------------------------------
st.sidebar.title("⚙️ Petrophysical Parameters")
st.sidebar.markdown("Configure core plug properties and SCAL lab measurements.")

# 1. Primary Core Properties
st.sidebar.subheader("1. Core Plug Petrophysics")
phi_pct = st.sidebar.slider(
    "Porosity φ (%)",
    min_value=1.0,
    max_value=35.0,
    value=18.5,
    step=0.1,
    help="Core helium porosity measured at ambient or reservoir stress conditions (1.0% to 35.0%)."
)
perm_md = st.sidebar.number_input(
    "Permeability k (mD)",
    min_value=0.01,
    max_value=5000.0,
    value=125.0,
    step=1.0,
    format="%.2f",
    help="Absolute air/brine Klinkenberg-corrected permeability (0.01 to 5000.0 mD)."
)
lithology = st.sidebar.selectbox(
    "Reservoir Lithology",
    options=["Sandstones", "Carbonates (Limestone/Dolomite)", "Shaly Sandstone", "Tight Carbonate"],
    index=0,
    help="Dominant mineralogical fabric influencing pore throat distribution and Archie exponents."
)

# 2. SCAL & Interfacial Properties
st.sidebar.subheader("2. SCAL & Fluid Interfacial")
sigma_dynes = st.sidebar.slider(
    "Interfacial Tension σ (dynes/cm)",
    min_value=10.0,
    max_value=80.0,
    value=48.0,
    step=1.0,
    help="Oil-brine or air-mercury interfacial tension (dynes/cm or mN/m)."
)
theta_deg = st.sidebar.slider(
    "Contact Angle θ (degrees)",
    min_value=0.0,
    max_value=85.0,
    value=30.0,
    step=1.0,
    help="Wettability contact angle (0° = strongly water-wet, higher = mixed-wet)."
)

# 3. Archie Saturation Exponents
st.sidebar.subheader("3. Electrical Archie Parameters")
col_arch1, col_arch2 = st.sidebar.columns(2)
with col_arch1:
    archie_a = st.number_input("Tortuosity (a)", min_value=0.5, max_value=2.0, value=1.0, step=0.05)
    archie_m = st.number_input("Cementation (m)", min_value=1.0, max_value=3.0, value=2.0, step=0.05)
    archie_n = st.number_input("Saturation Exp (n)", min_value=1.0, max_value=3.5, value=2.0, step=0.05)
with col_arch2:
    rw_val = st.number_input("Rw (ohm·m)", min_value=0.005, max_value=10.0, value=0.05, step=0.01, format="%.3f")
    rt_val = st.number_input("Rt (ohm·m)", min_value=0.1, max_value=1000.0, value=25.0, step=1.0)

# -----------------------------------------------------------------------------
# ERROR CHECKING & INPUT GUARDS
# -----------------------------------------------------------------------------
has_error = False
if phi_pct <= 0.0 or phi_pct >= 100.0:
    st.error("❌ Invalid Porosity: Porosity must be strictly between 0% and 100%. Please adjust sidebar inputs.")
    has_error = True

if perm_md <= 0.0:
    st.error("❌ Invalid Permeability: Permeability must be strictly greater than 0 mD.")
    has_error = True

if rt_val <= 0 or rw_val <= 0:
    st.warning("⚠️ Warning: Formation or brine resistivity must be positive non-zero values.")
    has_error = True

if has_error:
    st.stop()

# -----------------------------------------------------------------------------
# COMPUTATIONS
# -----------------------------------------------------------------------------
phi_frac = phi_pct / 100.0
rqi_val = calculate_rqi(perm_md, phi_frac)
phi_z_val = calculate_phi_z(phi_frac)
fzi_val = calculate_fzi(rqi_val, phi_z_val)
r35_val = calculate_winland_r35(perm_md, phi_pct)
rock_type, rock_desc, rock_color = classify_winland(r35_val)
sw_val = calculate_archie_sw(archie_a, rw_val, rt_val, phi_frac, archie_m, archie_n)

# -----------------------------------------------------------------------------
# MAIN DASHBOARD PRESENTATION
# -----------------------------------------------------------------------------
st.markdown('<div class="main-header">PetroRock: Reservoir Quality & SCAL Evaluator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Advanced Petrophysical Core Characterization, Hydraulic Flow Unit Zonation & Leverett Capillary Pressure Modeling</div>', unsafe_allow_html=True)

# Expandable User Instructions & Calculation Guide
with st.expander("📖 User Instructions & Petrophysical Calculation Guide", expanded=False):
    st.markdown("""
    ### Theoretical Foundation & Calculation Equations
    This application automates classic and modern reservoir engineering protocols for Core Analysis and Special Core Analysis (SCAL):
    
    1. **Reservoir Quality Index (RQI)**:
       $$RQI = 0.0314 \\times \\sqrt{\\frac{k}{\\phi}}$$
       *Where $k$ is permeability in mD, and $\\phi$ is fractional porosity. $RQI$ represents the hydraulic radius of the porous system (µm).*
       
    2. **Pore Volume-to-Grain Volume Ratio ($\\phi_z$)**:
       $$\\phi_z = \\frac{\\phi}{1 - \\phi}$$
       *Represents normalized pore volume available for fluid conductance per unit grain matrix volume.*

    3. **Flow Zone Indicator (FZI)**:
       $$FZI = \\frac{RQI}{\\phi_z}$$
       *A fundamental parameter characterizing hydraulic flow units (HFU) independent of depth and sample dimensions.*

    4. **Winland $R_{35}$ Pore Throat Radius**:
       $$\\log_{10}(R_{35}) = 0.732 + 0.588 \\log_{10}(k) - 0.864 \\log_{10}(\\phi_{\\%})$$
       *Empirical correlation establishing the pore aperture radius (microns) at 35% mercury saturation.*
       
    5. **Archie Water Saturation ($S_w$)**:
       $$S_w = \\left( \\frac{a \\cdot R_w}{R_t \\cdot \\phi^m} \\right)^{1/n}$$
       *Calculates uninvaded formation water saturation from resistivity logs and core electrical parameters.*

    6. **Leverett $J(S_w)$ Function & Capillary Pressure ($P_c$)**:
       $$J(S_w) = 0.21645 \\cdot \\frac{P_c}{\\sigma \\cos\\theta} \\sqrt{\\frac{k}{\\phi}} \\implies P_c(S_w) = \\frac{J(S_w) \\cdot \\sigma \\cos\\theta}{0.21645 \\sqrt{k/\\phi}}$$
       *Normalizes heterogeneous capillary pressure curves into a universal dimensionless function for reservoir-scale saturation-height modeling.*
    """)

# Metric Cards Row
m_col1, m_col2, m_col3, m_col4 = st.columns(4)
with m_col1:
    st.metric(
        label="Reservoir Quality Index (RQI)",
        value=f"{rqi_val:.3f} µm",
        delta=f"k/φ: {(perm_md/phi_frac):.1f}"
    )
with m_col2:
    st.metric(
        label="Flow Zone Indicator (FZI)",
        value=f"{fzi_val:.3f} µm",
        delta=f"φ_z: {phi_z_val:.3f}"
    )
with m_col3:
    st.metric(
        label="Winland R35 Radius",
        value=f"{r35_val:.2f} µm",
        delta=rock_type
    )
with m_col4:
    st.metric(
        label="Archie Water Saturation (Sw)",
        value=f"{sw_val * 100.0:.1f} %",
        delta=f"Hydrocarbon So: {((1.0 - sw_val) * 100.0):.1f}%"
    )

# Rock Classification Banner
st.markdown(f"""
<div style="background-color: {rock_color}15; border-left: 5px solid {rock_color}; padding: 12px 18px; border-radius: 6px; margin-top: 10px; margin-bottom: 20px;">
    <strong style="color: {rock_color}; font-size: 1.05rem;">Winland Rock Classification: {rock_type} (R35 = {r35_val:.2f} µm)</strong>
    <p style="margin: 4px 0 0 0; color: #334155; font-size: 0.92rem;">{rock_desc} | Lithology: <strong>{lithology}</strong></p>
</div>
""", unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# PLOTLY INTERACTIVE VISUALIZATION
# -----------------------------------------------------------------------------
st.subheader("📈 SCAL Capillary Pressure & Leverett J-Function Curves")

scal_df = generate_leverett_curve(perm_md, phi_frac, sigma_dynes, theta_deg)

# Create 2-subplot figure
fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=(
        "Capillary Pressure Pc (psi) vs Water Saturation Sw",
        "Dimensionless Leverett J(Sw) vs Water Saturation Sw"
    ),
    horizontal_spacing=0.1
)

# Plot 1: Pc vs Sw
fig.add_trace(
    go.Scatter(
        x=scal_df["Water_Saturation_Sw"] * 100,
        y=scal_df["Capillary_Pressure_Pc_psi"],
        mode="lines",
        name="Pc Curve (psi)",
        line=dict(color="#2563EB", width=3),
        hovertemplate="<b>Sw</b>: %{x:.1f}%<br><b>Pc</b>: %{y:.2f} psi<extra></extra>"
    ),
    row=1, col=1
)

# Highlight Current Archie Sw on Pc Curve
current_pc = np.interp(sw_val, scal_df["Water_Saturation_Sw"], scal_df["Capillary_Pressure_Pc_psi"])
fig.add_trace(
    go.Scatter(
        x=[sw_val * 100],
        y=[current_pc],
        mode="markers",
        name="Current Reservoir Operating Point",
        marker=dict(color="#DC2626", size=12, symbol="diamond"),
        hovertemplate="<b>In-situ Sw</b>: %{x:.1f}%<br><b>Equiv Pc</b>: %{y:.2f} psi<extra></extra>"
    ),
    row=1, col=1
)

# Plot 2: Leverett J(Sw) vs Sw
fig.add_trace(
    go.Scatter(
        x=scal_df["Water_Saturation_Sw"] * 100,
        y=scal_df["Leverett_J_Function"],
        mode="lines",
        name="Leverett J(Sw)",
        line=dict(color="#059669", width=3, dash="dash"),
        hovertemplate="<b>Sw</b>: %{x:.1f}%<br><b>J(Sw)</b>: %{y:.3f}<extra></extra>"
    ),
    row=1, col=2
)

fig.update_xaxes(title_text="Water Saturation Sw (%)", row=1, col=1, gridcolor="#E2E8F0")
fig.update_yaxes(title_text="Capillary Pressure Pc (psi)", row=1, col=1, gridcolor="#E2E8F0")
fig.update_xaxes(title_text="Water Saturation Sw (%)", row=1, col=2, gridcolor="#E2E8F0")
fig.update_yaxes(title_text="Leverett J-Function J(Sw) [-]", row=1, col=2, gridcolor="#E2E8F0")

fig.update_layout(
    height=480,
    margin=dict(l=40, r=40, t=60, b=40),
    template="plotly_white",
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
)

st.plotly_chart(fig, use_container_width=True)

# -----------------------------------------------------------------------------
# SYNTHETIC MULTI-CORE PLUG ANALYSIS & CSV EXPORT
# -----------------------------------------------------------------------------
st.subheader("🧪 Synthetic Core Plugs Suite & Flow Unit Zonation")
st.markdown("Multi-sample core plug batch evaluation simulating representative reservoir flow unit intervals.")

# Generate synthetic multi-plug core dataset around current operating point
np.random.seed(42)
num_plugs = 8
plug_ids = [f"CP-{101 + i}" for i in range(num_plugs)]
sample_depths = [2450.0 + (i * 2.5) for i in range(num_plugs)]

# Variance around current inputs
phi_variations = np.clip(np.random.normal(loc=phi_pct, scale=3.5, size=num_plugs), 3.0, 34.0)
k_multiplier = np.random.lognormal(mean=0, sigma=0.65, size=num_plugs)
perm_variations = np.clip(perm_md * k_multiplier, 0.05, 4500.0)

table_rows = []
for pid, depth, p_pct, k_val in zip(plug_ids, sample_depths, phi_variations, perm_variations):
    p_f = p_pct / 100.0
    r_qi = calculate_rqi(k_val, p_f)
    p_z = calculate_phi_z(p_f)
    f_zi = calculate_fzi(r_qi, p_z)
    r_35 = calculate_winland_r35(k_val, p_pct)
    w_type, _, _ = classify_winland(r_35)
    s_w = calculate_archie_sw(archie_a, rw_val, rt_val, p_f, archie_m, archie_n)
    
    table_rows.append({
        "Plug ID": pid,
        "Depth (m)": round(depth, 1),
        "Porosity φ (%)": round(p_pct, 2),
        "Permeability k (mD)": round(k_val, 2),
        "RQI (µm)": round(r_qi, 3),
        "φ_z": round(p_z, 4),
        "FZI (µm)": round(f_zi, 3),
        "Winland R35 (µm)": round(r_35, 2),
        "Pore Type": w_type,
        "Archie Sw (%)": round(s_w * 100.0, 1)
    })

core_df = pd.DataFrame(table_rows)

col_tbl, col_dl = st.columns([3, 1])
with col_tbl:
    st.dataframe(
        core_df,
        use_container_width=True,
        hide_index=True
    )
with col_dl:
    st.markdown("#### 📥 Data Export")
    st.markdown("Export complete petrophysical evaluation report for reservoir modeling.")
    csv_data = core_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Core Data (CSV)",
        data=csv_data,
        file_name="PetroRock_Core_Evaluation_Report.csv",
        mime="text/csv",
        use_container_width=True
    )
    st.info(f"Average Reservoir FZI: **{core_df['FZI (µm)'].mean():.2f} µm**\n\nDominant Flow Unit: **{core_df['Pore Type'].mode()[0]}**")
