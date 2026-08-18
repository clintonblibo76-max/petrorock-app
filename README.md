# PetroRock: Reservoir Quality & SCAL Evaluator

PetroRock is a specialized reservoir engineering and petrophysics web application designed for rapid core characterization, hydraulic flow unit (HFU) zonation, and special core analysis (SCAL) capillary pressure modeling. By integrating classic petrophysical frameworks—including Reservoir Quality Index (RQI), Flow Zone Indicator (FZI), Winland $R_{35}$ pore throat classification, Archie water saturation ($S_w$), and Leverett $J(S_w)$ capillary pressure normalization—PetroRock empowers subsurface engineers and geoscientists to evaluate reservoir rock quality, determine pore aperture distributions, and model fluid saturation distributions with rigorous numerical stability.

---

## 🚀 Live Application

- **Streamlit Community Cloud URL**: `https://petrorock-app-bcs2610061209.streamlit.app/` 

---

## 🛠️ Local Installation & Execution

Follow these steps to run PetroRock locally:

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/petrorock-evaluator.git
cd petrorock-evaluator
```

### 2. Create and Activate a Virtual Environment
```bash
# On macOS/Linux
python3 -m venv venv
source venv/bin/activate

# On Windows
python -m venv venv
venv\Scripts\activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Launch the Streamlit App
```bash
streamlit run app.py
```
The application will open automatically in your default browser at `http://localhost:8501`.

---

## 📝 Git Commit Strategy

Below is the recommended 3-stage Git commit sequence for this repository:

1. **`feat: initialize petrophysical calculation engine and core petrorock architecture`**
   - Implements mathematical formulations for RQI, $\phi_z$, FZI, Winland $R_{35}$, Archie $S_w$, and Leverett $J$-function.
   - Sets up parameter edge-case validations and numerical singularity guards.

2. **`feat: build interactive streamlit UI, plotly SCAL curves, and synthetic core plug suite`**
   - Adds responsive sidebar inputs for core properties, fluid interfacial tension, and Archie parameters.
   - Integrates dynamic dual-axis Plotly capillary pressure curves and batch core plug analysis table with CSV export.

3. **`docs: complete README instructions, calculation guides, and production requirements`**
   - Formulates petrophysical equations documentation and expandable calculation guide.
   - Adds `requirements.txt`, deployment instructions for Streamlit Community Cloud, and AI engineering methodology.

---

## 📚 Mathematical Equations Reference

| Parameter | Formula | Physical Meaning |
| :--- | :--- | :--- |
| **Reservoir Quality Index ($RQI$)** | $RQI = 0.0314 \sqrt{\frac{k}{\phi}}$ | Hydraulic radius of pore system ($\mu\text{m}$) |
| **Pore-to-Grain Ratio ($\phi_z$)** | $\phi_z = \frac{\phi}{1 - \phi}$ | Normalized pore volume per grain volume |
| **Flow Zone Indicator ($FZI$)** | $FZI = \frac{RQI}{\phi_z}$ | Hydraulic flow unit indicator ($\mu\text{m}$) |
| **Winland $R_{35}$** | $\log_{10}(R_{35}) = 0.732 + 0.588\log_{10}(k) - 0.864\log_{10}(\phi_{\\%})$ | Pore aperture radius at 35% Hg saturation ($\mu\text{m}$) |
| **Archie Water Saturation ($S_w$)** | $S_w = \left(\frac{a \cdot R_w}{R_t \cdot \phi^m}\right)^{1/n}$ | Uninvaded formation water saturation |
| **Leverett $J$-Function** | $J(S_w) = 0.21645 \frac{P_c}{\sigma \cos\theta}\sqrt{\frac{k}{\phi}}$ | Universal dimensionless capillary pressure scaling |
