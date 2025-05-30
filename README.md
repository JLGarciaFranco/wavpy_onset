# Midsummer Drought (MSD) Timing Toolkit  
*A Haar–wavelet approach for rainy‑season onset, retreat, and midsummer break detection*

---

This repository contains Python utilities to diagnose **rainy‑season onset/retreat** and the **Midsummer Drought (MSD)** in tropical precipitation. The method follows **García‑Franco (2021)**, who adapted a Haar covariance wavelet transform to isolate coherent shifts in daily rainfall intensity.

---
## 📌 Key features

| Feature | Description |
|---------|-------------|
| **Haar covariance transform** | Automated or fixed‑dilation Haar wavelet to highlight step‑like changes in rainfall. |
| **Rainy‑season timing** | Finds first sustained intensification (**onset**) and decline (**retreat**) from the seasonal transform. |
| **MSD timing** | Detects the dry spell between two seasonal peaks by analysing a shorter window within the rainy period. |
| **Diagnostics table** | For each year the code outputs nine metrics (see below). |

### Diagnostic fields (per year)
1. **Year**  
2. **Onset date** — first rainfall intensification  
3. **Retreat date** — final decline of rainy season  
4. **MSD onset date** — start of midsummer break  
5. **MSD end date** — re‑intensification after MSD  
6. **Coef min** — minimum wavelet coefficient before MSD  
7. **Coef max** — peak coefficient after MSD  
8. **Amplitude** = Coef max − Coef min  
9. **JJAS mean precip** — average rainfall between onset and retreat  

Results append to `MSD_obs_table.txt` (CSV‑style).

---
## 🔧 Quick start

### 1 · Install
```bash
pip install numpy pandas xarray matplotlib
```

### 2 · Run a simple analysis
```python
import xarray as xr
from msd_wavelet import calc_msd

# daily precipitation (time, lat, lon)
precip_ds = xr.open_dataset("CHIRPS_1981‑2020_daily.nc")
precip_da = precip_ds["precip"]

# single‑year MSD diagnostics
calc_msd(dataset_name="CHIRPS", year=2015, ds=precip_da)
```
The function appends a line to `MSD_obs_table.txt` and returns the same values as a Python list.

### 3 · Batch loop example
```python
for yr in range(1981, 2021):
    calc_msd("CHIRPS", yr, precip_da)
```

---
## 🗂️ Repository layout
| Path | Purpose |
|------|---------|
| `msd_wavelet.py` | Implementation of `calc_msd()` and lower‑level Haar utilities. |
| `example_data/` | Placeholder for CHIRPS or other daily precipitation datasets. |
| `example_notebook.ipynb` | Walk‑through of end‑to‑end MSD analysis (optional). |
| `MSD_obs_table.txt` | Diagnostics table generated by the script. |

---
## 📖 Reference
> **García‑Franco, J. L.** (2021). *A wavelet transform approach to detect the Midsummer Drought onset and retreat.* **International Journal of Climatology**, 41(6), 3862–3875. <https://doi.org/10.1002/joc.7061>

If you use this code in academic work, please cite the paper above and consider referencing this repository (Zenodo DOI forthcoming).

---
## 📬 Contact
Developed by **Jorge Luis García Franco**  
✉️ [jgcaspark@ciencias.unam.mx](mailto:jgcaspark@ciencias.unam.mx)


