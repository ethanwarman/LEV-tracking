# LEV-tracking
This code calculates leading edge vortex over a half cycle and calculates the centroid of the vortex and circulation strength.

- `LEV_tracking.m` expects input data to be organised in folders by spanwise location (`0.5m/`, `1m/`, `1.5m/` etc).  
- Within each folder:
  - `data-XX.X.csv` = raw vorticity fields at time `XX.X` seconds  
  - `airfoil-XX.X.csv` = airfoil coordinates at the same time  
## ⚙️ How the Script Works

1. **Load vorticity field** (`data-XX.X.csv`) and airfoil geometry (`airfoil-XX.X.csv`).
2. **Interpolate** the scattered vorticity values onto a uniform Cartesian grid.
3. **Mask** the airfoil region and filter weak vorticity values below a chosen threshold.
4. **Identify vortices** using connected regions and select the largest as the LEV.
5. **Compute:**
   - **Centroid position** (weighted by vorticity magnitude).  
   - **Circulation** (integrated vorticity strength).
6. **Track the LEV** across the chosen half-cycle.
7. **Export results** to `.csv` and generate plots for centroid trajectory and circulation evolution.

### 1. Centroid of the Vortex  
The LEV centroid \((x_c, y_c)\) is computed as the **vorticity-weighted average** of all valid vortex points:

\[
x_c = \frac{\sum_i x_i \, |\omega_i|}{\sum_i |\omega_i|}, 
\qquad
y_c = \frac{\sum_i y_i \, |\omega_i|}{\sum_i |\omega_i|}
\]

- \(x_i, y_i\): grid coordinates of the vortex region  
- \(\omega_i\): spanwise vorticity at \((x_i,y_i)\)  
- weights are taken as \(|\omega_i|\) to bias towards strong vortex cores  

