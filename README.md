# LEV-tracking
This code calculates leading edge vortex over a half cycle and calculates the centroid of the vortex and circulation strength.

- `LEV_tracking.m` expects input data to be organised in folders by spanwise location (`0.5m/`, `1m/`, `1.5m/` etc).  
- Within each folder:
  - `data-XX.X.csv` = raw vorticity fields at time `XX.X` seconds  
  - `airfoil-XX.X.csv` = airfoil coordinates at the same time  

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

