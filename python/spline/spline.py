import numpy as np
import scipy.interpolate as interpolate

# pip install numpy scipy matplotlib

# 1. Example data crossing the x-axis
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([1, -1, 2, -2, 1, -1])

# 2. Fit the spline
spline = interpolate.make_interp_spline(x, y, k=3)

# 3. Convert to PPoly and find exact roots
ppoly = interpolate.PPoly.from_spline(spline)
zeros = ppoly.roots()

print("The spline crosses zero at x positions:", zeros)
# Output: [0.09843471, 1.2305, 2.5000, 3.7694, 4.9015]

# To find where spline(x) == 5, solve: spline(x) - 5 == 0
ppoly_shifted = interpolate.PPoly.from_spline(spline)
ppoly_shifted.c[-1] -= 5  # Shift the constant term of the polynomial segments
target_points = ppoly_shifted.roots()
print("The spline crosses 5 at x positions:", target_points)

y_at_single_point = spline(-0.26587818)
print(f"The y-value at x=-0.26587818 is: {y_at_single_point}")

# At multiple specific points:
target_points = np.array([-0.26587818, 0.09843471])
y_values = spline(target_points)
print(f"The y-values at {target_points} are: {y_values}")

slope_at_point = spline(2.5, nu=1)
print(f"The slope at x=2.5 is: {slope_at_point}")

derivative_spline = spline.derivative(nu=1)

# 4. Convert the derivative to a piecewise polynomial to find exact zeros (where slope = 0)
ppoly_deriv = interpolate.PPoly.from_spline(derivative_spline)
critical_points = ppoly_deriv.roots()

# Filter critical points to ensure they stay strictly within your data bounds
critical_points = critical_points[
    (critical_points >= x.min()) & (critical_points <= x.max())
]

# 5. Separate them into High Peaks (maxima) and Low Valleys (minima)
peaks = []
valleys = []

# Get the second derivative to test if it's a max or min
second_deriv_spline = spline.derivative(nu=2)

for pt in critical_points:
    y_val = float(spline(pt))
    slope_change = float(second_deriv_spline(pt))
    
    if slope_change < 0:
        peaks.append((pt, y_val))     # Concave down = Peak
    elif slope_change > 0:
        valleys.append((pt, y_val))   # Concave up = Valley

# --- PRINT THE RESULTS ---
print(" LOCAL MAXIMA (HIGH PEAKS):")
for pt, y_val in peaks:
    print(f"  x = {pt:.4f}, y = {y_val:.4f}")

print("\n LOCAL MINIMA (VALLEYS):")
for pt, y_val in valleys:
    print(f"  x = {pt:.4f}, y = {y_val:.4f}")


