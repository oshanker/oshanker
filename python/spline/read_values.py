# Initialize an empty list to store the rows
import scipy.interpolate as interpolate
import numpy as np

def critical(x, spline):
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


def original():
    data_array = []

    # Open and read the file
    with open("zetaValues.txt", "r") as file:
        # Skip the first line (the header)
        next(file)
        
        # Process each subsequent line
        for line in file:
            # Clean up whitespace and check if the line isn't empty
            line = line.strip()
            if line:
                # Split the line by the comma delimiter
                col1, col2 = line.split(",")
                
                # Convert values to correct data types and append as a pair
                row = [int(col1), float(col2)]
                data_array.append(row)

    # Print the resulting 2D list
    print(data_array)


def load_zeta_data(filename):
    """Reads a comma-separated data file, skips the header,

    and returns separate NumPy arrays for X and Y.
    """
    # Load the 2D array (skips 1 line header, splits by comma)
    data_array = np.genfromtxt(filename, delimiter=",", skip_header=1)

    # Split into 1D arrays for X (column 0) and Y (column 1)
    n_indices = data_array[:, 0]
    y = data_array[:, 1]
    increment_val = 0.24359904590398668
    base_val = 244.02115917156451839965694310614387 - 2 * increment_val

    x_scaled = base_val + (n_indices * increment_val)
    # Loop through both arrays simultaneously and print line by line
    # f-strings ensure high-precision formatting for the double floats
    n = 0
    for x_val, y_val in zip(x_scaled, y):
        n = n + 1
        print(f"n {n} X: {x_val:<20.14f} | Y: {y_val:.16f}")
    return find_zeros(x_scaled, y)

def find_zeros(x_scaled, y):
    # 2. Fit the spline
    spline = interpolate.make_interp_spline(x_scaled, y, k=3)

    # 3. Convert to PPoly and find exact roots
    ppoly = interpolate.PPoly.from_spline(spline)
    zeros = ppoly.roots()

    print("The spline crosses zero at x positions:", zeros)
    slopes = spline(zeros, nu = 1)
    print("der at zero  positions:", slopes)
    for x_val, y_val in zip(zeros, slopes):
        print(f" X: {x_val:<20.14f} | Y: {y_val:.16f}")
    critical(x_scaled, spline)

    return x_scaled, y


if __name__ == "__main__":
    # Define your file name (or use an absolute path if needed)
    target_file = "zetaValues.txt"

    try:
        # Call the function and unpack the returned arrays
        # x_array, y_array = load_zeta_data(target_file)
        x = np.array([-2, -1, 0, 1, 2])
        y = np.array([-6,  0, 0, 0, 6])
        x_array, y_array = find_zeros(x, y)

        print("--- Data successfully loaded ---")
        print("X Array:", x_array)
        print("Y Array:", y_array)

    except FileNotFoundError:
        print(
            f"Error: Could not find '{target_file}'. Please check the file path."
        )
