# Initialize an empty list to store the rows
import scipy.interpolate as interpolate
import numpy as np


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

    return x_scaled, y


if __name__ == "__main__":
    # Define your file name (or use an absolute path if needed)
    target_file = "zetaValues.txt"

    try:
        # Call the function and unpack the returned arrays
        x_array, y_array = load_zeta_data(target_file)

        print("--- Data successfully loaded ---")
        print("X Array:", x_array)
        print("Y Array:", y_array)

    except FileNotFoundError:
        print(
            f"Error: Could not find '{target_file}'. Please check the file path."
        )
