def calculate_difference(num1, num2):
    """
    Calculate the difference between two numbers.

    Parameters:
    num1 (int or float): The first number.
    num2 (int or float): The second number.

    Returns:
    float: Difference between num1 and num2.
    """
    try:
        difference = num1 - num2
        return difference
    except TypeError:
        print("Please provide valid numbers (integers or floats).")


# Example usage:
def runsubract():
	number1 = 10
	number2 = 5
	result = calculate_difference(number1, number2)
	print("Difference:", result)

if __name__ == "__main__":
    runsubract()
