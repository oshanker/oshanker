

def add_numbers(numbers):
    """
    Calculate the sum of numbers in a list.

    Parameters:
    numbers (list of int or float): List of numbers to be summed.

    Returns:
    float: Sum of numbers in the list.
    """
    try:
        total = sum(numbers)
        return total
    except TypeError:
        print("Please provide a list of numbers (integers or floats).")

# Example usage:
def runadd():
	numbers = [1, 2, 3, 4, 5]
	total_sum = add_numbers(numbers)
	print("Sum of numbers:", total_sum)

if __name__ == "__main__":
    runadd()
    
# python3 base_transformer_shanker/add.py
# https://www.freecodecamp.org/news/build-your-first-python-package/
