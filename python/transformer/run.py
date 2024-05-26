import base_transformer_shanker

numbers = [1, 2, 3, 4, 5]
print( base_transformer_shanker.__doc__)
print(dir(base_transformer_shanker))
total_sum = base_transformer_shanker.add_numbers(numbers)
print("Sum of numbers:", total_sum)

number1 = 10
number2 = 5
result = base_transformer_shanker.calculate_difference(number1, number2)
print("Difference:", result)
