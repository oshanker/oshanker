#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 08:32:32 2024

@author: uorugant
"""

def bunch_points(points, start, interval_length):
    # Sort the points
    points.sort()

    # Initialize variables
    current_interval_start = start
    current_interval_end = current_interval_start + interval_length
    intervals = []
    count = 0

    # Iterate through points
    for point in points:
        # If point falls within current interval
        if point < current_interval_end:
            count = count + 1
            continue
        # If point is outside current interval, start a new interval
        intervals.append((current_interval_start, current_interval_end, count))
        current_interval_start = current_interval_start + interval_length
        current_interval_end = current_interval_start + interval_length
        count = 1

    # Append the last interval
    intervals.append((current_interval_start, current_interval_end, count))

    return intervals

# Example usage
def runpersonMain():
    points = [1, 2, 3, 5, 7, 9, 10, 12, 14, 16]
    interval_length = 3.0
    intervals = bunch_points(points, 1.0, interval_length)
    print("Intervals:", intervals)


if __name__ == "__main__":
    runpersonMain()
