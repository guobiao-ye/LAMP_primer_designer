def generate_distance_penalties(max_distance,a=20, b=1.01, max_penalty=100):
    """
    Generate distance penalties based on the maximum distance using the custom function 
    f(x) = a * (b^x - 1), and limit the penalty values to a maximum value.

    :param a: Scaling factor for the penalty values.
    :param b: Base of the exponent, b > 1 controls the growth rate.
    :param max_penalty: Maximum allowable penalty value (default is 25).
    :return: A list of distance penalties.
    """
    # Generate penalties using f(x) = a * (b^x - 1) and limit to max_penalty
    penalties = [min(a * (b ** x - 1), max_penalty) for x in range(max_distance)]

    return penalties