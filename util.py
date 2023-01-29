def float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0