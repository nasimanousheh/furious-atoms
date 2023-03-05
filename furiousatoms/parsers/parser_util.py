def float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0
    
def has_bond(bonds, new_bond):
    for i in range(len(bonds)):
        if bonds[i][0] == new_bond[0] and bonds[i][1] == new_bond[1] \
        or bonds[i][1] == new_bond[0] and bonds[i][0] == new_bond[1]:
            return True
    return False