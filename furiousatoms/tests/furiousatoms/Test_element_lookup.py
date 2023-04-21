from furiousatoms.element_lookup import lookup_num_bonds_by_atomic_number

def test_lookup_num_bonds_by_atomic_number():
    assert lookup_num_bonds_by_atomic_number(1) == 1 #H -> 1
    assert lookup_num_bonds_by_atomic_number(2) == 0 #He -> 0
    assert lookup_num_bonds_by_atomic_number(3) == 1 #Li -> 1
    assert lookup_num_bonds_by_atomic_number(4) == 2 #Be -> 2
    assert lookup_num_bonds_by_atomic_number(5) == 3 #B -> 3
    assert lookup_num_bonds_by_atomic_number(6) <= 4 #C -> 2 or 4
    assert lookup_num_bonds_by_atomic_number(7) <= 4 #N -> varies
    assert lookup_num_bonds_by_atomic_number(8) == 2 #O -> 2
    assert lookup_num_bonds_by_atomic_number(9) == 1 #F -> 1
    assert lookup_num_bonds_by_atomic_number(10) == 0 #Ne -> 0

    assert lookup_num_bonds_by_atomic_number(18) == 0  #Ar -> 0
    assert lookup_num_bonds_by_atomic_number(19) == 1  #K -> 1
    assert lookup_num_bonds_by_atomic_number(20) == 2 #Ca -> 2

    assert lookup_num_bonds_by_atomic_number(26) <= 4 #Fe -> 2,3,4
    assert lookup_num_bonds_by_atomic_number(30) == 4 #Zn -> 2 ... incorrect, but difficult to estimate
    assert lookup_num_bonds_by_atomic_number(32) <= 4 #Ge -> 2 or 4

    #Check entire periodic table for invalid results
    for atomic_number in range(118):
        max_bonds = lookup_num_bonds_by_atomic_number(atomic_number)
        print(atomic_number)
        assert max_bonds >= 0
        assert max_bonds <= 4