import csv
from io import StringIO

DEFAULT_MASS = 0.0

elems_csv = """\
number,symbol,name,mass,valency
1,H,hydrogen,1.0079,1
2,He,helium,4.0026,0
3,Li,lithium,6.941,1
4,Be,beryllium,9.0122,2
5,B,boron,10.811,3
6,C,carbon,12.0107,4
7,N,nitrogen,14.0067,-3
8,O,oxygen,15.9994,-2
9,F,fluorine,18.9984,-1
10,Ne,neon,20.1797,0
11,Na,sodium,22.9897,1
12,Mg,magnesium,24.305,2
13,Al,aluminium,26.9815,3
14,Si,silicon,28.0855,4
15,P,phosphorus,30.9738,5
16,S,sulfur,32.065,6
17,Cl,chlorine,35.453,-1
18,Ar,argon,39.948,0
19,K,potassium,39.0983,1
20,Ca,calcium,40.078,2
21,Sc,scandium,44.9559,3
22,Ti,titanium,47.867,4
23,V,vanadium,50.9415,5
24,Cr,chromium,51.9961,3
25,Mn,manganese,54.938,2
26,Fe,iron,55.845,2
27,Co,cobalt,58.9332,3
28,Ni,nickel,58.6934,2
29,Cu,copper,63.546,2
30,Zn,zinc,65.39,2
31,Ga,gallium,69.723,3
32,Ge,germanium,72.64,-4
33,As,arsenic,74.9216,-3
34,Se,selenium,78.96,-2
35,Br,bromine,79.904,-1
36,Kr,krypton,83.8,2
37,Rb,rubidium,85.4678,1
38,Sr,strontium,87.62,2
39,Y,yttrium,88.9059,3
40,Zr,zirconium,91.224,4
41,Nb,niobium,92.9064,5
42,Mo,molybdenum,95.94,4
43,Tc,technetium,98,4
44,Ru,ruthenium,101.07,3
45,Rh,rhodium,102.9055,3
46,Pd,palladium,106.42,2
47,Ag,silver,107.8682,1
48,Cd,cadmium,112.411,2
49,In,indium,114.818,3
50,Sn,tin,118.71,-4
51,Sb,antimony,121.76,0
52,Te,tellurium,127.6,-3
53,I,iodine,126.9045,-2
54,Xe,xenon,131.293,2
55,Cs,caesium,132.9055,1
56,Ba,barium,137.327,2
57,La,lanthanum,138.9055,3
58,Ce,cerium,140.116,3
59,Pr,praseodymium,140.9077,3
60,Nd,neodymium,144.24,3
61,Pm,promethium,145,3
62,Sm,samarium,150.36,3
63,Eu,europium,151.964,3
64,Gd,gadolinium,157.25,3
65,Tb,terbium,158.9253,3
66,Dy,dysprosium,162.5,3
67,Ho,holmium,164.9303,3
68,Er,erbium,167.259,3
69,Tm,thulium,168.9342,3
70,Yb,ytterbium,173.04,3
71,Lu,lutetium,174.967,3
72,Hf,hafnium,178.49,4
73,Ta,tantalum,180.9479,5
74,W,tungsten,183.84,4
75,Re,rhenium,186.207,4
76,Os,osmium,190.23,4
77,Ir,iridium,192.217,3
78,Pt,platinum,195.078,2
79,Au,gold,196.9665,3
80,Hg,mercury,200.59,1
81,Tl,thallium,204.3833,1
82,Pb,lead,207.2,2
83,Bi,bismuth,208.9804,3
84,Po,polonium,209,-2
85,At,astatine,210,-1
86,Rn,radon,222,0
87,Fr,francium,223,1
88,Ra,radium,226,2
89,Ac,actinium,227,3
90,Th,thorium,232.0381,4
91,Pa,protactinium,231.0359,5
92,U,uranium,238.0289,6
93,Np,Neptunium,237.064,5
94,Pu,Plutonium,244.064,3
95,Am,Americium,243.061,3
96,Cm,Curium,247.070,3
97,Bk,Berkelium,247.070,3
98,Cf,Californium,251.080,3
99,Es,Einsteinium,254,3
100,Fm,Fermium,257.0953,3
101,Md,Mendelevium,258.1,3
102,No,Nobelium,259.101,2
103,Lr,Lawrencium,262,3
104,Rf,Rutherfordium,261,4
105,Db,Dubnium,262,5
106,Sg,Seaborgium,266,6
107,Bh,Bohrium,264,7
108,Hs,Hassium,269,8
109,Mt,Meitnerium,278,0
110,Ds,Darmstadtium,281,0
111,Rg,Roentgenium,280,0
112,Cn,Copernicium,285,2
113,Nh,Nihonium,286,0
114,Fl,Fleronium,289,0
115,Mc,Moscovium,289,0
116,Lv,Livermorium,293,0
117,Ts,Tennessine,294,0
118,Og,Oganesson,294,0
"""


def lookup_symbol_by_mass(mass: float):
    for line in elems_csv.split()[1:]:
        fields = line.split(",") #Format: number,symbol,name,mass,valency
        m = float(fields[3])

        #Search for element with mass that's up to 0.1 different from input.
        #If it can't find one, try again with 1.0. 
        for tolerance in (0.1, 1.0):
            if abs(mass - m) < tolerance:
                return fields[1]
    
    raise ValueError("No element exists with mass", mass)


def lookup_mass_by_symbol(symbol: float):
    for line in elems_csv.split()[1:]:
        fields = line.split(",") #Format: number,symbol,name,mass,valency
        if fields[1] == symbol:
            return float(fields[3])
    
    raise ValueError("No element exists with symbol", symbol)


def lookup_num_bonds_by_atomic_number(num: int):
    '''
    Returns the maximum number of bonds that the element can form.
    '''
    #K, L, M, and N are the names of shells in chemistry, and these
    #are the max number of electrons each can hold.
    K = 2
    L = 8
    M = 8
    N = 18
    DEFAULT_NUM_BONDS = 4

    '''
    Fill up the smallest shells first.
    The number of electrons remaining tells the atom's valence.
    '''
    if num > K:
        num -= K
    else:
        if num == 2:
            return 0 #helium
        else:
            return 1 #hydrogen
    
    if num >= L:
        num -= L
    else:
        if num > 4:
            return L - num
        return num
    
    if num > M:
        num -= M
    else:
        if num > 4:
            return M - num
        return num
    
    if num > N:
        num -= N
    else:
        if num > 4:
            num = DEFAULT_NUM_BONDS #Uncertain for high atomic numbers
        return num
    
    return DEFAULT_NUM_BONDS


if __name__ == "__main__":
    #Check entire periodic table for invalid results
    for atomic_number in range(118):
        max_bonds = lookup_num_bonds_by_atomic_number(atomic_number)
        print(atomic_number)
        assert max_bonds >= 0
        assert max_bonds <= 4