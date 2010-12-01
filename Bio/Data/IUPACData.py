# Information about the IUPAC alphabets

protein_letters = "ACDEFGHIKLMNPQRSTVWY"
extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
#   B = "Asx";  aspartic acid or asparagine (D or N)
#   X = "Xxx";  unknown or 'other' amino acid
#   Z = "Glx";  glutamic acid or glutamine (E or Q)
#   J = "Xle";  leucine or isoleucine (L or I, used in mass-spec)
#   U = "Sec";  selenocysteine
#   O = "Pyl";  pyrrolysine
ambiguous_dna_letters = "GATCRYWSMKHBVDN"
unambiguous_dna_letters = "GATC"
ambiguous_rna_letters = "GAUCRYWSMKHBVDN"
unambiguous_rna_letters = "GAUC"

#   B == 5-bromouridine
#   D == 5,6-dihydrouridine
#   S == thiouridine
#   W == wyosine
extended_dna_letters = "GATCBDSW"

# are there extended forms?
#extended_rna_letters = "GAUCBDSW"

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
ambiguous_rna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
    }

ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }

ambiguous_rna_complement = {
    "A": "U",
    "C": "G",
    "G": "C",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


def _make_ranges(mydict):
    d = {}
    for key, value in mydict.iteritems():
        d[key] = (value, value)
    return d

# From bioperl's SeqStats.pm
unambiguous_dna_weights = {
    "A": 347.,
    "C": 323.,
    "G": 363.,
    "T": 322.,
    }
unambiguous_dna_weight_ranges = _make_ranges(unambiguous_dna_weights)

unambiguous_rna_weights = {
    "A": unambiguous_dna_weights["A"] + 16.,  # 16 for the oxygen
    "C": unambiguous_dna_weights["C"] + 16.,
    "G": unambiguous_dna_weights["G"] + 16.,
    "U": 340.,
}
unambiguous_rna_weight_ranges = _make_ranges(unambiguous_rna_weights)

def _make_ambiguous_ranges(mydict, weight_table):
    range_d = {}
    avg_d = {}
    for letter, values in mydict.iteritems():
        #Following line is a quick hack to skip undefined weights for U and O
        if len(values)==1 and values[0] not in weight_table : continue
        weights = map(weight_table.get, values)
        range_d[letter] = (min(weights), max(weights))
        total_w = 0.0
        for w in weights:
            total_w = total_w + w
        avg_d[letter] = total_w / len(weights)
    return range_d, avg_d

ambiguous_dna_weight_ranges, avg_ambiguous_dna_weights = \
               _make_ambiguous_ranges(ambiguous_dna_values,
                                      unambiguous_dna_weights)

ambiguous_rna_weight_ranges, avg_ambiguous_rna_weights = \
               _make_ambiguous_ranges(ambiguous_rna_values,
                                      unambiguous_rna_weights)

protein_weights = {
    "A": 89.09,
    "C": 121.16,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.18,
    "K": 146.19,
    "L": 131.18,
    "M": 149.21,
    "N": 132.12,
    #"O": 0.0, # Needs to be recorded!
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    #"U": 168.05, # To be confirmed
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19
    }

extended_protein_values = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "J": "IL",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "O": "O",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "U": "U",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    #TODO - Include U and O in the possible values of X?
    #This could alter the extended_protein_weight_ranges ...
    "Y": "Y",
    "Z": "QE",
}
    
protein_weight_ranges = _make_ranges(protein_weights)

extended_protein_weight_ranges, avg_extended_protein_weights = \
               _make_ambiguous_ranges(extended_protein_values,
                                      protein_weights)

# Joao Rodrigues @ 2010
# For Center of Mass Calculation.
# Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
atom_weigths = {
    'Ru': 101.070,
    'Re': 186.207,
    'Rf': 261.110,
    'Ra': 226.030,
    'Rb': 85.468,
    'Rn': 222.020,
    'Rh': 102.906,
    'Be': 9.012,
    'Ba': 137.327,
    'Bh': 264.120,
    'Bi': 208.980,
    'Bk': 247.070,
    'Br': 79.904,
    'H': 1.008,
    'P': 30.974,
    'Os': 190.230,
    'Ge': 72.640,
    'Gd': 157.250,
    'Ga': 69.723,
    'Pr': 140.908,
    'Pt': 195.078,
    'Pu': 244.060,
    'C': 12.011,
    'Pb': 207.200,
    'Pa': 231.036,
    'Pd': 106.420,
    'Xe': 131.293,
    'Po': 208.980,
    'Pm': 145.000,
    'Hs': 269.130,
    'Ho': 164.930,
    'Hf': 178.490,
    'Hg': 200.590,
    'He': 4.003,
    'Md': 258.100,
    'Mg': 24.305,
    'K': 39.098,
    'Mn': 54.938,
    'O': 15.999,
    'Mt': 268.140,
    'S': 32.065,
    'W': 183.840,
    'Zn': 65.390,
    'Eu': 151.964,
    'Es': 252.080,
    'Er': 167.259,
    'Ni': 58.693,
    'No': 259.100,
    'Na': 22.990,
    'Nb': 92.906,
    'Nd': 144.240,
    'Ne': 20.180,
    'Np': 237.050,
    'Fr': 223.020,
    'Fe': 55.845,
    'Fm': 257.100,
    'B': 10.811,
    'F': 18.998,
    'Sr': 87.620,
    'N': 14.007,
    'Kr': 83.800,
    'Si': 28.085,
    'Sn': 118.710,
    'Sm': 150.360,
    'V': 50.941,
    'Sc': 44.956,
    'Sb': 121.760,
    'Sg': 266.120,
    'Se': 78.960,
    'Co': 58.933,
    'Cm': 247.070,
    'Cl': 35.453,
    'Ca': 40.078,
    'Cf': 251.080,
    'Ce': 140.116,
    'Cd': 112.411,
    'Tm': 168.934,
    'Cs': 132.905,
    'Cr': 51.996,
    'Cu': 63.546,
    'La': 138.905,
    'Li': 6.941,
    'Tl': 204.383,
    'Lu': 174.967,
    'Lr': 262.110,
    'Th': 232.038,
    'Ti': 47.867,
    'Te': 127.600,
    'Tb': 158.925,
    'Tc': 98.000,
    'Ta': 180.948,
    'Yb': 173.040,
    'Db': 262.110,
    'Dy': 162.500,
    'At': 209.990,
    'I': 126.904,
    'U': 238.029,
    'Y': 88.906,
    'Ac': 227.030,
    'Ag': 107.868,
    'Ir': 192.217,
    'Am': 243.060,
    'Al': 26.982,
    'As': 74.922,
    'Ar': 39.948,
    'Au': 196.967,
    'Zr': 91.224,
    'In': 114.818,
    'Mo': 95.940
}