# Copyright (C) 2011, Joao Rodrigues (j.rodrigues@uu.nl)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrapper and for structural file formats.

This API follows the same semantics as Biopython's SeqIO and AlignIO.
"""

import os, warnings

# Write-related imports
from Bio.PDB import PDBIO

supported_o_formats = {
        'pdb': PDBIO 
        }

# Read-related imports 
from Bio.PDB import PDBParser

supported_i_formats = {
        'pdb': PDBParser 
        }

# Public Functions
def read(infile, format='pdb', id=None, **kwargs):
    """
    Wrapper for Bio.PDB parsers.
    
    Mandatory Arguments:
        file    - Path to File or File Handle

    Optional Arguments:
        format  - File Format (e.g. 'pdb')
        sid     - Structure Id to pass to Parser
    
    Other parser specific arguments can be passed as
    keyworded arguments.
    """
    
    # Get format and instanciate
    if format.lower() == 'mmcif':
        # Flex is an optional dependency
        try:
            from Bio.PDB import MMCIFParser
        except ImportError:
            message = "Missing module MMCIFlex. Please check the FAQ."
            raise MissingPythonDependencyError(message)
        else:
            supported_i_formats['mmcif'] = MMCIFParser
    
    p = supported_i_formats[format.lower()]
    Parser = p(**kwargs)
    
    if not id: # Take file name without the extension
        if hasattr(infile, 'name'):
            infile_name = infile.name
        elif isinstance(infile, basestring):
            infile_name = infile
        else: # StringIO I hope!
            infile_name = "structure.ext" # Ugly but avoids code redundancy
        ext = infile_name.split('.')[-1]
        id = os.path.basename(infile_name)[:-len(ext)-1]

    structure = getattr(Parser, 'get_structure')(id, infile)

    return structure

def write(structure, ofile, format='pdb', **kwargs):
    """
    Writes a structure to a file in a given format.
    
    Mandatory Arguments:
        structure   - Structure Data (SMCRA representation)
        ofile       - Path to File or File Handle

    Optional Arguments:
        format  - File Format (e.g. 'pdb')
    
    Other specific arguments can be passed as
    keyworded arguments.
    """

    w = supported_o_formats[format.lower()]
    Writer = w(**kwargs)
    
    if isinstance(ofile, basestring):
        of = open(ofile, 'w')
    else:
        of = ofile
    
    getattr(Writer, 'set_structure')(structure)
    n = getattr(Writer, 'save')(of)
    
    if isinstance(ofile, basestring):
        of.close()
    
    return n
    