# Copyright (C) 2011, Joao Rodrigues (j.rodrigues@uu.nl)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""High-level wrapper and extension for Bio.PDB"""

import os

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBIO

def read(fpath, format=None, PERMISSIVE=1):
    """
    Wrapper for Bio.PDB parsers.
    
    Arguments:
      fpath         - path to file
      format        - file format (pdb, mmcif)
      PERMISSIVE    - PERMISSIVE flag (for PDBParser only)
    """
    
    # Get format
    if not format or format.lower() == 'pdb':
        parser = PDBParser(PERMISSIVE=PERMISSIVE)
    elif format.lower() == 'mmcif':
        parser = MMCIFParser()
    else:
        raise KeyError("Unrecognized format: %s" %format)
    
    # Parse structure name from path
    name = os.path.basename(fpath).split('.')[-1]
    
    # Parse structure
    structure = parser.get_structure(name, fpath)
    
    return structure