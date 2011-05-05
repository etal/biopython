# Copyright (C) 2011, Joao Rodrigues (j.rodrigues@uu.nl)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""High-level wrapper and extension for Bio.PDB"""

import os, warnings

from Bio.PDB import PDBParser
try:
    from Bio.PDB import MMCIFParser
    has_mmcif = True
except ImportError:
    has_mmcif = False
    print "No module named MMCIFlex.\n" \
          "You will not be able to parse MMCIF files.\n" \
          "Please check the FAQ for details on the installation."

from Bio.PDB import PDBIO

def read(fpath, name=None, format=None, PERMISSIVE=1, quiet=False):
    """
    Wrapper for Bio.PDB parsers.
    
    Mandatory Arguments:
      fpath         - path to file

    Optional Arguments:
      name          - structure name (for Structure Object) [default: file name]
      format        - file format (pdb, mmcif) [default: PDB]
      PERMISSIVE    - PERMISSIVE flag (for PDBParser only) [default: 1]
      quiet         - Supresses warning messages while parsing [default: False]
    """
    
    # Allowed Formats
    allowed_formats = { 'pdb': PDBParser() }
    
    if has_mmcif:
        allowed_formats['mmcif'] = MMCIFParser()
    
    # Get format
    if not format:
        # Guess from extension or should it default to PDB?
        ext = fpath.split('.')[-1]
        format = ext

    format= format.lower()
    if format in allowed_formats:
        parser = allowed_formats[format]
    else:
        raise KeyError("Unrecognized format: %s" %format)
    
    # Parse structure name from path if not given
    if not name:
        ext = fpath.split('.')[-1]
        name = os.path.basename(fpath)[:-len(ext)-1]
    
    # Supress warnings
    if quiet:
        warnings.simplefilter('ignore')

    # Parse structure
    structure = parser.get_structure(name, fpath)

    # Reset warnings
    if quiet:
        warnings.resetwarnings()

    return structure