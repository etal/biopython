# Copyright (C) 2011, Joao Rodrigues (j.rodrigues@uu.nl)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrapper and for structural file formats.

This API follows the same semantics as Biopython's SeqIO and AlignIO.
"""
__docformat__ = "restructuredtext en"

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
    """Parse a molecular structure from the given file or handle.

    Wrapper for Bio.PDB parsers.

    Other parser-specific arguments can be passed as keyworded arguments.

    :Parameters:
        infile : string or file-like object
            Path to file, or file-like handle
        format : string (optional)
            File format. One of 'pdb' or 'mmcif'; default 'pdb'.
        id : string (optional)
            Structure Id. Typically the PDB ID, but can be anything.
    """
    # Get format and instantiate
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
    """Write a structure to a file or handle in the given format.

    Other writer-specific arguments can be passed as keyworded arguments.

    :Parameters:
        structure : Bio.PDB.Structure instance
            Structure Data (SMCRA representation)
        ofile : string or file-like object
            Path to file, or file-like handle
        format : string (optional)
            File format. Currently, only 'pdb' is supported.
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

