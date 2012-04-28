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

# Check for an mmCIF parser. It may be missing, either because of an unbuilt C
# extension (standard lex/yacc) or the pure-Python fallback which depends on the
# third-party package PLY.
try:
    from Bio.PDB import MMCIFParser
except ImportError:
    # mmCIF format isn't available, but don't bother the user about it until
    # it's requested as a format in the read() function (below).
    pass
else:
    supported_i_formats['mmcif'] = MMCIFParser


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
    if format.lower() == 'mmcif' and 'mmcif' not in supported_i_formats:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError(
            "Module MMCIParser is missing because of unmet dependencies "
            "(unbuilt C extension or missing PLY package).")

    # Get format and instantiate
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

