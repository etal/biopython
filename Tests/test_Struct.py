# Copyright 2011 by Joao Rodrigues.  All rights reserved.
# 
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.Struct module."""
import os
import tempfile
import unittest
import warnings
from StringIO import StringIO

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Struct.")

from Bio import Struct
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

class ParseTest(unittest.TestCase):
    
    def setUp(self):
        self.structure = 'PDB/a_structure.pdb'
    
    def test_read(self):
        """ Wrapper for Structural Parser (PDBParser, MMCIF)"""
        warnings.simplefilter('ignore', PDBConstructionWarning)
        s = Struct.read(self.structure)
        warnings.filters.pop()
        
        self.assertEqual('a_structure', s.id)
        self.assertEqual(167, len([r for r in s.get_residues()]))

    def test_write(self):
        """ Wrapper for Structural File Writers (PDBIO) """

        warnings.simplefilter('ignore', PDBConstructionWarning)        
        s = Struct.read(self.structure)
        warnings.filters.pop()

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            Struct.write(s, filename)
        finally:
            os.remove(filename)

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)