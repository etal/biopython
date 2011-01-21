# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

"""Chain class, used in Structure objects."""

from Bio.PDB.Entity import Entity
import warnings

class Chain(Entity):
    def __init__(self, id):
        self.level="C"
        Entity.__init__(self, id)

    # Private methods

    def _sort(self, r1, r2):
        """Sort function for residues in a chain

        Residues are first sorted according to their hetatm records.
        Protein and nucleic acid residues first, hetatm residues next, 
        and waters last. Within each group, the residues are sorted according
        to their resseq's (sequence identifiers). Finally, residues with the
        same resseq's are sorted according to icode.

        Arguments:
        o r1, r2 - Residue objects
        """
        hetflag1, resseq1, icode1=r1.id
        hetflag2, resseq2, icode2=r2.id
        if hetflag1!=hetflag2:
            return cmp(hetflag1[0], hetflag2[0])
        elif resseq1!=resseq2:
            return cmp(resseq1, resseq2)
        return cmp(icode1, icode2)

    def _translate_id(self, id):
        """
        A residue id is normally a tuple (hetero flag, sequence identifier, 
        insertion code). Since for most residues the hetero flag and the 
        insertion code are blank (i.e. " "), you can just use the sequence 
        identifier to index a residue in a chain. The _translate_id method
        translates the sequence identifier to the (" ", sequence identifier,
        " ") tuple. 

        Arguments:
        o id - int, residue resseq 
        """
        if isinstance(id, int):
            id=(' ', id, ' ')
        return id
            
    # Special methods   

    def __getitem__(self, id):
        """Return the residue with given id.

        The id of a residue is (hetero flag, sequence identifier, insertion code). 
        If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method. 

        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.__getitem__(self, id)

    def __delitem__(self, id):
        """
        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.__delitem__(self, id)

    def __repr__(self):
        return "<Chain id=%s>" % self.get_id()

    # Public methods

    def get_unpacked_list(self):
        """Return a list of undisordered residues.

        Some Residue objects hide several disordered residues
        (DisorderedResidue objects). This method unpacks them, 
        ie. it returns a list of simple Residue objects.
        """
        unpacked_list=[]
        for residue in self.get_list():
            if residue.is_disordered()==2:
                for dresidue in residue.disordered_get_list():
                    unpacked_list.append(dresidue)
            else:
                unpacked_list.append(residue)
        return unpacked_list

    def has_id(self, id):
        """Return 1 if a residue with given id is present.

        The id of a residue is (hetero flag, sequence identifier, insertion code).
                If id is an int, it is translated to (" ", id, " ") by the _translate_id
        method.

        Arguments:
        o id - (string, int, string) or int
        """
        id=self._translate_id(id)
        return Entity.has_id(self, id)


    # Public

    def get_atoms(self):
        for r in self:
            for a in r:
                yield a

    def renumber_residues(self, res_seed=1, het_seed=0):
        """ Renumbers the residues in a chain. Keeps gaps in the chain if they exist.
            
            Uses SEQRES records to correctly filter HETATMs. In its absence,
            checks residues for alpha carbons (modified residues have it)
            and tries to infer which should be considered as ATOMs.
                        
            Options:
            
            - res_seed [int] (default:1)
            First residue number
            - het_seed [int] (default:0)
            First HETATM number. Default means first HETATM will be numbered 1000.

        """
        
        # Derive list of regular residues from SEQRES
        # Safe way to catch MSE and alike.
        h = self.parent.parent.header
        if 'SEQRES' in h and self.id in h['SEQRES']:
            seqres = h['SEQRES'][self.id]
            magic = False
        else:
            warnings.warn("WARNING: SEQRES field could not be retrieved for chain %s\n"
                          "HETATM may be oddly renumbered."
                          %self.id)
            magic = True
        
        # Calculate displacement from 1st residue
        residue_list = self.child_list
        fresnum = residue_list[0].id[1]
        displace = res_seed - fresnum
        
        # To keep track of last residue number
        # Used when renumbering consecutive chains sequentially 
        last_num = 0
        
        for residue in residue_list:
            
            if residue.id[0] == ' ': # ATOM
                last_num = residue.id[1]+displace
                residue.id = (residue.id[0], residue.id[1]+displace, residue.id[2])
            else: # HETATMs
                if not magic and residue.resname in seqres:
                    last_num = residue.id[1]+displace
                    residue.id = (residue.id[0], residue.id[1]+displace, residue.id[2])
                # Filtering HETATM for modified residues
                # Assuming all have a CA
                elif magic and 'CA' in residue.child_dict.keys() and residue['CA'].fullname == ' CA ':
                    last_num = residue.id[1]+displace
                    residue.id = (residue.id[0], residue.id[1]+displace, residue.id[2])
                else:
                    if het_seed == 0:
                        het_seed = 1000*(int(last_num/1000)+1) if last_num/1000 > 1 else 1000
                    residue.id = (residue.id[0], het_seed, residue.id[2])
                    het_seed += 1

        return (last_num, het_seed)
            