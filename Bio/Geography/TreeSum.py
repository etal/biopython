"""
TreeSum.py: Objects and functions concerned with reading in and manipulating phylogenetic trees.
"""



# Used generally
#import handyfunctions
#from Bio.Geography.GeneralUtils import make_NaN_array, list1_items_in_list2, list1_items_not_in_list2

from GeneralUtils import make_NaN_array, list1_items_in_list2, list1_items_not_in_list2


import copy
import random
from numpy import average, nanmin, isnan, isfinite, std
from pylab import normpdf
from Bio.Nexus.Trees import Tree

class TreeSumError(Exception): pass


class TreeSum:
	"""
	TreeSum is a class for summarizing and manipulating tree objects.  Attempting to use Bio.Nexus.Trees with Brad Chapman's fix for ultrametric trees.
	"""

	# Info prints out the informative description of this class, if a user wants to see it.
	info = "Information about this class: TreeSum is a class for summarizing and manipulating tree objects."
	
	def __init__(self, treeobj=None, parent_tree=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		
		self.data = []
		self.ttlPD = 0.0 	#total phylogenetic distance
		self.half_distance_matrix = None   # distance matrix on one side of diagonal
		
		# Mean Phylogenetic Distance
		self.mpd = None
		self.mean_null_mpd = None
		self.sd_null_mpd = None
		self.numnulls = None
		self.nri = None
		self.mpd_pval = None
		
		# Mean Nearest-Neighbor Phylogenetic Distance
		self.mnpd = None
		self.mean_null_mnpd = None
		self.sd_null_mnpd = None
		self.num_mnpd_nulls = None
		self.nti = None
		self.mnpd_pval = None
		
		if treeobj:
			self.treeobj = treeobj
			self.newickstr = self.get_newickstr()
		else:
			print "TreeSum.__init__(): No treeobj, so creating empty object."
			self.treeobj = None
		
		if parent_tree:
			self.parent_tree = parent_tree
		else:
			self.parent_tree = treeobj
		
	def hello(self):
		return 'hello world'
	

	def get_newickstr(self):
		"""
		Returns a string with the tree in Newick format, with branch lengths
		"""
		newickstr = self.treeobj.to_string(support_as_branchlengths=False,branchlengths_only=True,plain=False,plain_newick=True,ladderize=None)
		return newickstr
		
	
	def print_tree_outline_format(self):
		"""
		Prints the tree out in "outline" format (daughter clades are indented, etc.)
		"""
		phylo_obj = self.treeobj
		
		print ""
		print "Printing out tree of", str(len(phylo_obj.get_terminals())), "taxa, hierarchical indented format."
		print ""
		
		nodenum = phylo_obj.root
		node = phylo_obj.chain[nodenum]
		rank = 1
		self.print_Node(node, rank)
		
		return
	
	
	
	def print_Node(self, node, rank):
		"""
		Prints the node in question, and recursively all daughter nodes, maintaining rank as it goes.
		"""
		
		tabstr = ""
		for count in range(1,rank):
			tabstr = tabstr + "	"
		
		
		printstr = ''.join(["Rank #", str(rank), ", label=", str(node.data.taxon), ", brlen=", str(node.data.branchlength)])
		
		printstr2 = tabstr + printstr
		
		print printstr2
		
		if len(node.get_succ()) > 0:
			rank = rank + 1
			for child_nodenum in node.get_succ():
				child = self.treeobj.chain[child_nodenum]
				self.print_Node(child, rank)
	
	
		
	
	def read_ultrametric_Newick(self, newickstr):
		"""
		Read a Newick file into a tree object (a series of node objects links to parent and daughter nodes), also reading node ages and node labels if any. 
		"""
		
		# Read in the tree (Bio.Nexus.Tree, version 2)
		self.treeobj = Tree(newickstr)
		return self.treeobj
	

	def list_leaves(self):
		"""
		Print out all of the leaves in a tree
		"""
		phylo_obj = self.treeobj
		
		# Go through all the leaves and print them
		for leaf_nodenum in phylo_obj.get_terminals():
			node = phylo_obj.chain[leaf_nodenum]
			print node.data.taxon
		return len(phylo_obj.get_terminals())
	

	def treelength(self):
		"""
		Gets the total branchlength in a tree by recursively adding through tree.
		"""
		
		PD = 0.0
		# Get all node ids
		#treeobj.all_ids() = treeobj.chain.keys()
		for nodeid in self.treeobj.chain.keys():
			tempnode = self.treeobj.chain[nodeid]
			PD = PD + tempnode.data.branchlength
			#print PD
		
		self.ttlPD = PD
		return PD


	def test_Tree(self):
		"""
		Run a series of operations on the TreeSum object (self), and contained tree (self.treeobj) in question, to see if all operations work.
		"""
		print ''
		print 'Running self.test_Tree()...'
		print '  This runs a series of operations on the TreeSum object (self), and contained tree (self.treeobj) in question, to see if all operations work.'
		print ''

		print 'self.treeobj:', self.treeobj
		print 'self.newickstr:', self.newickstr
		
		print ''
		print 'Printing tree in outline format'
		self.print_tree_outline_format()
		
		
		
		ts_reading_own_newick = TreeSum()
		ts_reading_own_newick.read_ultrametric_Newick(self.newickstr)
		
		print ''
		print 'ts_reading_own_newick.list_leaves():'
		ts_reading_own_newick.list_leaves()
		
		print ''
		print 'ts_reading_own_newick.treeobj.display():'
		ts_reading_own_newick.treeobj.display()
		
		print ''
		print 'Calculating ttlPD for tree...'
		print 'self.treelength():'
		self.treelength()
		
		leafs_ids = self.treeobj.get_terminals()
		leaf1_id = leafs_ids[0]
		leaf2_id = leafs_ids[len(leafs_ids) - 1]
		print 'Common ancestor of node ' + str(leaf1_id) + ' & ' + str(leaf2_id)
		print self.treeobj.common_ancestor(leaf1_id, leaf2_id)
		
		print 'PD between node ' + str(leaf1_id) + ' & ' + str(leaf2_id)
		print self.treeobj.distance(leaf1_id, leaf2_id)
		
		
		print ''
		print 'Calculating phylodistance array'
		dist_array = self.get_half_distance_matrix()
		print dist_array
		
		print ''
		print 'Calculating MPD'
		print 'original mpd:', self.mpd
		mpd = self.mean_phylo_dist(sample_ifbig = True)
		print 'Calculated MPD:', self.mpd
		nri = self.calc_nri()
		print 'Calculating NRI:'
		print 'mean_null_mpd:', self.mean_null_mpd
		print 'sd_null_mpd:', self.sd_null_mpd
		print 'numnulls:', self.numnulls
		print 'nri:', self.nri
		print 'mpd_pval:', self.mpd_pval

		
		print ''
		print 'Calculating MNPD'
		print 'original mnpd:', self.mnpd
		mnpd = self.mean_nn_dist(sample_ifbig = True)
		print 'Calculated MNPD:', self.mnpd
		print 'Calculating NTI:'
		nti = self.calc_nti()
		print 'mean_null_mnpd:', self.mean_null_mnpd
		print 'sd_null_mnpd:', self.sd_null_mnpd
		print 'num_mnpd_nulls:', self.num_mnpd_nulls
		print 'nti:', self.nti
		print 'mnpd_pval:', self.mnpd_pval
	


		
		print ''
		print 'self.test_Tree() completed without crashes.'
		return
	
	
	def get_half_distance_matrix(self):
		"""
		Get a matrix of all of the pairwise distances between the tips of a tree.
		(just gets the distances above the diagonal; diagonal is 0.0 here  (distance of tip to itself), but recorded as None to save computation later.)
		"""
		
		phylo_obj = self.treeobj	
		leaves = phylo_obj.get_terminals()
		
		numleaves = len(leaves)
		numdists = (numleaves^2) / 2 - numleaves/2
		if numdists > 1000:
			print "Number of leaves in tree =", numleaves
			print "Number of distances to calculate =", numdists

			print "This could be slow (minutes for >100 taxa), recommend you use the phylocom package (phydist command) to calculate the matrix efficiently)."

		
		# Get the PD between leaves for each pair
		nleaves = len(leaves)
		dist_array = make_NaN_array(nleaves, nleaves)
		
		for index1 in range(0, nleaves):
			# Staring with index1 skips the diagonals (which equal 0.0)
			for index2 in range(index1+1, nleaves):
				dist_array[index1, index2] = phylo_obj.distance(leaves[index1], leaves[index2])
		
		# Set the diagonals to None (saves steps later; the diagonals are 0.0, if anyone wants them
		#dist_array = set_diags_to_none(dist_array)
		
		self.half_distance_matrix = dist_array
		return dist_array
	
	
	
	def subset_Tree(self, list_leafids_to_keep):
		"""
		Given a list of tips and a tree, remove all other tips and resulting redundant nodes to produce a new smaller tree.
		"""
		
		# Get the old list of leafids
		old_leafids_list = self.treeobj.get_terminals()
		
		# Get the items shared (between BOTH lists)
		shared_between_lists = list1_items_in_list2(list_leafids_to_keep, old_leafids_list)
		
		# Hopefully all of list1 was found in list2, but check if not
		if len(shared_between_lists) != len(list_leafids_to_keep):
			print ''
			print 'Error!  Some leafids in list_leafids_to_keep are not found in the parent tree!'
			pass
		
		#print 'old_leafids_list:', old_leafids_list
		#print 'shared_between_lists', shared_between_lists
		
		# Get the parts of old_leafids_list to drop
		list_leafids_to_drop = list1_items_not_in_list2(old_leafids_list, shared_between_lists)
		
		#print 'list_leafids_to_drop', list_leafids_to_drop
		
		# Make a new subtree
		#print 'Making a new subtree (deepcopy)...'
		newtree = copy.deepcopy(self.treeobj)
		
		for leafid in list_leafids_to_drop:
			# prune works on taxon names, so they had better be unique
			newtree.prune(newtree.chain[leafid].data.taxon)
		
		return newtree




	def get_successor_leafids(self, nodeid):
		"""Return a list of all terminal nodes above a given node id."""
		treeobj = self.treeobj
		root = treeobj.root
		
		new_successors_list = []
		
		# See which of the leafids in the original tree is a descendent of nodeid
		for leafid in treeobj.get_terminals():
			# Include root in ancs_list, since trace doesn't for some reason
			ancs_list = [root] + self.treeobj.trace(root, leafid)
			for anc in ancs_list:
				#print treeobj.chain[anc].id, nodeid
				if treeobj.chain[anc].id == nodeid:
					new_successors_list.append(leafid)
					break
				else:
					pass
		
		return new_successors_list
		
	
	def mean_phylo_dist(self, sample_ifbig = False):
		"""
		Calculates the Mean Phylogenetic Distance (MPD) for the tree.  If the tree is huge, we could just sample the distances instead of doing the entire distance matrix.
		"""
		
		if sample_ifbig == True:
			numleaves = len(self.treeobj.get_terminals())
			if numleaves * numleaves > 1000:
				# Then we're going to sample the distance matrix instead
				print "Note: Taking 100 sample pairs of distance matrix instead of calculating all " +  str(numleaves * numleaves) + " pairwise distances."
				
				# 100 samples should do it
				# Sample distances directly off the tree, not the distance matrix
				numsamps = 100
				randdist_list = self.sample_distances(numsamps)
				
				self.mpd = average(randdist_list)
				return self.mpd
		
		# Otherwise, do it the standard way
		# Generate distance matrix if it doesn't exist
		if self.half_distance_matrix == None:
			self.get_half_distance_matrix()
		
		"""
		dists_list = []
		for row in self.half_distance_matrix:
			dists_list.extend(i for i in row if i > 0.0)
		"""
		
		# Average only over non-NaN values
		mask = isnan(self.half_distance_matrix) ==  False
		self.mpd = average(self.half_distance_matrix[mask])
		return self.mpd


	def mean_nn_dist(self, sample_ifbig = False):
		"""
		Calculates the Mean Nearest Neighbor Phylogenetic Distance (MNPD) for the tree.  If the tree is huge, we could just sample the distances instead of doing the entire distance matrix.
		"""

		if sample_ifbig == True:
			numleaves = len(self.treeobj.get_terminals())
			if numleaves * numleaves > 1000:
				# Then we're going to sample the distance matrix instead
				print "Note: Taking 100 sample pairs of distance matrix instead of calculating all " +  str(numleaves * numleaves) + " pairwise distances."
				
				# 100 samples should do it
				# Sample distances directly off the tree, not the distance matrix
				numsamps = 100
				randdist_list = self.sample_min_distances(numsamps)
				
				self.mnpd = average(randdist_list)
				return self.mnpd
		
		# Otherwise, do it the standard way
		
		# Generate distance matrix if it doesn't exist
		if self.half_distance_matrix == None:
			self.get_half_distance_matrix()
		
		# axis 0 = vertical, axis 1 = horizontal
		row_min_pds_list = nanmin(self.half_distance_matrix, 0)
		
		self.mnpd = average(row_min_pds_list[isfinite(row_min_pds_list)])
		return self.mnpd




	
	def sample_distances(self, numsamps):
		"""
		If the distance matrix is huge, it is slow to calculate it repeatedly as a null model.  So, we can randomly sample numsamps pairwise distances instead.
		"""
		leaf_ids = self.treeobj.get_terminals()
		randdist_list = []
		# Random generate numsamps pairs
		for i in range(0, numsamps):
			# Sample two leaves, without replacement
			randpair = random.sample(leaf_ids, 2)
			randdist = self.treeobj.distance(randpair[0], randpair[1])
			#print randdist
			randdist_list.append(randdist)

		return randdist_list

	def sample_min_distances(self, numsamps):
		"""
		If the distance matrix is huge, it is slow to calculate it repeatedly as a null model.  So, we can randomly sample numsamps pairwise distances instead.
		"""
		leaf_ids = self.treeobj.get_terminals()
		randdist_list = []
		# Random generate numsamps pairs
		for i in range(0, numsamps):
			# Sample two leaves, without replacement
			randtip = random.sample(leaf_ids, 1)[0]
			
			# Get immediate ancestor of randtip
			parent_id = self.treeobj.chain[randtip].get_prev()
			
			# Direct sister nodes could have daughters, so use get_successor_leafids
			sisters = []
			sisters = self.get_successor_leafids(parent_id)
			#print randtip, sisters
			
			# Remove the original tip
			sisters.remove(randtip)
			
			if len(sisters) < 1:
				print "Error: tip has no immediate sisters."

			temp_minslist = []
			for sister in sisters:
				mindist = self.treeobj.distance(sister, randtip)
				temp_minslist.append(mindist)
			
						
			randdist = nanmin(temp_minslist)
			#print randdist
			randdist_list.append(randdist)

		return randdist_list

		
	def randomize_leaves(self):
		"""
		Make a new tree that randomizes leaf positions. (For generation of a null model to generate standardized Net Related Index (NRI, standardized Mean Phylogenetic Distance) and Nearest Taxon Index (NTI, standardized Mean Minimum Phylogenetic Distance).
		"""
		oldtree_nodeids = self.treeobj.get_terminals()
		oldtree_nodeids.sort()
		
		newtree = copy.deepcopy(self.treeobj)
		
		# Deep copy list to new list
		newtree_nodeids = list(oldtree_nodeids)
		# Randomly shuffle in place
		random.shuffle(newtree_nodeids)
		
		for index,leafid_to_change in enumerate(oldtree_nodeids):
			#node_to_change = leaf
			#node_to_change.set_id(newtree_nodeids[index])
			
			node_to_change = newtree.chain[leafid_to_change]
			
			# Get the data from a random leaf of the old tree:
			data_for_newnode = self.treeobj.chain[newtree_nodeids[index]].data
			
			# Store the data in the new leaf on the new tree
			newtree.chain[leafid_to_change].set_data(data_for_newnode)
			
			# But KEEP the branchlengths from the old node
			newtree.chain[leafid_to_change].data.branchlength = self.treeobj.chain[oldtree_nodeids[index]].data.branchlength
			
		return newtree
		


	def calc_nri(self, numnulls = 100):
		"""
		Calculate the Net Related Index (NRI, standardized Mean Phylogenetic Distance) based on numnulls number of randomized trees (flipping node labels; this is not the only way to do it, see phylocom package). 
		"""	
		
		# Use the larger parent tree as the null
		parent_ts = TreeSum(self.parent_tree)
		
		# Size subtree
		size_nulltree = len(self.treeobj.get_terminals())
		pool = self.parent_tree.get_terminals()
		
		null_mpds_list = [None] * numnulls
		for i in range(0, numnulls):
			#nulltree = rand_ts.randomize_leaves()
			#nullts = TreeSum(nulltree)
			list_leafids_to_keep = random.sample(pool, size_nulltree)
			
			rand_subtree = parent_ts.subset_Tree(list_leafids_to_keep)
			#print rand_subtree.get_terminals()
			rand_ts = TreeSum(rand_subtree)

			null_mpds_list[i] = rand_ts.mean_phylo_dist(sample_ifbig = True)
		
	
		self.mean_null_mpd = average(null_mpds_list)
		self.sd_null_mpd = std(null_mpds_list)
		self.numnulls = numnulls
		self.nri = -1 * (self.mpd - self.mean_null_mpd) / self.sd_null_mpd
		self.mpd_pval = normpdf(self.nri, 0, 1)

		"""
		print ''
		print "mean_null_mpd:", self.mean_null_mpd
		print "sd_null_mpd:", self.sd_null_mpd
		print "nri:", self.nri
		print "mpd_pval:", self.mpd_pval	
		"""
		

	def calc_nti(self, numnulls = 100):
		"""
		Calculate the Nearest Taxon Index (NTI, standardized Mean Nearest-Neighbor Phylogenetic Distance) based on numnulls number of randomized trees (flipping node labels; this is not the only way to do it, see phylocom package). 
		"""

		# Use the larger parent tree as the null
		parent_ts = TreeSum(self.parent_tree)
		
		# Size subtree
		size_nulltree = len(self.treeobj.get_terminals())
		pool = self.parent_tree.get_terminals()

		null_mnpds_list = [None] * numnulls
		for i in range(0, numnulls):
			#nulltree = self.randomize_leaves()
			#nullts = TreeSum(nulltree)
			list_leafids_to_keep = random.sample(pool, size_nulltree)
			
			rand_subtree = parent_ts.subset_Tree(list_leafids_to_keep)
			#print rand_subtree.get_terminals()
			rand_ts = TreeSum(rand_subtree)

			null_mnpds_list[i] = rand_ts.mean_nn_dist(sample_ifbig = True)
		
	
		self.mean_null_mnpd = average(null_mnpds_list)
		self.sd_null_mnpd = std(null_mnpds_list)
		self.num_mnpd_nulls = numnulls
		self.nti = -1 * (self.mnpd - self.mean_null_mnpd) / self.sd_null_mnpd
		self.mnpd_pval = normpdf(self.nti, 0, 1)
		
		"""
		print ''
		print "mean_null_mpd:", self.mean_null_mpd
		print "sd_null_mpd:", self.sd_null_mpd
		print "nri:", self.nri
		print "mpd_pval:", self.mpd_pval
		"""

