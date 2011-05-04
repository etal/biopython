"""
Functions for accessing GBIF, downloading records, processing them into a class, and extracting information from the xmltree in that class.
"""

from xml.etree import ElementTree as ET

# Functions used generally (list, string functions)
import GeneralUtils
#from Bio.Geography import GeneralUtils
from GeneralUtils import fix_ASCII_line, fix_ASCII_lines

# Used by classification utility
#from Bio.Geography.GeogUtils import point_inside_polygon
from GeogUtils import point_inside_polygon


# Used by self._open
import urllib, time, warnings
import os.path
from Bio import File



class GbifObservationRecord(Exception): pass

class GbifObservationRecord():
	"""GbifObservationRecord is a class for holding an individual observation at an individual lat/long point."""
	
	# Info prints out the informative description of this class, if a user wants to see it.
	#info = "Information about this class: GbifObservationRecord is a class for holding an individual observation at an individual lat/long point."
	
	
	def __init__(self):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		self.gbifkey = None
		self.catalognum = None
		self.country = None
		self.lat = None
		self.long = None
		self.earlydate = None
		self.taxonconceptkey = None
		self.taxonname = None
		self.namecomplete = None
		self.taxon = None
		self.genus = None
		self.species = None
		self.scientific = None
		self.latedate = None
		self.area = None # This is not taken from XML, it is determined by classify_point_into_area
		
	
	def latlong_to_obj(self, line):
		"""
		Read in a string which consists of a simple table, read species/lat/long to GbifObservationRecord object
		# This can be slow, e.g. 10 seconds for even just ~1000 records.
		"""
		print line
		words = line.split("\t")
		numwords = len(words)
		self.gbifkey = int(words[0])
		self.lat = float(words[1])
		self.long = float(words[2])
		self.taxon = str(words[3])
		self.latedate = str(words[4])
		
		temptaxon_words = self.taxon.split()
		if len(temptaxon_words) == 2:
			self.genus = str(temptaxon_words[0])
			self.species = str(temptaxon_words[1])
		if len(temptaxon_words) == 1:
			self.genus = self.taxon

		
		return


	def classify_point_into_area(self, poly, polyname):
		"""
		Fill in the self.area attribute with polygon name "poly" if it falls within the polygon "poly".  Otherwise, don't change self.area (which will be "None" or a previously-determined value. Uses GeogUtils library.
		"""
		x = self.long
		y = self.lat
		inside = point_inside_polygon(x,y,poly)
		
		if inside == True:
			self.area = polyname
		
		return self.area
		

	def parse_occurrence_element(self, element):
		"""
		Parse a TaxonOccurrence element, store in OccurrenceRecord
		"""
		
		# Get the GBIF record key for the occurrence
		try:
			self.gbifkey = element.attrib['gbifKey']
		except:
			pass
		
		# Get the catalog number
		self.catalognum = self.fill_occ_attribute(element, 'catalogNumber', 'str')
		self.country = self.fill_occ_attribute(element, 'country', 'str')
		self.lat = self.fill_occ_attribute(element, 'decimalLatitude', 'float')
		self.long = self.fill_occ_attribute(element, 'decimalLongitude', 'float')
		self.earlydate = self.fill_occ_attribute(element, 'earliestDateCollected', 'str')
		
		try:
			matching_subel = self.find_1st_matching_subelement(element, 'TaxonConcept', None)
			if matching_subel is not None:
				self.taxonconceptkey = matching_subel.attrib['gbifKey']
		except:
			pass

		self.taxonname = self.fill_occ_attribute(element, 'TaxonName', 'str')			
		self.namecomplete = self.fill_occ_attribute(element, 'nameComplete', 'str')
		self.genus = self.fill_occ_attribute(element, 'genusPart', 'str')
		self.species = self.fill_occ_attribute(element, 'specificEpithet', 'str')
		self.scientific = self.fill_occ_attribute(element, 'scientific', 'str')
		self.latedate = self.fill_occ_attribute(element, 'latestDateCollected', 'str')
		
		return
		

	def fill_occ_attribute(self, element, el_tag, format='str'):
		"""
		Return the text found in matching element matching_el.text.
		"""
		return_element = None
		matching_el = self.find_1st_matching_subelement(element, el_tag, return_element)
		
		#print matching_el
		
		# Typing these variables makes it go much faster.
		if matching_el is not None:
			result = matching_el.text
			if result == '':
				return None
			elif result == 'true':
				return bool(True)
			elif result == 'false':
				return bool(False)
			elif result == None:
				#print "None check #1"
				return None
			else:
				if format == 'str':
					return str(fix_ASCII_line(result))
				else:
					textstr = format + '(' + result + ')'
					try: 
						return eval(textstr)
					except:
						textstr = "Expected type (" + format + ") not matched. String version: " + str(fix_ASCII_line(result))
						return textstr
		else:
			return None
		return None
		

	def find_1st_matching_subelement(self, element, el_tag, return_element):
		"""
		Burrow down into the XML tree, retrieve the first element with the matching tag.
		"""
		if element.tag.endswith(el_tag):
			return_element = element
			return return_element
		else:
			children = element.getchildren()
			if len(children) > 0:
				for child in children:
					return_element = self.find_1st_matching_subelement(child, el_tag, return_element)
					
					# Check if it was found
					if return_element is not None:
						#print return_element
						return return_element
			# If it wasn't found, return the empty element
			#print return_element
			return return_element
						


	def record_to_string(self):
		"""
		Print the attributes of a record to a string
		"""
		
		temp_attributes_list = [self.gbifkey, self.catalognum, self.country, self.lat, self.long, self.earlydate, self.taxonconceptkey, self.taxonname, self.namecomplete, self.taxon, self.genus, self.species, self.scientific, self.latedate, self.area]
		
		str_attributes_list = []
		for item in temp_attributes_list:
			str_attributes_list.append(str(item))
		
		return_string = '\t'.join(str_attributes_list)
		
		return return_string
		

		

		

class GbifDarwincoreXmlString(Exception): pass

class GbifDarwincoreXmlString(str):
	"""GbifDarwincoreXmlString is a class for holding the xmlstring returned by a GBIF search, & processing it to plain text, then an xmltree (an ElementTree).
	
	GbifDarwincoreXmlString inherits string methods from str (class String).
	"""
	
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "Information about this class: GbifDarwincoreXmlString is a class for holding the xmlstring returned by a GBIF search, & processing it to plain text, then an xmltree (an ElementTree)."
	

	def __init__(self, rawstring=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		
		if rawstring:
			self.plainstring = rawstring
		else:
			self.plainstring = None
		
		#return self.plainstring
		
						
	

class GbifXmlTreeError(Exception): pass

class GbifXmlTree():
	"""gbifxml is a class for holding and processing xmltrees of GBIF records."""
	
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "Information about this class: gbifxml is a class for holding and processing xmltrees of GBIF records."

	# Used in _open function to determine 3 second waiting time for server access
	

	def __init__(self, xmltree=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		
		self.data = []
		
		if xmltree:
			self.xmltree = xmltree
			self.root = xmltree.getroot()
		else:
			print "GbifXmlTree.__init__(): No xmltree, so creating empty object."
			self.xmltree = None
			self.root = None
	
		
	
	
	def print_xmltree(self):
		"""
		Prints all the elements & subelements of the xmltree to screen (may require 
		fix_ASCII to input file to succeed)
		"""
		xmltree = self.xmltree
		
		for element in xmltree.getroot():
			print element
			self.print_subelements(element)
			
	def print_subelements(self, element):
		"""
		Takes an element from an XML tree and prints the subelements tag & text, and
		the within-tag items (key/value or whatnot)
		"""	
		if element.__len__() == 0:
			# Print the text beginning the tag header, and between the tags
			print fix_ASCII_line(element.tag), fix_ASCII_line(element.text)
	
			# Check for any key/value pairs included in the tag header,
			# and print them if they exist
			if len(element.items()) > 0:
				print "Encoded items: ", fix_ASCII_line(repr(element.items()))
			return
		elif element.__len__() > 0:
			print fix_ASCII_line(element.tag), fix_ASCII_line(element.text), "#subelements =", element.__len__()
			if len(element.items()) > 0:
				print "Encoded items: ", fix_ASCII_line(repr(element.items()))
	
			for subelement in element.getchildren():
				self.print_subelements(subelement)
			return
	
	def _element_items_to_dictionary(self, element_items):
		"""
		If the XML tree element has items encoded in the tag, e.g. key/value or
		whatever, this function puts them in a python dictionary and returns 
		them.
		"""
		
		if len(element_items) < 2:
			print "_element_items_to_dictionary error: < 2 items in element"
			return None
		else:
			temp_dict = {}
			for item in element_items:
				temp_dict[item[0][1]] = item[1][1]
			return temp_dict

		"""
		temp_dict = _element_items_to_dictionary(element.keys())
		for key in temp_dict.keys():
			print key, ": ", temp_dict[key]
			# Check if there is more than a key/value pair
			if len(element.keys()) > 2:
				print "There is more than a key/value pair encoded in this element"
				print element.items()
		"""
	
	
	def extract_latlongs(self, element):
		"""
		Create a temporary pseudofile, extract lat longs to it, 
		return results as string.
		
		Inspired by: http://www.skymind.com/~ocrow/python_string/
		(Method 5: Write to a pseudo file)
		"""
		from cStringIO import StringIO
		file_str = StringIO()
		self._extract_latlong_datum(element, file_str)
		return file_str.getvalue()



	def _extract_latlong_datum(self, element, file_str):
		"""
		Searches an element in an XML tree for lat/long information, and the 
		complete name. Searches recursively, if there are subelements.
		
		file_str is a string created by StringIO in extract_latlongs() (i.e., a temp filestr)
		"""
		if element.__len__() == 0:
			return (element.tag, element.text)
		elif element.__len__() > 0:
			#print element.tag, element.text, "#subelements =", element.__len__()
			for subelement in element.getchildren():
				# Get the TaxonOccurrence
				if subelement.tag.endswith('TaxonOccurrence'):
					for item in subelement.items():
						if item[0] == 'gbifKey':
							#print item[0], item[1]
							file_str.write(item[1] + '\t')
				
				# Get the lat/long/name
				temptuple = self._extract_latlong_datum(subelement, file_str)
				if temptuple[0].endswith('decimalLatitude'):
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('decimalLongitude'):
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('taxonName'):
					#if temptuple[1] != '':
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('latestDateCollected'):
					#if temptuple[1] != '':
					file_str.write(str(temptuple[1]) + '\n')
			return ('tag: parent subelement', 'text: multiple subelements')
			

	
	
	def extract_all_matching_elements(self, start_element, el_to_match):
		"""
		Returns a list of the elements, picking elements by TaxonOccurrence; this should 
		return a list of elements equal to the number of hits.
		"""
		
		print ''
		print 'Running extract_all_matching_elements(start_element, ' + el_to_match + ')'
		
		output_list = []
		self._recursive_el_match(start_element, el_to_match, output_list)

		return output_list
	

	def _recursive_el_match(self, element, el_to_match, output_list):
		"""
		Search recursively through xmltree, starting with element, recording all instances of el_to_match.
		"""
		# Does the element match? If so, add to list
		if element.tag.endswith(el_to_match):
			output_list.append(element)

		# If there are NO subelements, then return...
		if element.__len__() == 0:
			return output_list
		# If there are some subelements, search them...
		elif element.__len__() > 0:
			for subelement in element.getchildren():
				output_list = self._recursive_el_match(subelement, el_to_match, output_list)
			return output_list			
		# If nothing is found ever, return the below
		return ('Error: no output_list of XML elements returned')
		
	
	
	def find_to_elements_w_ancs(self, el_tag, anc_el_tag):
		"""
		Burrow into XML to get an element with tag el_tag, return only those el_tags underneath a particular parent element parent_el_tag
		"""
		xmltree = self.xmltree
		
		match_el_list = []
	
		for element in xmltree.getroot():
			
			match_el_list = self.xml_recursive_search_w_anc(element, el_tag, anc_el_tag, match_el_list)
		
		
		# How many different ancs found?
		list_ancs = []
		for tupitem in match_el_list:
			list_ancs.append(tupitem[0])
		
		unique_ancs = unique(list_ancs)
		print unique_ancs
		
		print ""
		print "Number of elements found: " + str(len(list_ancs)) + " in " + str(len(unique_ancs)) + " XML 'ancestor' (higher group) categories."
		
		for unique_anc in unique_ancs:
			print ""
			print "Anc: " + str(unique_anc)
			for tupitem in match_el_list:
				if tupitem[0] == unique_anc:
					print "	el: " + str(tupitem[1])
		
		return match_el_list
	
	
	
	
	def xml_recursive_search_w_anc(self, element, el_tag, anc_el_tag, match_el_list):
		"""
		Recursively burrows down to find whatever elements with el_tag exist inside a parent_el_tag.
		"""
		xmltree = self.xmltree
		
		# If the element matches the tag you are looking for...
		if element.tag == el_tag:
	
			# Then check if the ancestor matches
			found_anc = None
			ancestor = self._xml_burrow_up(element, anc_el_tag, found_anc)
			
			if ancestor == None:
				pass
			else:
				if ancestor.tag == anc_el_tag:
					match_el = element
					match_el_list.append((ancestor, match_el))
	
		else:
			for child in element.getchildren():
				match_el_list = self.xml_recursive_search_w_anc(child, el_tag, anc_el_tag, match_el_list)
				
		return match_el_list
			

	def create_sub_xmltree(self, element):
		"""
		Create a subset xmltree (to avoid going back to irrelevant parents)
		"""
		
		xmltree = ET.ElementTree(element)
		
		return xmltree
	
	
			
	
	def _xml_burrow_up(self, element, anc_el_tag, found_anc):
		"""
		Burrow up xml to find anc_el_tag
		"""
		xmltree = self.xmltree
		
		if found_anc == None:
			# Just get the direct parent of child_to_search_for
			child_to_search_for = element
			parent_element = self._return_parent_in_xmltree(child_to_search_for)
			
			if parent_element == None:
				return found_anc
			
			# Does the parent match the searched-for ancestor?		
			if parent_element.tag == anc_el_tag:
				found_anc = parent_element
			else:
				# Move a level up and search again, return if found
	
				found_anc = self._xml_burrow_up(parent_element, anc_el_tag, found_anc)
	
			return found_anc
			
		else:
			return found_anc
			
		
	
	def _xml_burrow_up_cousin(element, cousin_el_tag, found_cousin):
		"""
		Burrow up from element of interest, until a cousin is found with cousin_el_tag
		"""
		xmltree = self.xmltree
		
		if found_cousin == None:
			# Just get the direct parent of child_to_search_for
			child_to_search_for = element
			parent_element = self._return_parent_in_xmltree(xmltree, child_to_search_for)
			
			if parent_element == None:
				return found_cousin
			
			grandparent_element = self._return_parent_in_xmltree(parent_element)
			if grandparent_element == None:
				return found_cousin
			
			# Does the parent or any cousins match the searched-for ancestor?		
			for aunt in grandparent_element.getchildren():
				if aunt.tag == cousin_el_tag:
					found_cousin = aunt
					return found_cousin
			
			if found_cousin == None:
				# Move a level up and search again, return if found
				found_cousin = self._xml_burrow_up_cousin(parent_element, cousin_el_tag, found_cousin)
	
			return found_cousin
			
		else:
			return found_cousin
		
	
	
	def _return_parent_in_xmltree(self, child_to_search_for):
		"""
		Search through an xmltree to get the parent of child_to_search_for
		"""
		xmltree = self.xmltree
		
		returned_parent = None
		for element in xmltree.getroot():
			
			potential_parent = element
			
			returned_parent = self._return_parent_in_element(potential_parent, child_to_search_for, returned_parent)
			
			return returned_parent
			
					
		
	def _return_parent_in_element(self, potential_parent, child_to_search_for, returned_parent):
		"""
		Search through an XML element to return parent of child_to_search_for
		"""
		
		if returned_parent == None:
			children = potential_parent.getchildren()
			if len(children) > 0:
				for child in potential_parent.getchildren():
					if child == child_to_search_for:
						returned_parent = potential_parent
		
					# If not found at this level, go down a level
					else:
						returned_parent = self._return_parent_in_element(child, child_to_search_for, returned_parent)
							
			return returned_parent
		else:
			return returned_parent
			
	
	def find_1st_matching_element(self, element, el_tag, return_element):
		"""
		Burrow down into the XML tree, retrieve the first element with the matching tag
		"""
		if return_element == None:
			if element.tag == el_tag:
				return_element = element
			else:
				children = element.getchildren()
				if len(children) > 0:
					for child in children:
						return_element = self.find_1st_matching_element(child, el_tag, return_element)
						if return_element is not None:
							return return_element
			return return_element	
		else:
			return return_element
		
	

	
	def extract_numhits(self, element):
		"""
		# Search an element of a parsed XML string and find the 
		# number of hits, if it exists.  Recursively searches, 
		# if there are subelements.
		# 
		"""
		# print "Running extract_numhits(element)..."
		if element.__len__() == 0:
			if len(element.items()) > 0:
				for item in element.items():
					for index, tupleitem in enumerate(item):
						if tupleitem == 'totalMatched':
							#print item
							#print int(item[index+1])
							return int(item[index+1])
						else:
							temp_return_item = None
		elif element.__len__() > 0:
			#print element.tag, element.text, "#subelements =", element.__len__()
			for subelement in element.getchildren():
				temp_return_item = self.extract_numhits(subelement)
			if temp_return_item != None:
				return_item = temp_return_item
				return return_item
			else:
				return return_item
				
	

	
		
	
	
	












	
class GbifSearchResults(Exception): pass

class GbifSearchResults():	
	"""
	GbifSearchResults is a class for holding a series of GbifObservationRecord records, and processing them e.g. into classified areas.
	"""
	
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "GbifSearchResults is a class for holding a series of GbifObservationRecord records, and processing them	e.g. into classified areas."
	
	
	def __init__(self, gbif_recs_xmltree=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""

		# These are used by the _open function, e.g. to set wait times for 
		# accessing GBIF
		self.time_last_gbif_access = 0
		self.email = None

		# xmltree of GBIF records (GbifXmlTree class)	
		if gbif_recs_xmltree:
			self.gbif_recs_xmltree = gbif_recs_xmltree
		else:
			self.gbif_recs_xmltree = GbifXmlTree()
		
		# list of xmltrees returned by get_all_records_by_increment
		self.gbif_xmltree_list = []
		
		# xmltree of GBIF count query
		self.gbif_count_xmltree = None
		
		# String (tempfile-ish) list of records
		# An GbifDarwincoreXmlString object, with methods for cleaning the string
		# This could a list if multiple downloads are required
		self.obs_recs_xmlstring = None
		self.count_recs_xmlstring = None
		
		# Summary statistics for search
		self.gbif_hits_count = None
		self.gbif_count_params = None
		self.count_params_str = None
		
		# List of GbifObservationRecord objects
		self.obs_recs_list = []



	def print_records(self):
		"""
		Print all records in tab-delimited format to screen.
		"""
		
		for index, record in enumerate(self.obs_recs_list):
			# Get string for record
			print_string = record.record_to_string()
			
			# Fix ASCII
			print_string2 = fix_ASCII_line(print_string)
			
			# Print it
			print str(index+1) + '\t' + print_string
		
		return


	def print_records_to_file(self, fn):
		"""
		Print the attributes of a record to a file with filename fn
		"""
		
		# Open the file
		fh = open(fn, 'w')
		
		for index, record in enumerate(self.obs_recs_list):
			print_string = record.record_to_string()
		
			# Write to the file
			fh.write(str(index+1) + '\t' + print_string + '\n')

		# Close the file
		fh.close()
		
		return fn


	def classify_records_into_area(self, poly, polyname):
		"""
		Take all of the records in the GbifSearchResults object, fill in their area attribute if they fall within the polygon poly.
		"""
		
		if self.obs_recs_list == []:
			print 'Error: No records stored in self.obs_recs_list.'
			return
		
		matching_count = 0
		for record in self.obs_recs_list:
			record.classify_point_into_area(poly, polyname)
			if record.area == polyname:
				matching_count = matching_count + 1
		
		print str(matching_count) + " of " + str(len(self.obs_recs_list)) + ' fell inside area "' + polyname + '".'
		
		return


	def latlongs_to_obj(self):
		"""
		Takes the string from extract_latlongs, puts each line into a
		GbifObservationRecord object.
		
		Return a list of the objects
		"""
		
		if self.obs_recs_xmlstring == None:
			# If there is no string of records yet, extract the string,
			# using the gbif_recs_xmltree root (getroot from xmltree) as the input
			# element
			print "Filling in self.obs_recs_xmlstring"
			self.obs_recs_xmlstring = self.gbif_recs_xmltree.extract_latlongs(self.gbif_recs_xmltree.root)
			#print type(self.obs_recs_xmlstring)
			# Print the header
			for index, line in enumerate(self.obs_recs_xmlstring.splitlines()):
				print "Line" + str(index) + ": " + line
				if index > 4:
					print '...'
					break
		
		for line in self.obs_recs_xmlstring.splitlines():
			temprec = GbifObservationRecord()
			temprec.latlong_to_obj(line)
			self.obs_recs_list.append(temprec)
			#del temprec
		
		return
	

	# Functions devoted to accessing/downloading GBIF records
		
	def access_gbif(self, url, params):
		"""
		# Helper function to access various GBIF services
		# 
		# choose the URL ("url") from here:
		# http://data.gbif.org/ws/rest/occurrence
		#
		# params are a dictionary of key/value pairs
		#
		# "self._open" is from Bio.Entrez.self._open, online here: 
		# http://www.biopython.org/DIST/docs/api/Bio.Entrez-pysrc.html#self._open
		#
		# Get the handle of results
		# (looks like e.g.: <addinfourl at 75575128 whose fp = <socket._fileobject object at 0x48117f0>> )
		
		# (open with results_handle.read() )
		"""
		print 'Accessing GBIF with access_gbif...'
		
	
		results_handle = self._open(url, params)
		return results_handle
	
	
	def _get_hits(self, params):
		"""
		Get the actual hits that are be returned by a given search
		(this allows parsing & gradual downloading of searches larger 
		than e.g. 1000 records)
	
		It will return the LAST non-none instance (in a standard search result there
		should be only one, anyway).
		"""
		print '  Running _get_hits(params)...'
	
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/list'
	
		cmd = url + self._paramsdict_to_string(params)
		results_handle = self.access_gbif(url, params)
		
		self.obs_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())
				
		print 'XML search results stored via: "self.obs_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())"'
		
		return self.obs_recs_xmlstring

	
	def get_xml_hits(self, params):
		"""
		Returns hits like _get_hits, but returns a parsed XML tree.
		"""
		
		print ''
		print 'Running get_xml_hits(params)...'
	
		# Returns GbifDarwincoreXmlString object
		# (stored in self.obs_recs_xmlstring)
		self._get_hits(params)
		
		# Fix the xmlstring:
		# returns a plain string
		plain_xmlstring_fixed = self.obs_recs_xmlstring.plainstring
		#plain_xmlstring_fixed = fix_ASCII_lines(self.obs_recs_xmlstring)


		# (temp_xmlstring.splitlines() works because GbifDarwincoreXmlString inherits from string class)
		
		# Store the plain string in an GbifDarwincoreXmlString object, store as method of
		# GbifXmlTree object
		self.obs_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)
				
		# make sure the self.obs_recs_xmlstring (an GbifDarwincoreXmlString object) parses into an ElementTree
		try:
			self.gbif_recs_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(self.obs_recs_xmlstring)))
			print_it = 0
			if print_it == 1:
				print ''
				print "Printing xmlstring_fixed..."
				fixed_string = fix_ASCII_lines(self.obs_recs_xmlstring)
				print fixed_string
				print ''
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('"self.obs_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)"', inst)
		
		
		return 



	def get_record(self, key):
		"""
		Given the key, get a single record, return xmltree for it.
		"""
		
		print ''
		print 'Running get_record(params)...'
	
		# URL for the record utility
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/get'
	
		params = {'format': 'darwin', 'key' : key}
	
		cmd = url + self._paramsdict_to_string(params)
		results_handle = self.access_gbif(url, params)
		
		#print results_handle.read()
		xmlstring = GbifDarwincoreXmlString(results_handle.read())
		
		# Save to a tempfile
		"""
		fn_unfixed ='tempxml_unfixed2.xml'
		fh = open(fn_unfixed, 'w')
		fh.write(xmlstring)
		fh.close()
		"""
		
		# returns plain string, with linebreaks for parsing with ET.fromstring
		xmlstring2 = fix_ASCII_lines(xmlstring.plainstring)
		
		"""
		fh = open(fn, 'r')
		xmlstring2 = fh.read()
		fh.close()
		"""
		
		#print ''
		#print 'xmlstring2'
		#print xmlstring2
		
		xmltree = ET.ElementTree(ET.fromstring(xmlstring))
		temp = GbifXmlTree(xmltree)
		#temp.print_xmltree()
		
		return xmltree


	
	
	def get_numhits(self, params):
		"""
		Get the number of hits that will be returned by a given search
		(this allows parsing & gradual downloading of searches larger 
		than e.g. 1000 records)
	
		It will return the LAST non-none instance (in a standard search result there
		should be only one, anyway).
		"""
		self.gbif_count_params = params

		print ''
		print 'Running get_numhits(params)...'
	
		# URL for the count utility
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/count'
		
		self.count_params_str = self._paramsdict_to_string(params)
		cmd = url + self.count_params_str
		results_handle = self.access_gbif(url, params)
		
		self.count_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())
		
		
		
		# Fix xmlstring
		#plain_xmlstring_fixed = fix_ASCII_lines(self.count_recs_xmlstring)
		plain_xmlstring_fixed = self.count_recs_xmlstring
		#temp_xmlstring.splitlines() works because GbifDarwincoreXmlString inherits from string class
		
		# Store the plain string in an GbifDarwincoreXmlString object, store as method of
		# GbifXmlTree object
		self.count_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)
				
		# make sure the self.obs_recs_xmlstring (an GbifDarwincoreXmlString object) parses into an ElementTree
		try:
			self.gbif_count_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(self.count_recs_xmlstring)))
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('"self.count_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)"', inst)

		
		#print_xmltree(xmltree)

		# Get the element from the xmltree attribute of the GbifXmlTree object stored in self
		for element in self.gbif_count_xmltree.xmltree.getroot():
			temp_numhits = self.gbif_count_xmltree.extract_numhits(element)
			if temp_numhits != None:
				self.gbif_hits_count = int(temp_numhits)
		
		
		print '# hits on "' + self.count_params_str + '" = ' + str(self.gbif_hits_count)
		
		return self.gbif_hits_count
	
	

	
	def xmlstring_to_xmltree(self, xmlstring):
		"""
		Take the text string returned by GBIF and parse to an XML tree using ElementTree.  
		Requires the intermediate step of saving to a temporary file (required to make
		ElementTree.parse work, apparently)
		"""
		
		"""
		tempfn = 'tempxml.xml'
		fh = open(tempfn, 'w')
		fh.write(xmlstring)
		fh.close()
		"""
		
		# instructions for ElementTree:
		# http://docs.python.org/library/xml.etree.elementtree.html
	
		# make sure the file is findable
		try:
			xmltree = ET.fromstring(xmlstring)
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('file_str', inst)
	
		# make sure the file is parsable
		try:
			xmltree.getroot()
		except Exception, inst:
			print "Unexpected error running getroot() on text in file %s: %s" % (tempfn, inst)
	
		return xmltree
	
	




	
	def get_all_records_by_increment(self, params, inc):
		"""
		Download all of the records in stages, store in list of elements.
		Increments of e.g. 100 to not overload server
		"""

		print ''
		print "Running get_all_records_by_increment(params, inc)"

		# Set up the list of chunks to get
		numhits = self.get_numhits(params)
		print "#hits = ", numhits
		list_of_chunks = range(0, numhits-1, inc)

		# Set up list of xmlstring results
		"""
		from cStringIO import StringIO
		file_str = StringIO()
		"""
		self.gbif_xmltree_list = []
		
		
		# Set up url to parse
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/list'

		# download them by increment	
		for index, startindex in enumerate(list_of_chunks):
			if startindex + inc  > numhits:
				print "Downloading records #" + str(startindex) + "-" + str(numhits-1) + " of " + str(numhits) +" (-1)."
				params['startindex'] = str(startindex)
				params['maxresults'] = str(numhits - startindex)
			else:
				print "Downloading records #" + str(startindex) + "-" + str(startindex+inc-1) + " of " + str(numhits) +" (-1)."
				params['startindex'] = str(startindex)
				params['maxresults'] = str(inc)
	
			print params
			results_handle = self.access_gbif(url, params)
		
			# returns GbifDarwincoreXmlString
			temp_xmlstring = GbifDarwincoreXmlString(results_handle.read())
			
			# returns plain string
			#xmlstring2 = temp_xmlstring.fix_ASCII_lines('\n')
			xmlstring2 = temp_xmlstring
			


			#Bug-checking; make True if you want to see the xmlstring as a text file
			xmlstring_tofile = False
			if xmlstring_tofile == True:
				error_fh = open('error2.txt', 'w')
				error_fh.write(xmlstring2)
				error_fh.close()
				
			
			"""
			import re
			from GbifXml import GbifDarwincoreXmlString
		
			url = 'http://data.gbif.org/ws/rest/occurrence/list'
			results_handle = recs3.access_gbif(url, params)
			temp_xmlstring = GbifDarwincoreXmlString(results_handle.read())
			xmlstring2 = temp_xmlstring.fix_ASCII_lines('\n')
			
			error_fh = open('error.txt', 'w')
			error_fh.write(xmlstring2)
			error_fh.close()
			
			error_fh = open('error.txt', 'r')
			lines = error_fh.readlines()
			
			
			lines = error_fh.readlines()
			temp_xmlstring = GbifDarwincoreXmlString(lines)
			xmlstring2 = temp_xmlstring.fix_ASCII_lines('\n')
			"""
			
			# Close results_handle
			results_handle.close()

			#try:				
			# Add these latlongs to filestr
			gbif_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(xmlstring2)))
			
			# Append the xmltree
			self.gbif_xmltree_list.append(gbif_xmltree)
			
			# Search through it to get the occurrences to add
			self.extract_occurrences_from_gbif_xmltree(gbif_xmltree)
			
			"""
			except ExpatError:
			error_fh = open('error.txt', 'w')
			error_fh.write(xmlstring2)
			error_fh.close()
			print ''
			print '...failure to parse xmlstring2 into ElementTree, xmlstring2 written to error.txt...'
			#print ExpatError
			"""
		return self.gbif_xmltree_list


	
	def extract_occurrences_from_gbif_xmltree(self, gbif_xmltree):
		"""
		Extract all of the 'TaxonOccurrence' elements to a list, store them in a GbifObservationRecord.
		"""
		
		
		"""
		# Debugging:
		gbif_xmltree = recs3.gbif_xmltree_list[0]
		start_element = gbif_xmltree.xmltree.getroot()
		el_to_match = 'TaxonOccurrence'
		occurrences_list = gbif_xmltree.extract_all_matching_elements(start_element, el_to_match)
		element = occurrences_list[len(occurrences_list)-1]
		element.getchildren()
		from GbifXml import GbifObservationRecord
		self = GbifObservationRecord()
		el_tag = 'nameComplete'
		self.taxon = self.fill_occ_attribute(element, el_tag, 'str')
		return_element = None
		x=self.find_1st_matching_subelement(element, el_tag, return_element)
		"""
		
		start_element = gbif_xmltree.xmltree.getroot()
		
		# Element to pull out:
		el_to_match = 'TaxonOccurrence'


		# Get all occurrences:
		occurrences_list = gbif_xmltree.extract_all_matching_elements(start_element, el_to_match)
		
		#print occurrences_list
		
		# For each one, extract info
		for element in occurrences_list:
			# Make a temporary observation record
			temp_observation = GbifObservationRecord()
			
			# Populate it from the element
			temp_observation.parse_occurrence_element(element)
			
			#print temp_observation.record_to_string()
			
			# Add the observation to the list
			self.obs_recs_list.append(temp_observation)
		
		return


	def _paramsdict_to_string(self, params):
		"""
		# Converts the python dictionary of search parameters into a text 
		# string for submission to GBIF
		"""
		temp_outstring_list = []
		for key in params.keys():
			temp_outstring_list.append(str(key) + '=' + str(params[key]))
		outstring = '&'.join(temp_outstring_list)
		return outstring

	

	
	def _open(self, cgi, params={}):
		"""
		Function for accessing online databases.
		
		Modified from: 
		http://www.biopython.org/DIST/docs/api/Bio.Entrez-module.html

		Helper function to build the URL and open a handle to it (PRIVATE).
	
		Open a handle to GBIF.  cgi is the URL for the cgi script to access.
		params is a dictionary with the options to pass to it.  Does some
		simple error checking, and will raise an IOError if it encounters one.
	
		This function also enforces the "three second rule" to avoid abusing
		the GBIF servers (modified after NCBI requirement).
		"""
		# NCBI requirement: At least three seconds between queries
		delay = 3.0
		current = time.time()

		wait = self.time_last_gbif_access + delay - current
		if wait > 0:
			time.sleep(wait)
			self.time_last_gbif_access = current + wait
		else:
			self.time_last_gbif_access = current
		# Remove None values from the parameters
		for key, value in params.items():
			if value is None:
				del params[key]
		
		# Eliminating this bit; irrelevant in GBIF
		# Tell Entrez that we are using Biopython
		"""
		if not "tool" in params:
			params["tool"] = "biopython"
		"""
		# Tell Entrez who we are
		if not "email" in params:
			if self.email != None:
				params["email"] = email
		# Open a handle to Entrez.
		options = urllib.urlencode(params, doseq=True)
		cgi += "?" + options
		handle = urllib.urlopen(cgi)
	
		# Wrap the handle inside an UndoHandle.
		uhandle = File.UndoHandle(handle)
	
		# Check for errors in the first 5 lines.
		# This is kind of ugly.
		lines = []
		for i in range(5):
			lines.append(uhandle.readline())
		for i in range(4, -1, -1):
			uhandle.saveline(lines[i])
		data = ''.join(lines)
					   
		if "500 Proxy Error" in data:
			# Sometimes Entrez returns a Proxy Error instead of results
			raise IOError("500 Proxy Error (NCBI busy?)")
		elif "502 Proxy Error" in data:
			raise IOError("502 Proxy Error (NCBI busy?)")
		elif "WWW Error 500 Diagnostic" in data:
			raise IOError("WWW Error 500 Diagnostic (NCBI busy?)")
		elif data.startswith("Error:") :
			#e.g. 'Error: Your session has expired. Please repeat your search.\n'
			raise IOError(data.strip())
		elif data.startswith("The resource is temporarily unavailable") :
			#This can occur with an invalid query_key
			#Perhaps this should be a ValueError?
			raise IOError("The resource is temporarily unavailable")
		elif data.startswith("download dataset is empty") :
			#This can occur when omit the identifier, or the WebEnv and query_key
			#Perhaps this should be a ValueError?
			raise IOError("download dataset is empty")
		elif data[:5] == "ERROR":
			# XXX Possible bug here, because I don't know whether this really
			# occurs on the first line.  I need to check this!
			raise IOError("ERROR, possibly because id not available?")
		# Should I check for 404?  timeout?  etc?
		return uhandle
	
	




