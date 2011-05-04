#!/usr/bin/python

print ""
print ""
print "================================================"
print "Read & parse a GBIF XML file (DarwinCore format)"
print "================================================"
print ""


import os
from Bio.Geography.GbifXml import GbifXmlTree, GbifSearchResults
from Bio.Geography.GenUtils import fix_ASCII_file
#from GbifXml import GbifXmlTree, GbifSearchResults
#from GeneralUtils import fix_ASCII_file




# Example filename
# You can find this file in the biopython directory Tests/Geography
# Example output in Tests/output/test_GbifXml
xml_fn = 'utric_search_v2.xml'

print ""
print "Reading a GBIF XML file (DarwinCore format), parsing it, and extracting lat/longs."
print "Example filename is:", xml_fn


print "Converting XML file to pure ASCII (stops crashes later during print-to-screen etc.)..."
xml_fn_new = fix_ASCII_file(xml_fn)




from xml.etree import ElementTree as ET
#from geogUtils import print_subelements, extract_latlong

# Name of file to output lats/longs to
#outfilename = 'latlongs.txt'

try:
	xmltree = ET.parse(xml_fn_new)
except Exception, inst:
	print "Unexpected error opening %s: %s" % (xml_fn, inst)


# Store results in an object of Class GbifXmlTree:
gbif_recs_xmltree = GbifXmlTree(xmltree)



# ======================
# Testing print_xmltree
# ======================
# Print the object:
#  (also uses gbif_recs_xmltree.print_sublelements)

#print ''
#print 'Printing the GBIF object xmltree with print_xmltree...'
#gbif_recs_xmltree.print_xmltree()


"""
# ======================
# Testing extract_latlongs
# ======================
# Initiate GbifSearchResults object, containing gbif_recs_xmltree
recs = GbifSearchResults(gbif_recs_xmltree)

outstr = recs.gbif_recs_xmltree.extract_latlongs(gbif_recs_xmltree.root)
print ''
print 'Printing outstr'
print outstr
#recs.latlongs_to_obj()
print ''
print 'Printing all records'
recs.print_records()
"""



# Make recs object hold all of the observation records
recs = GbifSearchResults(gbif_recs_xmltree)
recs.extract_occurrences_from_gbif_xmltree(recs.gbif_recs_xmltree)
print ''
print 'Printing all records'
recs.print_records()

print ''
print "Printing first five record objects..."
print recs.obs_recs_list[0:4], '...'


print ''
print '# ========================'
print '# Testing get_record'
print '# ========================'
key = 175067484
recs3 = GbifSearchResults()
xmlrec = recs3.get_record(key)
print 'xmltree result:'
print xmlrec

print ''
print '# ========================'
print '# Testing get_numhits()'
print '# ========================'
params = {'format': 'darwin', 'scientificname': 'Utricularia*'}
recs4 = GbifSearchResults()
recs4.get_numhits(params)



