#!/usr/bin/python

print ""
print ""
print "================================================"
print "Remote access to GBIF; classifying points into areas"
print "================================================"
print ""

from Bio.Geography.GbifXml import GbifXmlTree, GbifSearchResults
#from GbifXml import GbifXmlTree, GbifSearchResults

params = {'format': 'darwin', 'scientificname': 'Genlisea*'}
inc = 100
recs3 = GbifSearchResults()
gbif_xmltree_list = recs3.get_all_records_by_increment(params, inc)

#recs3.print_records()



# Set up polygon
# (long, lat) = x, y
ul = (-180, 90)
ur = (180, 90)
ll = (-180, 0)
lr = (180, 0)
poly = [ul, ur, ll, lr]
polyname = "NorthernHemisphere"

recs3.classify_records_into_area(poly, polyname)
recs3.print_records()

print ''
recs3.classify_records_into_area(poly, polyname)
