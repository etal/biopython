"""
Some simple functions using free libraries for reading shapefiles, performing point-in-polygon operations, etc.
"""

def readshpfile(fn):
	
	# Get the filename
	# fn = shapefilename
	
	# Try to open the file, print error if error	
	try:
		f = open(fn)
	# some code here
	except IOError, reason:
		print "couldn't open file, because of", reason
		return

	# load the shapefile, populating a list of dictionaries
	shpRecords = shpUtils.loadShapefile(fn)
	
	return shpRecords
	

def summarize_shapefile(fn, output_option, outfn):
	shpRecords = readshpfile(fn)
	
	if output_option == 'toscreen':
		print "Feature#	Name	X	Y"
		
		for index, record in enumerate(shpRecords):
			
			# Make an output string
			
			# consecutive ount
			pt1 = str(index+1)
			
			# Feature name
			pt2 = str(record['dbf_data']['NAME'])
			
			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])
	
			printstr = '	'.join([pt1, pt2, pt3, pt4])
			print printstr
		
		return "printed string to screen"
			

	elif output_option == 'tofile':
		outfile = open(outfn, 'w')
		
		header_str = '	'.join(['Feature#', 'Name', 'X', 'Y'])
		outfile.write(header_str)

		
		for index, record in enumerate(shpRecords):
			
			# Make an output string
			
			# consecutive ount
			pt1 = str(index+1)
			
			# Feature name
			pt2 = str(record['dbf_data']['NAME'])
			
			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])
	
			printstr = '	'.join([pt1, pt2, pt3, pt4])
			outfile.write(printstr)
		
		outfile.close()
		return outfn


	elif output_option == 'tolist':
		
		# Create blank list
		outlist = []
		
		# Determine header row
		header_row = ['Feature#', 'Name', 'X', 'Y']
		outlist.append(header_row)
		
		for index, record in enumerate(shpRecords):
			
			# Make an output string
			
			# consecutive ount
			pt1 = str(index+1)
			
			# Feature name
			pt2 = str(record['dbf_data']['NAME'])
			
			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])
	
			line = [pt1, pt2, pt3, pt4]
			outlist.append(line)
		
		return outlist

	return



def point_inside_polygon(x,y,poly):
	"""
	# Code from here:
	# http://www.ariel.com.au/a/python-point-int-poly.html
	#
	# NOTE: does not presently deal with the case of polygons that
	# cross the international dateline!
	#
	# determine if a point is inside a given polygon or not
	#
	# Polygon is a list of (x,y) pairs.
	"""
	n = len(poly)
	inside = False

	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y

	return inside



def shapefile_points_in_poly(pt_records, poly):
	num_inside = 0
	for index,pt in enumerate(pt_records):
		x = pt['shp_data']['x']
		y = pt['shp_data']['y']
		inside = point_inside_polygon(x, y, poly)
		
		if inside == True:
			num_inside = num_inside + 1
		print index +1, inside
	
	print num_inside, " out of ", len(pt_records), " inside specified polygon."
	return



def tablefile_points_in_poly(fh, ycol, xcol, namecol, poly):
	num_inside = 0
	for index,line in enumerate(fh):
		words = line.split()
		
		if len(words) < 3:
			continue
		
		x = float(words[xcol])
		y = float(words[ycol])
		inside = point_inside_polygon(x, y, poly)
		
		if inside == True:
			num_inside = num_inside + 1
		
		if len(words) > 3:
			name = words[namecol] + ' ' + words[namecol+1]
		else:
			name = words[namecol]
		print index +1, name, inside
	
	print num_inside, " out of ", index+1, " inside specified polygon."
	return









# DBF utilities (dbfs are tables holding attribute data for shapefiles)
# Source: http://code.activestate.com/recipes/362715/
# dbfUtils.py
# By Raymond Hettinger
# Also used in Google Summer of Code here:
# http://code.google.com/p/primary-maps-2008/source/browse/trunk/dbfUtils.py

import struct, datetime, decimal, itertools

def dbfreader(f):
    """Returns an iterator over records in a Xbase DBF file.

    The first row returned contains the field names.
    The second row contains field specs: (type, size, decimal places).
    Subsequent rows contain the data records.
    If a record is marked as deleted, it is skipped.

    File should be opened for binary reads.

    """
    # See DBF format spec at:
    #     http://www.pgts.com.au/download/public/xbase.htm#DBF_STRUCT

    numrec, lenheader = struct.unpack('<xxxxLH22x', f.read(32))    
    numfields = (lenheader - 33) // 32

    fields = []
    for fieldno in xrange(numfields):
        name, typ, size, deci = struct.unpack('<11sc4xBB14x', f.read(32))
        name = name.replace('\0', '')       # eliminate NULs from string   
        fields.append((name, typ, size, deci))
    yield [field[0] for field in fields]
    yield [tuple(field[1:]) for field in fields]

    terminator = f.read(1)
    assert terminator == '\r'

    fields.insert(0, ('DeletionFlag', 'C', 1, 0))
    fmt = ''.join(['%ds' % fieldinfo[2] for fieldinfo in fields])
    fmtsiz = struct.calcsize(fmt)
    for i in xrange(numrec):
        record = struct.unpack(fmt, f.read(fmtsiz))
        if record[0] != ' ':
            continue                        # deleted record
        result = []
        for (name, typ, size, deci), value in itertools.izip(fields, record):
            if name == 'DeletionFlag':
                continue
            if typ == "N":
                value = value.replace('\0', '').lstrip()
                if value == '':
                    value = 0
                elif deci:
                    value = decimal.Decimal(value)
                else:
                    value = int(value)
            elif typ == 'D':
                y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                value = datetime.date(y, m, d)
            elif typ == 'L':
                value = (value in 'YyTt' and 'T') or (value in 'NnFf' and 'F') or '?'
            result.append(value)
        yield result


def dbfwriter(f, fieldnames, fieldspecs, records):
    """ Return a string suitable for writing directly to a binary dbf file.

    File f should be open for writing in a binary mode.

    Fieldnames should be no longer than ten characters and not include \x00.
    Fieldspecs are in the form (type, size, deci) where
        type is one of:
            C for ascii character data
            M for ascii character memo data (real memo fields not supported)
            D for datetime objects
            N for ints or decimal objects
            L for logical values 'T', 'F', or '?'
        size is the field width
        deci is the number of decimal places in the provided decimal object
    Records can be an iterable over the records (sequences of field values).
    
    """
    # header info
    ver = 3
    now = datetime.datetime.now()
    yr, mon, day = now.year-1900, now.month, now.day
    numrec = len(records)
    numfields = len(fieldspecs)
    lenheader = numfields * 32 + 33
    lenrecord = sum(field[1] for field in fieldspecs) + 1
    hdr = struct.pack('<BBBBLHH20x', ver, yr, mon, day, numrec, lenheader, lenrecord)
    f.write(hdr)
                      
    # field specs
    for name, (typ, size, deci) in itertools.izip(fieldnames, fieldspecs):
        name = name.ljust(11, '\x00')
        fld = struct.pack('<11sc4xBB14x', name, typ, size, deci)
        f.write(fld)

    # terminator
    f.write('\r')

    # records
    for record in records:
        f.write(' ')                        # deletion flag
        for (typ, size, deci), value in itertools.izip(fieldspecs, record):
            if typ == "N":
                value = str(value).rjust(size, ' ')
            elif typ == 'D':
                value = value.strftime('%Y%m%d')
            elif typ == 'L':
                value = str(value)[0].upper()
            else:
                value = str(value)[:size].ljust(size, ' ')
            assert len(value) == size
            f.write(value)

    # End of file
    f.write('\x1A')












# Shapefile utilities
# Original source is apparently:
# http://code.activestate.com/recipes/362715/
# 
# ...but got via here:
# http://indiemaps.com/blog/2008/03/easy-shapefile-loading-in-python/

from struct import unpack
#import dbfUtils, math
import math
XY_POINT_RECORD_LENGTH = 16
db = []

def loadShapefile(file_name):
	global db
	shp_bounding_box = []
	shp_type = 0
	file_name = file_name
	records = []
	# open dbf file and get records as a list
	dbf_file = file_name[0:-4] + '.dbf'
	dbf = open(dbf_file, 'rb')
	db = list(dbfUtils.dbfreader(dbf))
	dbf.close()
	fp = open(file_name, 'rb')
	
	# get basic shapefile configuration
	fp.seek(32)
	shp_type = readAndUnpack('i', fp.read(4))		
	shp_bounding_box = readBoundingBox(fp)
	
	# fetch Records
	fp.seek(100)
	while True:
		shp_record = createRecord(fp)
		if shp_record == False:
			break
		records.append(shp_record)
	
	return records

record_class = {0:'RecordNull', 1:'RecordPoint', 8:'RecordMultiPoint', 3:'RecordPolyLine', 5:'RecordPolygon'}

def createRecord(fp):
	# read header
	record_number = readAndUnpack('>L', fp.read(4))
	if record_number == '':
		print 'doner'
		return False
	content_length = readAndUnpack('>L', fp.read(4))
	record_shape_type = readAndUnpack('<L', fp.read(4))

	shp_data = readRecordAny(fp,record_shape_type)
	dbf_data = {}
	for i in range(0,len(db[record_number+1])):
		dbf_data[db[0][i]] = db[record_number+1][i]
	
	return {'shp_data':shp_data, 'dbf_data':dbf_data} 
	
# Reading defs

def readRecordAny(fp, type):
	if type==0:
		return readRecordNull(fp)
	elif type==1:
		return readRecordPoint(fp)
	elif type==8:
		return readRecordMultiPoint(fp)
	elif type==3 or type==5:
		return readRecordPolyLine(fp)
	else:
		return False

def readRecordNull(fp):
	data = {}
	return data

point_count = 0
def readRecordPoint(fp):
	global point_count
	data = {}
	data['x'] = readAndUnpack('d', fp.read(8))
	data['y'] = readAndUnpack('d', fp.read(8))
	point_count += 1
	return data


def readRecordMultiPoint(fp):
	data = readBoundingBox(fp)
	data['numpoints'] = readAndUnpack('i', fp.read(4))	
	for i in range(0,data['numpoints']):
		data['points'].append(readRecordPoint(fp))
	return data


def readRecordPolyLine(fp):
	data = readBoundingBox(fp)
	data['numparts'] = readAndUnpack('i', fp.read(4))
	data['numpoints'] = readAndUnpack('i', fp.read(4))
	data['parts'] = []
	for i in range(0, data['numparts']):
		data['parts'].append(readAndUnpack('i', fp.read(4)))
	points_initial_index = fp.tell()
	points_read = 0
	for part_index in range(0, data['numparts']):
		point_index = data['parts'][part_index]
		
		# if(!isset(data['parts'][part_index]['points']) or !is_array(data['parts'][part_index]['points'])):
		data['parts'][part_index] = {}
		data['parts'][part_index]['points'] = []
		
		# while( ! in_array( points_read, data['parts']) and points_read < data['numpoints'] and !feof(fp)):
		checkPoint = []
		while (points_read < data['numpoints']):
			currPoint = readRecordPoint(fp)
			data['parts'][part_index]['points'].append(currPoint)
			points_read += 1
			if points_read == 0 or checkPoint == []:
				checkPoint = currPoint
			elif currPoint == checkPoint:
				checkPoint = []
				break
			
	fp.seek(points_initial_index + (points_read * XY_POINT_RECORD_LENGTH))
	return data

# General defs

def readBoundingBox(fp):
	data = {}
	data['xmin'] = readAndUnpack('d',fp.read(8))
	data['ymin'] = readAndUnpack('d',fp.read(8))
	data['xmax'] = readAndUnpack('d',fp.read(8))
	data['ymax'] = readAndUnpack('d',fp.read(8))
	return data

def readAndUnpack(type, data):
	if data=='': return data
	return unpack(type, data)[0]


####
#### additional functions
####

def getCentroids(records, projected=False):
	# for each feature
	if projected:
		points = 'projectedPoints'
	else:
		points = 'points'
		
	for feature in records:
		numpoints = cx = cy = 0
		for part in feature['shp_data']['parts']:
			for point in part[points]:
				numpoints += 1
				cx += point['x']
				cy += point['y']
		cx /= numpoints
		cy /= numpoints
		feature['shp_data']['centroid'] = {'x':cx, 'y':cy}
				
		
def getBoundCenters(records):
	for feature in records:
		cx = .5 * (feature['shp_data']['xmax']-feature['shp_data']['xmin']) + feature['shp_data']['xmin']
		cy = .5 * (feature['shp_data']['ymax']-feature['shp_data']['ymin']) + feature['shp_data']['ymin']
		feature['shp_data']['boundCenter'] = {'x':cx, 'y':cy}
	
def getTrueCenters(records, projected=False):
	#gets the true polygonal centroid for each feature (uses largest ring)
	#should be spherical, but isn't

	if projected:
		points = 'projectedPoints'
	else:
		points = 'points'
		
	for feature in records:
		maxarea = 0
		for ring in feature['shp_data']['parts']:
			ringArea = getArea(ring, points)
			if ringArea > maxarea:
				maxarea = ringArea
				biggest = ring
		#now get the true centroid
		tempPoint = {'x':0, 'y':0}
		if biggest[points][0] != biggest[points][len(biggest[points])-1]:
			print "mug", biggest[points][0], biggest[points][len(biggest[points])-1]
		for i in range(0, len(biggest[points])-1):
			j = (i + 1) % (len(biggest[points])-1)
			tempPoint['x'] -= (biggest[points][i]['x'] + biggest[points][j]['x']) * ((biggest[points][i]['x'] * biggest[points][j]['y']) - (biggest[points][j]['x'] * biggest[points][i]['y']))
			tempPoint['y'] -= (biggest[points][i]['y'] + biggest[points][j]['y']) * ((biggest[points][i]['x'] * biggest[points][j]['y']) - (biggest[points][j]['x'] * biggest[points][i]['y']))
			
		tempPoint['x'] = tempPoint['x'] / ((6) * maxarea)
		tempPoint['y'] = tempPoint['y'] / ((6) * maxarea)
		feature['shp_data']['truecentroid'] = tempPoint
		

def getArea(ring, points):
	#returns the area of a polygon
	#needs to be spherical area, but isn't
	area = 0
	for i in range(0,len(ring[points])-1):
		j = (i + 1) % (len(ring[points])-1)
		area += ring[points][i]['x'] * ring[points][j]['y']
		area -= ring[points][i]['y'] * ring[points][j]['x']
			
	return math.fabs(area/2)
	

def getNeighbors(records):
	
	#for each feature
	for i in range(len(records)):
		#print i, records[i]['dbf_data']['ADMIN_NAME']
		if not 'neighbors' in records[i]['shp_data']:
			records[i]['shp_data']['neighbors'] = []
		
		#for each other feature
		for j in range(i+1, len(records)):
			numcommon = 0
			#first check to see if the bounding boxes overlap
			if overlap(records[i], records[j]):
				#if so, check every single point in this feature to see if it matches a point in the other feature
				
				#for each part:
				for part in records[i]['shp_data']['parts']:
					
					#for each point:
					for point in part['points']:
						
						for otherPart in records[j]['shp_data']['parts']:
							if point in otherPart['points']:
								numcommon += 1
								if numcommon == 2:
									if not 'neighbors' in records[j]['shp_data']:
										records[j]['shp_data']['neighbors'] = []
									records[i]['shp_data']['neighbors'].append(j)
									records[j]['shp_data']['neighbors'].append(i)
									#now break out to the next j
									break
						if numcommon == 2:
							break
					if numcommon == 2:
						break
				
									
								
								
def projectShapefile(records, whatProjection, lonCenter=0, latCenter=0):
	print 'projecting to ', whatProjection
	for feature in records:
		for part in feature['shp_data']['parts']:
			part['projectedPoints'] = []
			for point in part['points']:
				tempPoint = projectPoint(point, whatProjection, lonCenter, latCenter)
				part['projectedPoints'].append(tempPoint)

def projectPoint(fromPoint, whatProjection, lonCenter, latCenter):
	latRadians = fromPoint['y'] * math.pi/180
	if latRadians > 1.5: latRadians = 1.5
	if latRadians < -1.5: latRadians = -1.5
	lonRadians = fromPoint['x'] * math.pi/180
	lonCenter = lonCenter * math.pi/180
	latCenter = latCenter * math.pi/180
	newPoint = {}
	if whatProjection == "MERCATOR":
		newPoint['x'] = (180/math.pi) * (lonRadians - lonCenter)
		newPoint['y'] = (180/math.pi) * math.log(math.tan(latRadians) + (1/math.cos(latRadians)))
		if newPoint['y'] > 200:
			newPoint['y'] = 200
		if newPoint['y'] < -200:
			newPoint['y'] = 200
		return newPoint
	if whatProjection == "EQUALAREA":
		newPoint['x'] = 0
		newPoint['y'] = 0
		return newPoint
		

def overlap(feature1, feature2):
	if (feature1['shp_data']['xmax'] > feature2['shp_data']['xmin'] and feature1['shp_data']['ymax'] > feature2['shp_data']['ymin'] and feature1['shp_data']['xmin'] < feature2['shp_data']['xmax'] and feature1['shp_data']['ymin'] < feature2['shp_data']['ymax']):
		return True
	else:
		return False



	