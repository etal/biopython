"""
Generic functions for simple operations.
"""
# Used by unescape
import unicodedata
import re, htmlentitydefs
from numpy import array, reshape, NaN


# Incorporates:
# Public domain module to convert diverse characters to ASCII (dammit!)
# http://newsbruiser.tigris.org/source/browse/~checkout~/newsbruiser/nb/lib/AsciiDammit.py
# from AsciiDammit import asciiDammit
"""ASCII, Dammit

Stupid library to turn MS chars (like smart quotes) and ISO-Latin
chars into ASCII, dammit. Will do plain text approximations, or more
accurate HTML representations. Can also be jiggered to just fix the
smart quotes and leave the rest of ISO-Latin alone.

Sources:
 http://www.cs.tut.fi/~jkorpela/latin1/all.html
 http://www.webreference.com/html/reference/character/isolat1.html

1.0 Initial Release (2004-11-28)

The author hereby irrevocably places this work in the public domain.
To the extent that this statement does not divest the copyright,
the copyright holder hereby grants irrevocably to every recipient
all rights in this work otherwise reserved under copyright.
"""

__author__ = "Leonard Richardson (leonardr@segfault.org)"
__version__ = "$Revision: 1.3 $"
__date__ = "$Date: 2009/04/28 10:45:03 $"
__license__ = "Public domain"

import re
import string
import types

CHARS = { '\x80' : ('EUR', 'euro'),
          '\x81' : ' ',
          '\x82' : (',', 'sbquo'),
          '\x83' : ('f', 'fnof'),
          '\x84' : (',,', 'bdquo'),
          '\x85' : ('...', 'hellip'),
          '\x86' : ('+', 'dagger'),
          '\x87' : ('++', 'Dagger'),
          '\x88' : ('^', 'caret'),
          '\x89' : '%',
          '\x8A' : ('S', 'Scaron'),
          '\x8B' : ('<', 'lt;'),
          '\x8C' : ('OE', 'OElig'),
          '\x8D' : '?',
          '\x8E' : 'Z',
          '\x8F' : '?',
          '\x90' : '?',
          '\x91' : ("'", 'lsquo'),
          '\x92' : ("'", 'rsquo'),
          '\x93' : ('"', 'ldquo'),
          '\x94' : ('"', 'rdquo'),
          '\x95' : ('*', 'bull'),
          '\x96' : ('-', 'ndash'),
          '\x97' : ('--', 'mdash'),
          '\x98' : ('~', 'tilde'),
          '\x99' : ('(TM)', 'trade'),
          '\x9a' : ('s', 'scaron'),
          '\x9b' : ('>', 'gt'),
          '\x9c' : ('oe', 'oelig'),
          '\x9d' : '?',
          '\x9e' : 'z',
          '\x9f' : ('Y', 'Yuml'),
          '\xa0' : (' ', 'nbsp'),
          '\xa1' : ('!', 'iexcl'),
          '\xa2' : ('c', 'cent'),
          '\xa3' : ('GBP', 'pound'),
          '\xa4' : ('$', 'curren'), #This approximation is especially lame.
          '\xa5' : ('YEN', 'yen'),
          '\xa6' : ('|', 'brvbar'),
          '\xa7' : ('S', 'sect'),
          '\xa8' : ('..', 'uml'),
          '\xa9' : ('', 'copy'),
          '\xaa' : ('(th)', 'ordf'),
          '\xab' : ('<<', 'laquo'),
          '\xac' : ('!', 'not'),
          '\xad' : (' ', 'shy'),
          '\xae' : ('(R)', 'reg'),
          '\xaf' : ('-', 'macr'),
          '\xb0' : ('o', 'deg'),
          '\xb1' : ('+-', 'plusmm'),
          '\xb2' : ('2', 'sup2'),
          '\xb3' : ('3', 'sup3'),
          '\xb4' : ("'", 'acute'),
          '\xb5' : ('u', 'micro'),
          '\xb6' : ('P', 'para'),
          '\xb7' : ('*', 'middot'),
          '\xb8' : (',', 'cedil'),
          '\xb9' : ('1', 'sup1'),
          '\xba' : ('(th)', 'ordm'),
          '\xbb' : ('>>', 'raquo'),
          '\xbc' : ('1/4', 'frac14'),
          '\xbd' : ('1/2', 'frac12'),
          '\xbe' : ('3/4', 'frac34'),
          '\xbf' : ('?', 'iquest'),          
          '\xc0' : ('A', "Agrave"),
          '\xc1' : ('A', "Aacute"),
          '\xc2' : ('A', "Acirc"),
          '\xc3' : ('A', "Atilde"),
          '\xc4' : ('A', "Auml"),
          '\xc5' : ('A', "Aring"),
          '\xc6' : ('AE', "Aelig"),
          '\xc7' : ('C', "Ccedil"),
          '\xc8' : ('E', "Egrave"),
          '\xc9' : ('E', "Eacute"),
          '\xca' : ('E', "Ecirc"),
          '\xcb' : ('E', "Euml"),
          '\xcc' : ('I', "Igrave"),
          '\xcd' : ('I', "Iacute"),
          '\xce' : ('I', "Icirc"),
          '\xcf' : ('I', "Iuml"),
          '\xd0' : ('D', "Eth"),
          '\xd1' : ('N', "Ntilde"),
          '\xd2' : ('O', "Ograve"),
          '\xd3' : ('O', "Oacute"),
          '\xd4' : ('O', "Ocirc"),
          '\xd5' : ('O', "Otilde"),
          '\xd6' : ('O', "Ouml"),
          '\xd7' : ('*', "times"),
          '\xd8' : ('O', "Oslash"),
          '\xd9' : ('U', "Ugrave"),
          '\xda' : ('U', "Uacute"),
          '\xdb' : ('U', "Ucirc"),
          '\xdc' : ('U', "Uuml"),
          '\xdd' : ('Y', "Yacute"),
          '\xde' : ('b', "Thorn"),
          '\xdf' : ('B', "szlig"),
          '\xe0' : ('a', "agrave"),
          '\xe1' : ('a', "aacute"),
          '\xe2' : ('a', "acirc"),
          '\xe3' : ('a', "atilde"),
          '\xe4' : ('a', "auml"),
          '\xe5' : ('a', "aring"),
          '\xe6' : ('ae', "aelig"),
          '\xe7' : ('c', "ccedil"),
          '\xe8' : ('e', "egrave"),
          '\xe9' : ('e', "eacute"),
          '\xea' : ('e', "ecirc"),
          '\xeb' : ('e', "euml"),
          '\xec' : ('i', "igrave"),
          '\xed' : ('i', "iacute"),
          '\xee' : ('i', "icirc"),
          '\xef' : ('i', "iuml"),
          '\xf0' : ('o', "eth"),
          '\xf1' : ('n', "ntilde"),
          '\xf2' : ('o', "ograve"),
          '\xf3' : ('o', "oacute"),
          '\xf4' : ('o', "ocirc"),
          '\xf5' : ('o', "otilde"),
          '\xf6' : ('o', "ouml"),
          '\xf7' : ('/', "divide"),
          '\xf8' : ('o', "oslash"),
          '\xf9' : ('u', "ugrave"),
          '\xfa' : ('u', "uacute"),
          '\xfb' : ('u', "ucirc"),
          '\xfc' : ('u', "uuml"),
          '\xfd' : ('y', "yacute"),
          '\xfe' : ('b', "thorn"),
          '\xff' : ('y', "yuml"),
          }
def _makeRE(limit):
    """Returns a regular expression object that will match special characters
    up to the given limit."""
    return re.compile("([\x80-\\x%s])" % limit, re.M)
ALL = _makeRE('ff')
ONLY_WINDOWS = _makeRE('9f')

def _replHTML(match):
    "Replace the matched character with its HTML equivalent."
    return _repl(match, 1)
          
def _repl(match, html=0):
    "Replace the matched character with its HTML or ASCII equivalent."
    g = match.group(0)
    a = CHARS.get(g,g)
    if type(a) == types.TupleType:
        a = a[html]
        if html:
            a = '&' + a + ';'
    return a

def _dammit(t, html=0, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into an ASCII representation, dammit."

    r = ALL
    if fixWindowsOnly:
        r = ONLY_WINDOWS
    m = _repl
    if html:
        m = _replHTML

    return re.sub(r, m, t)

def asciiDammit(t, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into a plain ASCII approximation, dammit."
    return _dammit(t, 0, fixWindowsOnly)

def htmlDammit(t, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into plain ASCII with HTML codes, dammit."
    return _dammit(t, 1, fixWindowsOnly=fixWindowsOnly)

def demoronise(t):
    """Helper method named in honor of the original smart quotes
    remover, The Demoroniser:

    http://www.fourmilab.ch/webtools/demoroniser/"""
    return asciiDammit(t, 1)

if __name__ == '__main__':

    french = '\x93Sacr\xe9 bleu!\x93'
    print "First we mangle some French."
    print asciiDammit(french)
    print htmlDammit(french)

    print
    print "And now we fix the MS-quotes but leave the French alone."
    print demoronise(french)
    print htmlDammit(french, 1)




# ==========================
# non-AsciiDammit functions (by Nick Matzke) follow...

def make_NaN_array(xdim, ydim):
	"""
	Make an empty floating-point array with the specified dimensions
	"""
	temparray = array( [NaN] * (xdim * ydim), dtype=float)
	temparray = reshape(temparray, (xdim, ydim))
	return temparray
	
	
def make_None_list_array(xdim, ydim):
	"""
	Make a list of lists ("array") with the specified dimensions
	"""
	
	# Create empty array
	temp_row_array = [None] * xdim
	temp_cols_array = []
	for i in range(0, ydim):
		temp_cols_array.append(temp_row_array)
	
	return temp_cols_array

def set_diags_to_none(list_array):
	"""
	Given a square list of lists (rows), set the diagonal to None.
	"""
	ydim = len(list_array)
	xdim = len(list_array[0])
	
	if xdim != ydim:
		print ''
		print 'ERROR: This list_array is not a square!'
	
	for index1 in range(0, ydim-1):
		for index2 in range(0, xdim-1):
			list_array[index1, index2] = None
	
	return list_array


def list1_items_in_list2(list1, list2):
	"""
	Returns the list of list1 items which are found in list2
	http://vermeulen.ca/python-techniques.html
	"""
	intersection = filter(lambda x:x in list1, list2)
	return intersection

def list1_items_not_in_list2(list1, list2):
	"""
	Returns the list of list1 items which are NOT found in list2
	http://vermeulen.ca/python-techniques.html
	"""
	difference=filter(lambda x:x not in list2,list1)
	return difference


def remove_nonxml_brackets(line):
	"""
	When a GBIF DarwinCore-formatted XML record is converted to ASCII, links in some of the metadata (e.g. "<rap.conservation.org>") can get interpreted as XML tags, which are then mismatched, causing a crash.  This converts those to e.g. "[rap.conservation.org]".
	"""
	

	
	match_strings = re.findall(r'<.*?>', line)
	
	list_of_valid_xml = ['<gbif', '</gbif', '<?xml', '<to:', '</to:', '<tc:', '</tc:', '<tn:', '</tn:', '<tcom:', '</tcom:', '<tax', '</tax', '<data', '</data', '<name', '</name', '<occurrence', '</occurrence' ]
	
	
	# If no match, continue to next loop
	if match_strings != []:		
		# If there is a match:
		for tempmatch in match_strings:
			# Skip legit tags
			valid_tag_found = False
			for valid_tag in list_of_valid_xml:
				if tempmatch.startswith(valid_tag):
					# if it matches, go to next tempmatch
					valid_tag_found = True
			if valid_tag_found == False:
				# if not, replace
				replacement = tempmatch.replace('<', '[')
				replacement = replacement.replace('>', ']')
				print 'Replacing "' + tempmatch + '" with "' + replacement + '"'
				line = line.replace(tempmatch, replacement)
				continue
			
		
	return line
				
		

def unescape(text):
	"""
	##
	# Removes HTML or XML character references and entities from a text string.
	#
	# @param text The HTML (or XML) source text.
	# @return The plain text, as a Unicode string, if necessary.
	# source: http://effbot.org/zone/re-sub.htm#unescape-html
	"""
	def fixup(m):
		text = m.group(0)
		if text[:2] == "&#":
			# character reference
			try:
				if text[:3] == "&#x":
					return unichr(int(text[3:-1], 16))
				else:
					return unichr(int(text[2:-1]))
			except ValueError:
				pass
		else:
			# named entity
			try:
				text = unichr(htmlentitydefs.name2codepoint[text[1:-1]])
			except KeyError:
				pass
		return text # leave as is
	return re.sub("&#?\w+;", fixup, text)
	

def fix_ampersand(line):
	"""
	Replaces "&" with "&amp;" in a string; this is otherwise 
	not caught by the unescape and unicodedata.normalize functions.
	"""
	return line.replace('&', '&amp;')



def fix_ASCII_file(fn_unfixed):
	"""
	# Search-replace to fix annoying 
	# non-ASCII characters in search results
	# 
	# inspiration:
	# http://www.amk.ca/python/howto/unicode
	# http://www.peterbe.com/plog/unicode-to-ascii
	"""
	if fn_unfixed.count('.') < 1:
		return None
	
	# Make a new filename, insert '_fixed' before the rightmost period
	fn_fixed_tuple = fn_unfixed.rpartition('.')
	fn_fixed = ''.join( [fn_fixed_tuple[0], '_fixed.', fn_fixed_tuple[2]] )

	# Open unfixed file and fix all non-ASCII characters
	fh = open(fn_unfixed, 'r')
	fh_fixed = open(fn_fixed, 'w')

	lines = fh.readlines()
	newlines = fix_ASCII_lines(lines)
	
	fh_fixed.write(newlines)
	fh_fixed.close()
	fh.close()
	return fn_fixed




def fix_ASCII_lines(lines):
	"""
	Convert each line in an input string into pure ASCII
	(This avoids crashes when printing to screen, etc.)
	"""
	# If type is string:
	if lines.__class__ == 'string'.__class__:
		lines2 = lines.split('\n')
	# If type is list:
	elif lines.__class__ == [1,2].__class__:
		lines2 = lines
	else:
		print ''
		print 'Error: fix_ASCII_lines takes only strings with endlines, or lists.'
	
	newstr_list = []
	for line in lines2:
		ascii_content5 = fix_ASCII_line(line)
		
		newstr_list.append(ascii_content5 + '\n')
	return ''.join(newstr_list)


def fix_ASCII_line(line):
	"""
	Run several functions that fix common problems found when printing e.g. GBIF's XML results to screen: non-ASCII characters, ampersands, and <links in brackets> which can get misinterpreted as unmatched XML tags.
	"""
	
	# Catch error of an empty string.
	if line == None:
		line = "None"
	if line == '':
		line = ''
	
	# library from here: http://effbot.org/zone/re-sub.htm#unescape-html
	ascii_content1 = unescape(line)
	
	# Public domain module to convert diverse characters to ASCII
	# http://newsbruiser.tigris.org/source/browse/~checkout~/newsbruiser/nb/lib/AsciiDammit.py
	ascii_content2 = asciiDammit(ascii_content1)
	
	# inspiration: http://www.amk.ca/python/howto/unicode
	ascii_content3 = unicodedata.normalize('NFKC', unicode(ascii_content2)).encode('ascii','ignore')
	
	# Fix the ampersand
	ascii_content4 = fix_ampersand(ascii_content3)
	ascii_content5 = remove_nonxml_brackets(ascii_content4)
	
	return ascii_content5

def element_text_to_string(txt):
	if txt == None:
		txt = ""
	return str(txt).strip()


def element_items_to_string(items):
	"""
	Input a list of items, get string back
	"""
	s = ""
	for item in items:
		s = s + " " + str(item)
		s = s.strip()
	
	return s




def get_str_subset(start, end, seq):
	
	index1 = start-1
	index2 = end
	
	newstring = str()
	
	for i in range(index1, index2):
		newstring = newstring + seq[i]
	
	#print 'len(newstring)=', str(len(newstring))
	#print 'start - stop', str(end - start + 1) 
	
	return newstring



def find_1st_match(list1, list2):
	"""
	Find the first match in two ordered lists.
	"""
	anc_list1 = list1
	anc_list2 = list2
	
	for anc1 in anc_list1:
		for anc2 in anc_list2:
			if anc1==anc2:
				return(anc1)
	
	return None



	


"""ASCII, Dammit

Stupid library to turn MS chars (like smart quotes) and ISO-Latin
chars into ASCII, dammit. Will do plain text approximations, or more
accurate HTML representations. Can also be jiggered to just fix the
smart quotes and leave the rest of ISO-Latin alone.

Sources:
 http://www.cs.tut.fi/~jkorpela/latin1/all.html
 http://www.webreference.com/html/reference/character/isolat1.html

1.0 Initial Release (2004-11-28)

The author hereby irrevocably places this work in the public domain.
To the extent that this statement does not divest the copyright,
the copyright holder hereby grants irrevocably to every recipient
all rights in this work otherwise reserved under copyright.
"""

__author__ = "Leonard Richardson (leonardr@segfault.org)"
__version__ = "$Revision: 1.3 $"
__date__ = "$Date: 2009/04/28 10:45:03 $"
__license__ = "Public domain"

import re
import string
import types

CHARS = { '\x80' : ('EUR', 'euro'),
          '\x81' : ' ',
          '\x82' : (',', 'sbquo'),
          '\x83' : ('f', 'fnof'),
          '\x84' : (',,', 'bdquo'),
          '\x85' : ('...', 'hellip'),
          '\x86' : ('+', 'dagger'),
          '\x87' : ('++', 'Dagger'),
          '\x88' : ('^', 'caret'),
          '\x89' : '%',
          '\x8A' : ('S', 'Scaron'),
          '\x8B' : ('<', 'lt;'),
          '\x8C' : ('OE', 'OElig'),
          '\x8D' : '?',
          '\x8E' : 'Z',
          '\x8F' : '?',
          '\x90' : '?',
          '\x91' : ("'", 'lsquo'),
          '\x92' : ("'", 'rsquo'),
          '\x93' : ('"', 'ldquo'),
          '\x94' : ('"', 'rdquo'),
          '\x95' : ('*', 'bull'),
          '\x96' : ('-', 'ndash'),
          '\x97' : ('--', 'mdash'),
          '\x98' : ('~', 'tilde'),
          '\x99' : ('(TM)', 'trade'),
          '\x9a' : ('s', 'scaron'),
          '\x9b' : ('>', 'gt'),
          '\x9c' : ('oe', 'oelig'),
          '\x9d' : '?',
          '\x9e' : 'z',
          '\x9f' : ('Y', 'Yuml'),
          '\xa0' : (' ', 'nbsp'),
          '\xa1' : ('!', 'iexcl'),
          '\xa2' : ('c', 'cent'),
          '\xa3' : ('GBP', 'pound'),
          '\xa4' : ('$', 'curren'), #This approximation is especially lame.
          '\xa5' : ('YEN', 'yen'),
          '\xa6' : ('|', 'brvbar'),
          '\xa7' : ('S', 'sect'),
          '\xa8' : ('..', 'uml'),
          '\xa9' : ('', 'copy'),
          '\xaa' : ('(th)', 'ordf'),
          '\xab' : ('<<', 'laquo'),
          '\xac' : ('!', 'not'),
          '\xad' : (' ', 'shy'),
          '\xae' : ('(R)', 'reg'),
          '\xaf' : ('-', 'macr'),
          '\xb0' : ('o', 'deg'),
          '\xb1' : ('+-', 'plusmm'),
          '\xb2' : ('2', 'sup2'),
          '\xb3' : ('3', 'sup3'),
          '\xb4' : ("'", 'acute'),
          '\xb5' : ('u', 'micro'),
          '\xb6' : ('P', 'para'),
          '\xb7' : ('*', 'middot'),
          '\xb8' : (',', 'cedil'),
          '\xb9' : ('1', 'sup1'),
          '\xba' : ('(th)', 'ordm'),
          '\xbb' : ('>>', 'raquo'),
          '\xbc' : ('1/4', 'frac14'),
          '\xbd' : ('1/2', 'frac12'),
          '\xbe' : ('3/4', 'frac34'),
          '\xbf' : ('?', 'iquest'),          
          '\xc0' : ('A', "Agrave"),
          '\xc1' : ('A', "Aacute"),
          '\xc2' : ('A', "Acirc"),
          '\xc3' : ('A', "Atilde"),
          '\xc4' : ('A', "Auml"),
          '\xc5' : ('A', "Aring"),
          '\xc6' : ('AE', "Aelig"),
          '\xc7' : ('C', "Ccedil"),
          '\xc8' : ('E', "Egrave"),
          '\xc9' : ('E', "Eacute"),
          '\xca' : ('E', "Ecirc"),
          '\xcb' : ('E', "Euml"),
          '\xcc' : ('I', "Igrave"),
          '\xcd' : ('I', "Iacute"),
          '\xce' : ('I', "Icirc"),
          '\xcf' : ('I', "Iuml"),
          '\xd0' : ('D', "Eth"),
          '\xd1' : ('N', "Ntilde"),
          '\xd2' : ('O', "Ograve"),
          '\xd3' : ('O', "Oacute"),
          '\xd4' : ('O', "Ocirc"),
          '\xd5' : ('O', "Otilde"),
          '\xd6' : ('O', "Ouml"),
          '\xd7' : ('*', "times"),
          '\xd8' : ('O', "Oslash"),
          '\xd9' : ('U', "Ugrave"),
          '\xda' : ('U', "Uacute"),
          '\xdb' : ('U', "Ucirc"),
          '\xdc' : ('U', "Uuml"),
          '\xdd' : ('Y', "Yacute"),
          '\xde' : ('b', "Thorn"),
          '\xdf' : ('B', "szlig"),
          '\xe0' : ('a', "agrave"),
          '\xe1' : ('a', "aacute"),
          '\xe2' : ('a', "acirc"),
          '\xe3' : ('a', "atilde"),
          '\xe4' : ('a', "auml"),
          '\xe5' : ('a', "aring"),
          '\xe6' : ('ae', "aelig"),
          '\xe7' : ('c', "ccedil"),
          '\xe8' : ('e', "egrave"),
          '\xe9' : ('e', "eacute"),
          '\xea' : ('e', "ecirc"),
          '\xeb' : ('e', "euml"),
          '\xec' : ('i', "igrave"),
          '\xed' : ('i', "iacute"),
          '\xee' : ('i', "icirc"),
          '\xef' : ('i', "iuml"),
          '\xf0' : ('o', "eth"),
          '\xf1' : ('n', "ntilde"),
          '\xf2' : ('o', "ograve"),
          '\xf3' : ('o', "oacute"),
          '\xf4' : ('o', "ocirc"),
          '\xf5' : ('o', "otilde"),
          '\xf6' : ('o', "ouml"),
          '\xf7' : ('/', "divide"),
          '\xf8' : ('o', "oslash"),
          '\xf9' : ('u', "ugrave"),
          '\xfa' : ('u', "uacute"),
          '\xfb' : ('u', "ucirc"),
          '\xfc' : ('u', "uuml"),
          '\xfd' : ('y', "yacute"),
          '\xfe' : ('b', "thorn"),
          '\xff' : ('y', "yuml"),
          }

def _makeRE(limit):
    """Returns a regular expression object that will match special characters
    up to the given limit."""
    return re.compile("([\x80-\\x%s])" % limit, re.M)
ALL = _makeRE('ff')
ONLY_WINDOWS = _makeRE('9f')

def _replHTML(match):
    "Replace the matched character with its HTML equivalent."
    return _repl(match, 1)
          
def _repl(match, html=0):
    "Replace the matched character with its HTML or ASCII equivalent."
    g = match.group(0)
    a = CHARS.get(g,g)
    if type(a) == types.TupleType:
        a = a[html]
        if html:
            a = '&' + a + ';'
    return a

def _dammit(t, html=0, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into an ASCII representation, dammit."

    r = ALL
    if fixWindowsOnly:
        r = ONLY_WINDOWS
    m = _repl
    if html:
        m = _replHTML

    return re.sub(r, m, t)

def asciiDammit(t, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into a plain ASCII approximation, dammit."
    return _dammit(t, 0, fixWindowsOnly)

def htmlDammit(t, fixWindowsOnly=0):
    "Turns ISO-Latin-1 into plain ASCII with HTML codes, dammit."
    return _dammit(t, 1, fixWindowsOnly=fixWindowsOnly)

def demoronise(t):
    """Helper method named in honor of the original smart quotes
    remover, The Demoroniser:

    http://www.fourmilab.ch/webtools/demoroniser/"""
    return asciiDammit(t, 1)

if __name__ == '__main__':

    french = '\x93Sacr\xe9 bleu!\x93'
    print "First we mangle some French."
    print asciiDammit(french)
    print htmlDammit(french)

    print
    print "And now we fix the MS-quotes but leave the French alone."
    print demoronise(french)
    print htmlDammit(french, 1)
