def complement(char):
	if char == 'A':
		return 'T'
	elif char == 'T':
		return 'A'
	elif char == 'C':
		return 'G'
	elif char == 'G':
		return 'C'
	elif char == 'N':
		return 'N'
	else:
		exit("ERROR: Invalid DNA string")

def complement_string(string):
	return ''.join([complement(x) for x in string])

def reverse_complement(string):
	return complement_string(string)[::-1]
