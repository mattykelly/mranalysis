import MapReduce
import sys

"""
Matrix Multiply Example in the Simple Python MapReduce Framework
"""

mr = MapReduce.MapReduce()

# =============================
# Do not modify above this line

n = 5

def mapper(record):
	# key: position in nxn matrix
	# value: row
	if record[0] == 'a': 
		for i in range(n):
			mr.emit_intermediate((record[1], i), record)
	if record[0] == 'b':
		for i in range(n):
			mr.emit_intermediate((i, record[2]), record)

def reducer(key, list_of_values):
	# key: position on nxn matrix
	# value: list of records that may have an effect on the position
	sum = 0
	a = []
	b = []
	for record in list_of_values:
		a.append(record) if record[0] == 'a' else b.append(record)
	for ra in a:
		for rb in b:
			if ra[2] == rb[1]:
				sum += ra[3] * rb[3]; 
	mr.emit((key[0], key[1], sum))

# Do not modify below this line
# =============================
if __name__ == '__main__':
	inputdata = open('data/matrix.json')
	mr.execute(inputdata, mapper, reducer)
