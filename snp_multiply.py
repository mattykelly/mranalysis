import MapReduce
import sys

"""
Matrix Multiply Example in the Simple Python MapReduce Framework
"""

mr = MapReduce.MapReduce()

# =============================
# Do not modify above this line

n = 4

def mapper(record):
	# key: position in nxn matrix
	# value: row
	for i in range(1, n+1):
		for j in range(len(record[1])):
			mr.emit_intermediate((record[0], i), ("a", record[0], j, record[1][j]))
			mr.emit_intermediate((i, record[0]), ("b", j, record[0], record[1][j]))

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
				sum += int(ra[3]) * int(rb[3]); 
	mr.emit((key[0], key[1], sum))

# Do not modify below this line
# =============================
if __name__ == '__main__':
	inputdata = open('data/wide_data.json')
	mr.execute(inputdata, mapper, reducer)
