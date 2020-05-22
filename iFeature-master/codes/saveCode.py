#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys

def savetsv(encodings, file = 'encoding.tsv'):
	print("hello world\n")
	with open(file, 'w') as f:
		if encodings == 0:
			f.write('Descriptor calculation failed.')
		else:
			for i in range(len(encodings[0]) - 1):
				f.write(encodings[0][i] + '\t')
			f.write(encodings[0][-1] + '\n')
			for i in encodings[1:]:
				f.write(i[0] + '\t')
				for j in range(1, len(i) - 1):
					f.write(str(float(i[j])) + '\t')
				try:
					#print(str(float(i[len(i)-1])))
					f.write(str(float(i[len(i)-1])) + '\n')
				except:
					#print(str(float(i[len(i)-1])))
					f.write("\n")
					continue
					f.write("why" + "\n")
				
	return None
