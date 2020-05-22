#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def TPC(fastas, **kw):
	AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWYX'
	encodings = []
	triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
	header = ['#'] + triPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i
	print(AADict)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 9261
		for j in range(len(sequence) - 3 + 1):
			print('j = %d' %j)
			print(sequence)
			tmpCode[AADict[sequence[j]] * 441 + AADict[sequence[j+1]]*21 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 441 + AADict[sequence[j+1]]*21 + AADict[sequence[j+2]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings
