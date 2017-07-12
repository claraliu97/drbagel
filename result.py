#parse result
#Clara Liu
#2017/7/6

from parameters import *

def parse_result(str,num):
	counter = 0
	result = {}
	while counter < num:
		gene_name = str[str.find('\n')+1:str.find(':')]
		a1 = str.find('= ')
		a2 = str.find('(')
		b1 = str.find('= ',a2)
		b2 = str.find(')',a2)
		c1 = str.find('= ',b2)
		c2 = str.find('(',b2)
		d1 = str.find('= ',c2)
		d2 = str.find(')',c2)
		a = int(str[a1+1:a2])
		b = int(str[b1+1:b2])
		c = int(str[c1+1:c2])
		d = int(str[d1+1:d2])
		counter +=1
		result[gene_name] = [a,b,c,d]
		str = str[str.find('\n',1):]
	return result

def print_dict(dict):
	for item in dict:
		print "%s %s %s %s %s" %(item,dict[item][0],dict[item][1],dict[item][2],dict[item][3])

with open(species+'/Second_result.txt','r') as myfile:
	data = myfile.read()
	dict = parse_result(data,21)
	print_dict(dict)