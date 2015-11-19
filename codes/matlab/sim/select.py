import sys

filename1=sys.argv[1]
filename2=sys.argv[2]
f1=open(filename1)
f2=open(filename2)

headers = []
for line2 in f2:
	if line2.strip() != '':
		headers.append(line2.strip())

dummy = False
for line1 in f1:
	
	line1 = line1.strip()

	if line1.startswith('>'):
		t1=line1[1:].strip().split('|')

		if t1[0] in headers:
			print('>%s_unspliced cdna:known chromosome:GRCh37:%s:%s:%s:%s gene:%s' % (t1[0],t1[1],t1[2],t1[3],t1[4],t1[0]))
			dummy = True
		else:
			dummy = False

	else:

		if dummy:
			print line1



f1.close()
f2.close()
