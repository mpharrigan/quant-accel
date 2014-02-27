
with open('redo_all.sh', 'w') as f:
	for i in range(1,101):
		f.write('./redo-one.sh %d\n' % i)
