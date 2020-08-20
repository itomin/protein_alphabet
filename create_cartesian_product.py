import csv


with open('cath.csv') as my_csv_file:
	reader = csv.reader(my_csv_file, delimiter = ';')
	cath = list(reader)

	
for i in range(0, len(cath)):
	iter = ["{}:{};{}:{}".format(cath[i][4],cath[i][5],row[4],row[5]) for row in cath[i:] if cath[i][2] != row[2]]
	with open('cath_pairs.csv', 'a') as myfile:
		myfile.write("\n".join(iter))
		myfile.write("\n")
	if(i%100==0):
		print(i)
		