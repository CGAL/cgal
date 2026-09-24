from sys import argv
from sys import exit
from sys import stderr

try:
  f=file(argv[1])
  x=float(argv[2])
  y=float(argv[3])
  z=float(argv[4])
except:
  stderr.write("Usage: "+argv[0]+" file.off x y z\n")
  exit()


print f.readline(),
line=f.readline()
print line,
print f.readline(),
for i in range(0,int(line.split()[0])):
  line=f.readline()
  line=line.split()
  print float(line[0])+x, float(line[1])+y, float(line[2])+z
  
line=f.readline()
while line:
  print line,
  line=f.readline()