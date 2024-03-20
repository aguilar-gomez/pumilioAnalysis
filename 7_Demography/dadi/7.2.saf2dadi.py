#!/usr/bin/env python3
# Author :Debora Brandt, modified by Diana Aguilar
import sys

angsdfile=sys.argv[1]
n1=int(sys.argv[2])
n2=int(sys.argv[3])
p1=str(sys.argv[4])
p2=str(sys.argv[5])

print(2*n1+1, 2*n2+1, 'folded', f'"{p1}"',f'"{p2}"')
with open(angsdfile, 'r') as fh:
    for line in fh:
        print(line.strip())
#print(' '.join(['0' for i in range((2*n1+1) * (2*n2+1))]))
mask=[str(int(0==float(x))) for x in line.split()]
print(' '.join(mask))
