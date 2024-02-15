import sys

f1 = open('out.txt', 'w')
f2 = open("warning.txt", "w")
sys.stdout = f1
sys.stderr = f2

f1.close()
f2.close()
