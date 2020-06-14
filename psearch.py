import sys

for n in range(int(0.001*1000),int(9.999*1000),int(0.001*1000)):
    print(str(0.714+(float(n)/1000000)))
    sys.argv = ['simpleReader_cli.py', '-fJ0332+5434_194_chB.sdf', '-p'+str(0.714+(float(n)/1000000)), '-n91', '-d26.7787']
    #print(sys.argv)
    execfile('simpleReader_cli.py')
