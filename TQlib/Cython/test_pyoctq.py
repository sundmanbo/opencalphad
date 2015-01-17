from pyoctq import tqini, tqrfil, tqgpn, tqsetc, tqce, tqgetv, ctuple, cname, cnel
import numpy as np
import pyoctq as tq

n = 0
xf = np.array([.5,.5])
pxf = np.zeros((10,),dtype=np.double)
mu = np.zeros((10,),dtype=np.double)
phnames = [""]
tqini(n)
tqrfil('FENI.TDB')
formatter = "System with %i elements: " 
print formatter % cnel() +str([ cname(i) for i in range(0,cnel()) ])
print "and %i phases: " % ctuple()
phnames = [""]
for i in range(1,ctuple()+1):
    phnames.insert(i, tqgpn(i))
print phnames[1:]
n1 = 0
n2 = 0
# numerical values of conditions T, P, N
tp = [1e3,1e5,1.]
# set conditions for temperature, pressure and amount
cond = ["T", "P", "N"]
cnum = []
cnum.append(tqsetc(cond[0], n1,n2, tp[0]))
cnum.append(tqsetc(cond[1], n1, n2, tp[1]))
cnum.append(tqsetc(cond[2], n1, n2, tp[2]))
#for i in range(1,3):
#    xf[i] = 1.0/2.0

for i in range(1,2):
    cond = "X"
    cnum.append(tqsetc(cond,i,n2,xf[i]))
# calculate equilibria
target = " "
n1 = 0
n2 = 0
value = .0
tqce(target, n1, n2, value) 
statvar = "NP"
npf = np.array([0.1,0.2,0.3,0.5])
n1 = -1
n2 = 0
stable_ph = tqgetv(statvar,n1,n2, npf.size, npf)
print "Amount of %i phases:" % (stable_ph)
for i in npf[:stable_ph]:
    print "%.9f" % i

for i in range(1, ctuple()+1):
    print "Phase : %s " % phnames[i]
    if npf[i] > 0.0:
        print "Stable phase : %s, amount: %lf" % (phnames[i], npf[i])
        # use phase tuple index i
        statvar = "X"
        n2 = -1
        n4 = tqgetv(statvar, i, n2, pxf.size, pxf)
        for k in range(0,n4):
            print " %s : %lf, " % (cname(k), pxf[k])

print "Component, mole fraction and chemical potential"
for i in range(1, cnel()+1):
    statvar = "MU"
    n2 = 0
    n4 = pxf.size
    tqgetv(statvar, i, n2, n4, pxf)
    mu[i] = pxf[0]
    statvar = "X"
    tqgetv(statvar, i, n2, n4, pxf)
    print "%s        %lf        %lf" % (cname(i-1), pxf[0], mu[i]) 

