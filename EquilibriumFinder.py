#Henry Heffan. This program finds the equilibrium of a system of equilibrium equations for chemical reactions
#Copyright (C) 2016  Henry Heffan

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import random
import itertools
import numpy as np

COMPLEATION_EXP=100
COMP_MAX=10**COMPLEATION_EXP
COMP_MIN=10**(-COMPLEATION_EXP)


def Q(equilibriumEquation,molarityMap,x=0):#if dividing by zero, return Q=1E200, no K's above this please
    top=1
    bottom=1
    reactants, products, Keq = equilibriumEquation
    for reac in reactants:
        bottom=bottom*(molarityMap[reac[0]]-x*reac[1])**reac[1]
    for prod in products:
        top=top*(molarityMap[prod[0]]+x*prod[1])**prod[1]
    Q=1E200*(1 if top==0 else top) if bottom==0 else top/bottom
    #print "(bottom,top,x,molarityMap,Q)="+str((bottom,top,x,molarityMap,Q))
    return Q

def getShiftToEq(e,M,xMin,xMax,maxSerchRounds=1000):
    products,reactants,Keq=e
    #since Q monitonicly increasing or decreasing on xMin to xMax
    xSmall=xMin*.75+xMax*.25
    xBig=xMin*.25+xMax*.75
    #y=Q=\prod_{all products p}(M[p]+C_px)^C_p
    #print "calcing QxS then QxB="+str((xSmall,xBig))
    QxS=Q(e,M,xSmall)
    QxB=Q(e,M,xBig)
    #print "(xS,xB,qS,qB)="+str((xSmall,xBig,QxS,QxB,e,M))
    increasing=QxS<=QxB#probably triavial. will assume so for now

    Qmin=Q(e,M,xMin)
    Qmax=Q(e,M,xMax)
    assert increasing, "not increasing  (xS,xB,qS,qB)="+str((xSmall,xBig,QxS,QxB,Keq,e,M))
    assert Qmin<=Keq<=Qmax or (Qmin==Qmax and Qmin==1E200), "Keq not in range (xS,xB,qS,qB)="+str((xSmall,xBig,QxS,QxB,Keq,e,M))
    if(Keq==Qmin or Keq<=COMP_MIN):
        return xMin
    if(Keq==Qmax or Keq>=COMP_MAX):
        return xMax
    xMinL=1E250
    xMaxL=-1E250
    for i in range(0,maxSerchRounds): #change to distance from Keq?
        xTest=(xMin+xMax)/2.0
        Qt=Q(e,M,xTest)
        if(xMin==xMinL and xMax==xMaxL):
            return xTest
        xMinL=xMin
        xMaxL=xMax
        if(Keq==Qt):
            return xTest
        elif(Keq<Qt):
            xMax=xTest
            Qmax=Qt
        elif(Keq>Qt):
            xMin=xTest
            Qmin=Qt
        else:
            print "error"
        #print "step LS= "+str(i)+"(xMin,xMax)= "+str((xMin,xMax,(xMin+xMax)/2.0))
    return (xMin+xMax)/2.0



def shiftMap(equilibriumEquation,molarityMap,N,maxSerchRounds=1000):
    #print equilibriumEquation
    reactants,products,Keq=equilibriumEquation
    xMinL=[-molarityMap[prod[0]]*1.0/prod[1] for prod in products]
    xMaxL=[molarityMap[reac[0]]*1.0/reac[1] for reac in reactants]

    xMin=max(xMinL) if len(xMinL) else -1000
    xMax=min(xMaxL) if len(xMaxL) else  1000

    #print "X-RANGE="+str((xMin,xMax))

    X=getShiftToEq(equilibriumEquation,molarityMap,xMin,xMax,maxSerchRounds)

    shift=X*N

    for reac in reactants:
        molarityMap[reac[0]]=molarityMap[reac[0]]-reac[1]*shift
    for prod in products:
        molarityMap[prod[0]]=molarityMap[prod[0]]+prod[1]*shift

    #print "at end of round, M="+str(molarityMap)

def checkEQ(equilibriumEquationList,molarityMap):
    rMax=1.0
    for equation in equilibriumEquationList:
        Qr=Q(equation,molarityMap)
        r=Qr/equation[2]
        if r<1:
            r=1.0/r
        rMax=max(rMax,r)
    print "max error Q from Keq = "+str((rMax-1.0)*100.0)+"%"

def findEquilibriumFinal(equilibriumEquationList,initialMolarityMap,maxApplyRounds=1000,maxSerchRounds=1000,debug=False):
    molarityMap = {}
    for key in initialMolarityMap:
        molarityMap[key]=initialMolarityMap[key]
    
    oldMolarityMap = {}
    N=.5
    for key in initialMolarityMap:
        oldMolarityMap[key]=0
    for i in range(0,maxApplyRounds):
        if(debug):print "ROUND: "+str(i)+"\t"+str(molarityMap)
        #first check no change
        if all(oldMolarityMap[key]==molarityMap[key] for key in initialMolarityMap):
            break
        #else do another round
        for key in initialMolarityMap:
            oldMolarityMap[key]=molarityMap[key]
        for key in initialMolarityMap:
            oldMolarityMap[key]=molarityMap[key]
        #update rConstants
        for equation in equilibriumEquationList: #in the form ((reactants),(products),Keq)
            shiftMap(equation,molarityMap,N,maxSerchRounds)
    return molarityMap

def makeEquation(reactants,products,Keq):
    return (reactants,products,Keq)

def fillInitListWithMissing0(m,Es):
    for e in Es:
        for reac in e[0]:
            if not reac[0] in initialMolarityMap:
                initialMolarityMap[reac[0]]=0
        for prod in e[1]:
            if not prod[0] in initialMolarityMap:
                initialMolarityMap[prod[0]]=0

def anyOverlap(lA,lB):
    for rA in lA:
        if any(rA == rB for rB in lB):
            return True
    return False

def makeCombEq(a,b,fliped):
    r1=a[0]
    p1=a[1]
    r2=b[1] if fliped else b[0]
    p2=b[0] if fliped else b[1]
    r=[]
    r.extend(r1)
    r.extend(r2)
    p=[]
    p.extend(p1)
    p.extend(p2)
    while anyOverlap(r,p):
        for reac in r:
            if reac in p:
                r.remove(reac)
                p.remove(reac)
    res=makeEquation(r,p,a[2]*(1.0/b[2] if fliped else b[2]))
    return res


def genCombEqsN_is_2(es):
    re=[]
    for a,b in itertools.combinations(es, 2):
        if anyOverlap(a[0],b[0]) or anyOverlap(a[1],b[1]):
            re.append(makeCombEq(a,b,True))
        if anyOverlap(a[0],b[1]) or anyOverlap(a[0],b[1]):
            re.append(makeCombEq(a,b,False))
    return re

def genCombEqsN_is_more_than_2(es,setn_m_1):
    re=[]
    for a in es:
        for b in setn_m_1:
            if anyOverlap(a[0],b[0]) or anyOverlap(a[1],b[1]):
                re.append(makeCombEq(a,b,True))
            if anyOverlap(a[0],b[1]) or anyOverlap(a[0],b[1]):
                re.append(makeCombEq(a,b,False))
    return re

def forOrBackSame(a,b):
    return (a[0]==b[0] and a[1]==b[1]) or (a[0]==b[1] and a[1]==b[0])


def removeDoubles(eqs):
    re=[]
    for a in eqs:
        if not any(forOrBackSame(a,r) for r in re):
            re.append(a)
    if re==eqs:
        return list(re)
    else:
        return removeDoubles(re)


def genCombEqs(es):
    last=genCombEqsN_is_2(es)
    re=last
    for n in range(len(es)):
        last=genCombEqsN_is_more_than_2(es,last)
        re.extend(last)

    re = removeDoubles(re)
    return re

def findEquilibrium(equilibriumEquationList,initialMolarityMap,maxApplyRounds=1000,maxSerchRounds=1000,debug=False):
    fillInitListWithMissing0(initialMolarityMap,equilibriumEquationList)
    extraEquations=genCombEqs(equilibriumEquationList)
    equilibriumEquationList.extend(extraEquations)
    mMap=findEquilibriumFinal(equilibriumEquationList,initialMolarityMap,maxApplyRounds,maxSerchRounds,debug)
    print "\n\n"
    print mMap
    print "ph="+str(-np.log10(mMap["H+"]))
    checkEQ(equilibriumEquationList,mMap)

#TODO add text equation parser


e1=makeEquation([],[("H+",1),("OH-",1)],1E-14)
e2=makeEquation([("HC3H5O2",1)],[("C3H5O2-",1),("H+",1)],1.3E-5)
e3=makeEquation([("KOH",1)],[("K+",1),("OH-",1)],1E30)
equilibriumEquationList=[e1,e2,e3]
#add step that combines oposit equations with overlaping parts, i.e. weak acid weak base slow convertion.

initialMolarityMap={"H+":1,"OH-":1,"HC3H5O2":1,"C3H5O2-":0,"KOH":1}
findEquilibrium(equilibriumEquationList,initialMolarityMap)


