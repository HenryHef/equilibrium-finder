import random
import numpy as np

COMPLEATION_EXP=100
COMP_MAX=10**COMPLEATION_EXP
COMP_MIN=10**(-COMPLEATION_EXP)


def Q(equilibriumEquation,molarityMap,x):#if dividing by zero, return Q=1E200, no K's above this please
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

def getShiftToEq(e,M,xMin,xMax):
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
    for i in range(0,300): #change to distance from Keq?
        xTest=(xMin+xMax)/2.0
        Qt=Q(e,M,xTest)
        if(Keq==Qt):
            return xTest
        elif(Keq<Qt):
            xMax=xTest
            Qmax=Qt
        elif(Keq>Qt):
            xMin=xTest
            Qmin=Qt
        else:
            print "out of range error"
        #print "step= "+str(i)+"(xMin,xMax)= "+str((xMin,xMax,(xMin+xMax)/2.0))
    return (xMin+xMax)/2.0



def shiftMap(equilibriumEquation,molarityMap,N):
    #print equilibriumEquation
    reactants,products,Keq=equilibriumEquation
    xMinL=[-molarityMap[prod[0]]*1.0/prod[1] for prod in products]
    xMaxL=[molarityMap[reac[0]]*1.0/reac[1] for reac in reactants]

    xMin=max(xMinL) if len(xMinL) else -1000
    xMax=min(xMaxL) if len(xMaxL) else  1000

    #print "X-RANGE="+str((xMin,xMax))

    X=getShiftToEq(equilibriumEquation,molarityMap,xMin,xMax)

    shift=X*N

    for reac in reactants:
        molarityMap[reac[0]]=molarityMap[reac[0]]-reac[1]*shift
    for prod in products:
        molarityMap[prod[0]]=molarityMap[prod[0]]+prod[1]*shift

    #print "at end of round, M="+str(molarityMap)

def findEquilibrium(equilibriumEquationList,initialMolarityMap):
    molarityMap = {}
    for key in initialMolarityMap:
        molarityMap[key]=initialMolarityMap[key]
    
    oldMolarityMap = {}
    N=2.0/3
    rounds=1000
    for key in initialMolarityMap:
        oldMolarityMap[key]=0
    for i in range(0,rounds):
        print "ROUND: "+str(i)+"\t"+str(molarityMap)
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
            shiftMap(equation,molarityMap,N)
    return molarityMap

def makeEquation(reactants,products,Keq):
    return (reactants,products,Keq)

def fillInitListWithMissing0(m,Es):
    for e in Es:
        for reac in e[0]:
            if not reac[0] in initialMolarityMap:
                initialMolarityMap[reac[0]]=0
        for prod in e[0]:
            if not prod[0] in initialMolarityMap:
                initialMolarityMap[prod[0]]=0


e1=makeEquation([],[("H+",1),("OH-",1)],1E-14)
e2=makeEquation([("HC3H5O2",1)],[("C3H5O2-",1),("H+",1)],1.3E-5)
#e3=makeEquation([("KOH",1)],[("K+",1),("OH-",1)],1E30)
equilibriumEquationList=[e1,e2]

for i in range(0,1):
    initialMolarityMap={"H+":1,"OH-":1,"HC3H5O2":1,"C3H5O2-":0}
    fillInitListWithMissing0(initialMolarityMap,equilibriumEquationList)
    mMap=findEquilibrium(equilibriumEquationList,initialMolarityMap)
    print mMap
    print "ph="+str(-np.log10(mMap["H+"]))


