import numpy as np
import matplotlib.pyplot as plt
import math
NA=6.023*1e23

VZr=(((217)**2-(215)**2)+np.pi*((5)**2-(4.5)**2)*196)*1e-6
VMod=((215)**2-196*np.pi*(5)**2)*1e-6
VUO2=(196*np.pi*(4.5)**2)*1e-6

DZr=6800
DUO2=18700
DH2O=1000
DD2O=1000
DNa=820


MZr=DZr*VZr
MH2O=DH2O*VMod
MD2O=DD2O*VMod
MNa=DNa*VMod
MUO2=DUO2*VUO2

MolZr=91*1e-3
MolH2O=18*1e-3
MolD2O=20*1e-3
MolNa=23*1e-3
MolUO2=270.04*1e-3

NZr=MZr*NA/MolZr
NH2O=MH2O*NA/MolH2O
ND2O=MD2O*NA/MolD2O
NNa=MNa*NA/MolNa
NUO2=MUO2*NA/MolUO2


NU235=(NUO2)
NU238=(NUO2)
NO=NUO2*2

muU235=2.42
muU238=2.9

sigfU235=1.02
sigfU238=2.68

sigaU238=2.68
sigaNa=0.53
sigaU235=530
sigaO=0.19
sigaZr=0.18

sigsNa=4.1
sigsO=4.2
sigsU235=7.2
sigsU238=4.3
sigsZr=4.6

zhiNa=2/(23+2/3)


def kinfH20(x):
    epsln=(1+(0.690*(NU238/NH2O)))/(1+(0.563*(NU238/NH2O)))
    print("epsln = ",epsln)

    eta=2.42*NU235*579/(NU235*(101+579)+2.72*NU238)
    print("eta = ",eta)

    f=(NU235*(579+101)+NU238*2.72)/(NU235*(579+101)+NU238*2.72+NO*0+NH2O*0.66+NZr*0.182)
    print("f = ",f)

    #zhi_bar=(NH2O*50*2/(18+2/3)+NU235*10*2/(235+2/3)+NU238*8.3*2/(238+2/3)+NO*3.8*2/(16+2/3)+NZr*8*2/(91+2/3))/(NH2O*50+NU235*10+NU238*8.3+NO*3.8+NZr*8)
    zhi_bar=(NH2O*50*0.92+NU235*10*0.0085+NU238*8.3*.00084+NO*3.8*0.121+NZr*8*0.012)/(NH2O*50+NU235*10+NU238*8.3+NO*3.8+NZr*8)

    p=np.exp(-2.73*(((NH2O*50+NU235*10+NU238*8.3+NO*3.8+NZr*8)/NU238)**-0.514)/zhi_bar)
    print("p = ",p)
    print("kinfiH20 = ",epsln*eta*f*p)

    return epsln*eta*f*p

def kinfD2O(x):
    epsln=(1+(0.690*(NU238/ND2O)))/(1+(0.563*(NU238/ND2O)))
    print("epsln = ",epsln)


    eta=2.42*NU235*579/(NU235*(101+579)+2.72*NU238)
    print("eta = ",eta)

    f=(NU235*(579+101)+NU238*2.72)/(NU235*(579+101)+NU238*2.72+NO*0+ND2O*0.001+NZr*0.182)
    print("f = ",f)

    zhi_bar=(ND2O*10.6*0.509+NU235*10*0.0085+NU238*8.3*.00084+NO*3.8*0.121+NZr*8*0.012)/(ND2O*10.6+NU235*10+NU238*8.3+NO*3.8+NZr*8)

    p=np.exp(-2.73*(((ND2O*10.6+NU235*10+NU238*8.3+NO*3.8+NZr*8)/NU238)**-0.514)/zhi_bar)
    print("p = ",p)
    print("kinfiH20 = ",epsln*eta*f*p)

    return epsln*eta*f*p

def kinfNa(x):
    epsln=(1+(0.690*(NU238/NNa)))/(1+(0.563*(NU238/NNa)))
    print("epsln = ",epsln)

    eta=(muU235*sigfU235*NU235+muU238*sigfU238*NU238)/(NU235*sigaU235+NU238*sigaU238)
    print("eta = ",eta)

    f=(NU235*sigaU235+NU238*sigaU238)/(NU235*sigaU235+NU238*sigaU238+NO*sigaO+NNa*sigaNa+NZr*sigaZr)
    print("f = ",f)

    #zhi_bar=(NH2O*50*2/(18+2/3)+NU235*10*2/(235+2/3)+NU238*8.3*2/(238+2/3)+NO*3.8*2/(16+2/3)+NZr*8*2/(91+2/3))/(NH2O*50+NU235*10+NU238*8.3+NO*3.8+NZr*8)
    zhi_bar=(NNa*sigsNa*zhiNa+NU235*sigsU235*0.0085+NU238*sigsU238*.0084+NO*sigsO*0.121+NZr*sigsZr*0.0218)/(NNa*sigsNa+NU235*sigsU235+NU238*sigsU238+NO*sigsO+NZr*sigsZr)

    p=np.exp(-2.73*(((NNa*sigsNa+NU235*sigsU235+NU238*8.3+NO*3.8+NZr*8)/NU238)**-0.514)/zhi_bar)
    p=1
    print(epsln*eta*f*p)
    return epsln*eta*f*p

while True :
    #print("Data for H2O as moderator")
    print("----------------------------------------------------------------------------------------------------------------------------------------------")
    print("Data for Na as moderator")
    x=float(input("x = "))
    NU235=(NUO2)*(x/(x+1))
    NU238=(NUO2)/(x+1)
    kinfNa(x)
