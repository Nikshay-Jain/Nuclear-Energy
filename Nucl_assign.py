import numpy as np
import matplotlib.pyplot as plt

def kh2o(x):
    #x = float(input("Enter x if U235 : U238 = 1:x\n"))
    n235 = 5.2011*1e26/(x+1)
    n238 = 5.2011*1e26*x/(x+1)
    nh2o = 1.032*1e27
    nzr = 1.699*1e26
    no = 1.04022*1e27

    eta = (2.42*579*n235)/(n235*680 + n238*2.72)
    f = (n235*680 + n238*2.72)/(n235*680 + n238*2.72 + nh2o*0.66 + nzr*0.182)
    eps = (1+(0.690*n238/nh2o))/(1+(0.563*n238/nh2o))
    zi_bar = (n235*10*0.0085 + n238*8.3*0.0084 + nh2o*50*0.92 + nzr*8*0.02182 + no*0.121*3.8)/(n235*10 + n238*8.3 + nh2o*50 + nzr*8 + no*3.8)
    sig = n235*10 + n238*8.3 + nh2o*50 + nzr*8 + no*3.8
    p = np.exp((-2.73/zi_bar)*((sig/n238)**-0.514))
    k = eps*p*f*eta
    print(eps,eta,f,p,k)
    return (k)

def kd2o(x):
    #x = float(input("Enter x if U235 : U238 = 1:x\n"))
    n235 = 5.2011*1e26/(x+1)
    n238 = 5.2011*1e26*x/(x+1)
    nd2o = 9.288*1e26
    nzr = 1.699*1e26
    no = 1.04022*1e27

    eta = (2.42*579*n235)/(n235*680 + n238*2.72)
    f = (n235*680 + n238*2.72)/(n235*680 + n238*2.72 + nd2o*0.001 + nzr*0.182)
    eps = (1+(0.690*n238/nd2o))/(1+(0.563*n238/nd2o))
    zi_bar = (n235*10*0.0085 + n238*8.3*0.0084 + nd2o*10.6*0.509 + nzr*8*0.02182 + no*3.8*0.121)/(n235*10 + n238*8.3 + nd2o*10.6 + nzr*8 + no*3.8)
    sig = (n235*10 + n238*8.3 + nd2o*10.6 + nzr*8 + no*3.8)
    p = np.exp((-2.73/zi_bar)*((sig/n238)**-0.514))
    k = eps*p*f*eta
    print(eps,eta,f,p,k)
    return k

def kna(x):
    #x = float(input("Enter x if U235 : U238 = 1:x\n"))
    sigf235 = 1.02
    sigs235 = 7.2
    siga235 = 530

    sigf238 = 2.68
    sigs238 = 4.3
    siga238 = 2.68

    sigazr = 0.18
    sigszr = 4.6

    sigao = 0.19
    sigso = 4.2

    sigana = 0.53
    sigsna = 4.1

    n235 = 5.2011*1e26/(x+1)
    n238 = 5.2011*1e26*x/(x+1)
    nna = 0.03084*820*6.023*1e23/0.023
    nzr = 1.699*1e26
    no = 1.04022*1e27

    eta = (2.42*sigf235*n235 + 2.9*sigf238*n238)/(n235*(siga235) + n238*(siga238))
    f = (n235*(siga235) + n238*(siga238))/(n235*(siga235) + n238*(siga238) + no*sigao + nna*sigana + nzr*sigazr)
    eps = (1+(0.690*n238/nna))/(1+(0.563*n238/nna))
    # zi_bar = (n235*sigs235*0.0085 + n238*sigs238*0.0084 + nna*sigsna*0.8 + nzr*sigszr*0.02182 + no*sigso*0.121)/(n235*sigs235 + n238*sigs238 + nna*sigsna + nzr*sigszr + no*sigso)
    # sig = (n235*sigs235 + n238*sigs238 + nna*sigsna + nzr*sigszr + no*sigso)
    # p = np.exp((-2.73/zi_bar)*((sig/n238)**-0.514))
    k = eps*f*eta
    print(eps,eta,f,k)
    return k

x = np.arange(0.001,200,0.001)
for i in x:
    if abs(kh2o(i)-1.3) <= 0.001:
        plt.plot(i,kh2o(i), marker = 'o')
        plt.text(i, kh2o(i), '({}, {})'.format(i, kh2o(i)))
        break

for i in x:
    if abs(kd2o(i)-1.2) <= 0.001:
        plt.plot(i,kd2o(i), marker = 'o')
        plt.text(i, kd2o(i), '({}, {})'.format(i, kd2o(i)))
        break

for i in x:
    if abs(kna(i)-1.25) <= 0.001:
        plt.plot(i,kna(i), marker = 'o')
        plt.text(i, kna(i), '({}, {})'.format(i, kna(i)))
        break

plt.xlabel("x")
plt.ylabel("k")
plt.plot(x,kh2o(x))
plt.plot(x,kd2o(x))
plt.plot(x,kna(x))
plt.show()

# x = float(input("Enter x if U235 : U238 = 1:x\n"))
# kna(x)