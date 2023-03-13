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
    n235 = 5.2011*1e26/(x+1)
    n238 = 5.2011*1e26*x/(x+1)
    nna = 0.03084*820*6.023*1e23/0.023
    nzr = 1.699*1e26
    no = 1.04022*1e27

    eta = (2.42*2.3*n235 + 2.9*2.87*n238)/(n235*3.6 + n238*3.47)
    f = (n235*680 + n238*2.72)/(n235*680 + n238*2.72 + nna*0.001 + nzr*0.182)
    eps = (1+(0.690*n238/nna))/(1+(0.563*n238/nna))
    zi_bar = (n235*10*0.0085 + n238*8.3*0.0084 + nna*10.6*0.509 + nzr*8*(2/(91.33+0.667)) + no*3.8*0.121)/(n235*10 + n238*8.3 + nna*10.6 + nzr*8 + no*3.8)
    sig = (n235*10 + n238*8.3 + nna*10.6 + nzr*8 + no*3.8)
    p = np.exp((-2.73/zi_bar)*((sig/n238)**-0.514))
    print(sig,eta,f,eps,p,zi_bar)
    k = eps*p*f*eta
    return k

# x = np.arange(0.001,150,0.001)
# for i in x:
#     if kh2o(i)-1.3<=0.0001:
#         plt.plot(i,kh2o(i), marker = 'o')
#         plt.text(i, kh2o(i), '({}, {})'.format(i, kh2o(i)))
#         break

# for i in x:
#     if kd2o(i)-1.2<=0.0001:
#         plt.plot(i,kd2o(i), marker = 'o')
#         plt.text(i, kd2o(i), '({}, {})'.format(i, kd2o(i)))
#         break

# for i in x:
#     if kna(i)-1.3<=0.0001:
#         plt.plot(i,kna(i), marker = 'o')
#         plt.text(i, kna(i), '({}, {})'.format(i, kna(i)))
#         break

# plt.xlabel("x")
# plt.ylabel("k")
# plt.plot(x,kh2o(x))
# plt.plot(x,kd2o(x))
# plt.plot(x,kna(x))
# plt.show()
x = float(input("Enter x if U235 : U238 = 1:x\n"))
kh2o(x)