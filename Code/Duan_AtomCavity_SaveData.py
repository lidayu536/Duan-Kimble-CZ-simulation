# Mainly discuss Duan PHYSICAL REVIEW A 67, 032305 ~2003
# equation 19-23
from ast import Delete
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.fftpack import fft
from scipy.fftpack import fftshift
from scipy.integrate import odeint
import os

def f(t,T):
    return np.exp(-((t-T/2)/(T/10))**2)

def f1(omega,s):
    return np.exp(-(omega/s)**2)

def Normalized2(x,f,N):
    I=float(0)
    for i in range(N-1):
        I+=(abs((f[i]+f[i+1])/2)**2)*(x[i+1]-x[i])
    return f/np.sqrt(I)

def Normalized1(x,f,N):
    I=float(0)
    for i in range(N-1):
        I+=(abs((f[i]+f[i+1])/2))*(x[i+1]-x[i])
    return f/I

def Update_dt(t_T, dt, omega, C_e0, C_g0, C_gomegap0, C_somegap0, C_s0, g, kappap, gammae, gamma, Delta): 
    #to calculate the 'C's the next dt 
    C_e_t= C_e0*(-1j*Delta-gammae/2) - 1j*g*C_g0
    C_g_t= C_g0*(-1j*Delta-gamma/2) - 1j*g*C_e0 + kappap*np.sum(C_gomegap0*np.exp(-1j*omega*t_T))
    C_s_t= C_s0*(-1j*Delta-gamma/2) + kappap*np.sum(C_somegap0*np.exp(-1j*omega*t_T))
    C_gomegap_t= - kappap*C_g0*np.exp(1j*omega*t_T)
    C_somegap_t= - kappap*C_s0*np.exp(1j*omega*t_T)

    C_e1=C_e0+C_e_t*dt
    C_g1=C_g0+C_g_t*dt
    C_gomegap1=C_gomegap0+C_gomegap_t*dt
    C_somegap1=C_somegap0+C_somegap_t*dt
    C_s1=C_s0+C_s_t*dt

    #Improved Euler
    for k in range(2):
        C_e_t= (C_e0+C_e1)/2*(-1j*Delta-gammae/2) - 1j*g*(C_g0+C_g1)/2
        C_g_t= (C_g0+C_g1)/2*(-1j*Delta-gamma/2) - 1j*g*(C_e0+C_e1)/2 + kappap*(np.sum(C_gomegap0*np.exp(-1j*omega*t_T))+np.sum(C_gomegap1*np.exp(-1j*omega*(t_T+dt))))/2
        C_s_t= (C_s0+C_s1)/2*(-1j*Delta-gamma/2) + kappap*(np.sum(C_somegap0*np.exp(-1j*omega*t_T))+np.sum(C_somegap1*np.exp(-1j*omega*(t_T+dt))))/2
        C_gomegap_t= - kappap*(C_g0*np.exp(1j*omega*t_T)+C_g1*np.exp(1j*omega*(t_T+dt)))/2
        C_somegap_t= - kappap*(C_s0*np.exp(1j*omega*t_T)+C_s1*np.exp(1j*omega*(t_T+dt)))/2
        
        C_e1=C_e0+C_e_t*dt
        C_g1=C_g0+C_g_t*dt
        C_gomegap1=C_gomegap0+C_gomegap_t*dt
        C_somegap1=C_somegap0+C_somegap_t*dt
        C_s1=C_s0+C_s_t*dt

    return [C_e1, C_g1, C_gomegap1, C_somegap1, C_s1]


def sum_of_squares(lst):
    if len(lst)==0:
        return 0
    return sum([i**2 for i in lst])

def SaveData(path, name, data):
    N=data.size
    fp=open(path+f"/{name}.txt", "w")
    for i in range(N):
        fp.write(f"{data[i]}\n")
    fp.close()

def main(N=int(2**15), Nt=2**20, T=float(10), Omega=0, g=1, kappa=1, Delta=0 , omega_b=10, fun=f, gammae=0.001, gamma=0, somega=0.1, qubitx=1):

    path0 = os.getcwd()
    path = path0+f"/Duan_SaveData/g={g},kappa={kappa},gammae={gammae},gamma={gamma},somega={somega},qubitx={qubitx}"
    path1 = path+f"/fomega_series"
    if not os.path.exists(path1):
        os.makedirs(path1)

    fp0=open(path + f'/Condition.txt', 'w')
    fp0.write(f"g={g}\n")
    fp0.write(f"kappa={kappa}\n")
    fp0.write(f"gammae={gammae}\n")
    fp0.write(f"gamma={gamma}\n")
    fp0.write(f"somega={somega}\n")
    fp0.write(f"qubitx={qubitx}\n")
    fp0.close()
    
    dT=T/Nt
    t_T=np.linspace(0,T,Nt)
    domega=omega_b*2/N
    omega=np.linspace(-omega_b+domega/2, omega_b-domega/2, N)
    dt=(2*np.pi/Nt)/domega
    t=np.linspace(dt/2,(N-1/2)*dt, N)
    
    fomega_in=(1+0j)*Normalized2(omega,fun(omega,somega),N)
    
    SaveData(path, "omega", omega)
    SaveData(path, "f_in(omega)",fomega_in)

    C_g=(0+0j)
    C_s=(0+0j)
    C_e=(0+0j)
    C_gomegap0=(np.sqrt(qubitx)+0j)*fomega_in*np.sqrt(domega)
    C_somegap0=(np.sqrt(1-qubitx)+0j)*fomega_in*np.sqrt(domega)
    C_gomegap=(np.sqrt(qubitx)+0j)*fomega_in*np.sqrt(domega)
    C_somegap=(np.sqrt(1-qubitx)+0j)*fomega_in*np.sqrt(domega)
    kappap=np.sqrt(kappa*domega/2*np.pi)
    C_g_series=np.linspace(0+0j,0+0j,Nt)
    C_s_series=np.linspace(0+0j,0+0j,Nt)
    C_e_series=np.linspace(0+0j,0+0j,Nt)
    C_gomegap_series=np.linspace(0+0j,0+0j,Nt)
    C_somegap_series=np.linspace(0+0j,0+0j,Nt)
    for i in range(Nt):
        [C_e, C_g, C_gomegap, C_somegap, C_s]=Update_dt(t_T[i]-T/4, dT, omega, C_e, C_g, C_gomegap, C_somegap, C_s, g, kappap, gammae, gamma, Delta)
        C_g_series[i]=C_g
        C_s_series[i]=C_s
        C_e_series[i]=C_e
        C_gomegap_series[i]=np.linalg.norm(C_gomegap)
        C_somegap_series[i]=np.linalg.norm(C_somegap)
        if i%int(Nt/(2**7))==0 :
            SaveData(path1, f"fg(omega)|t={t_T[i]}", C_gomegap/np.sqrt(domega))
            SaveData(path1, f"fs(omega)|t={t_T[i]}", C_somegap/np.sqrt(domega))


    
    fgomegap_out=C_gomegap/np.sqrt(domega)
    fsomegap_out=C_somegap/np.sqrt(domega)
    SaveData(path, "fg_out(omega)", fgomegap_out)
    SaveData(path, "fs_out(omega)", fsomegap_out)

    SaveData(path, "t_series", t_T)
    SaveData(path, "C_g(t)", C_g_series)
    SaveData(path, "C_s(t)", C_s_series)
    SaveData(path, "C_e(t)", C_e_series)
    
    SaveData(path, "n_internal(t)", abs(C_g_series)**2+abs(C_s_series)**2)
    SaveData(path, "n_external(t)", abs(C_gomegap_series)**2+abs(C_somegap_series)**2)
    SaveData(path, "n_total(t)", abs(C_somegap_series)**2+abs(C_gomegap_series)**2+abs(C_g_series)**2+abs(C_s_series)**2)
    
    fp=open(path+f'/LossAnalysis.txt', 'w')
    Loss=1-(np.linalg.norm(C_somegap)**2+np.linalg.norm(C_gomegap)**2+abs(C_s)**2+abs(C_g)**2)
    WaveDot=np.dot(C_gomegap, C_gomegap0) - np.dot(C_somegap, C_somegap0)
    Fidelity=abs(WaveDot)**2
    fp.write(f"WaveDot={WaveDot}\n")
    fp.write(f"Loss={Loss}\n")
    fp.write(f"Fidelity={Fidelity}\n")
    fp.close()
    print(Fidelity, Loss)

    fp=open(path+"/Fidelity.txt", "w")
    fp.write(Fidelity+'\n')
    fp.close()
    fp=open(path+"/Loss.txt", "w")
    fp.write(Loss+'\n')
    fp.close()



main(g=10, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**15)

#for i in range(10):
#    main(g=10+i, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**17)

#Thr(g=10, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**15)
