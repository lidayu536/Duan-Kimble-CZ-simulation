# Mainly discuss Duan PHYSICAL REVIEW A 67, 032305 ~2003
# equation 19-23
from ast import Delete
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.fftpack import fft
from scipy.fftpack import fftshift
from scipy.integrate import odeint

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


def main(N=int(2**15), Nt=2**20, T=float(10), Omega=0, g=1, kappa=1, Delta=0 , omega_b=10, fun=f, gammae=0.001, gamma=0, somega=0.1, qubitx=1):
    
    tag=int(gammae)
    
    fp0=open(f'./Duan_AnotherMethod/image{tag}/Condition.txt', 'w')
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
    plt.figure(17)
    ftin=fftshift(fft(fomega_in))
    plt.plot(t[int(493*N/1000):int(507*N/1000)], abs(ftin[int(493*N/1000):int(507*N/1000)]), linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/fin(t).png')
    plt.clf()
    plt.figure(1)
    plt.plot(omega, abs(fomega_in))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/abs(f_in(omega)).png')
    plt.clf()
    plt.figure(2)
    plt.plot(omega, np.angle(fomega_in))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/arg(f_in(omega)).png')
    plt.clf()
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
            plt.figure(13)
            plt.clf()
            plt.plot(omega, abs(C_gomegap)/np.sqrt(domega), linewidth=0.3)
            plt.savefig(f'./Duan_AnotherMethod/image{tag}/fomegaseries/abs(fg_{int(i/int(Nt/(2**7)))}(omega)).png', dpi=600)
            plt.clf()
            plt.plot(omega, abs(C_somegap)/np.sqrt(domega), linewidth=0.3)
            plt.savefig(f'./Duan_AnotherMethod/image{tag}/fomegaseries/abs(fs_{int(i/int(Nt/(2**7)))}(omega)).png', dpi=600)
            plt.clf()
            plt.plot(omega, np.angle(C_gomegap), linewidth=0.3)
            plt.savefig(f'./Duan_AnotherMethod/image{tag}/fomegaseries/arg(fg_{int(i/int(Nt/(2**7)))}(omega)).png', dpi=600)
            plt.clf()
            plt.plot(omega, np.angle(C_somegap), linewidth=0.3)
            plt.savefig(f'./Duan_AnotherMethod/image{tag}/fomegaseries/arg(fs_{int(i/int(Nt/(2**7)))}(omega)).png', dpi=600)
            plt.clf()


    
    fgomegap_out=C_gomegap/np.sqrt(domega)
    plt.figure(3)
    plt.plot(omega, abs(fgomegap_out),linewidth=0.1)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/abs(fg_out(omega)).png', dpi=600)
    plt.clf()
    plt.figure(4)
    plt.plot(omega, np.angle(fgomegap_out),linewidth=0.1)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/arg(fg_out(omega)).png', dpi=600)
    plt.clf()
    fsomegap_out=C_somegap/np.sqrt(domega)
    plt.figure(3)
    plt.plot(omega, abs(fsomegap_out),linewidth=0.1)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/abs(fs_out(omega)).png', dpi=600)
    plt.clf()
    plt.figure(4)
    plt.plot(omega, np.angle(fsomegap_out),linewidth=0.1)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/arg(fs_out(omega)).png', dpi=600)
    plt.clf()

    plt.figure(6)
    plt.plot(t_T, abs(C_g_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_g(t).png')
    plt.clf()
    plt.figure(7)
    plt.plot(t_T, abs(C_s_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_s(t).png')
    plt.clf()
    plt.figure(8)
    plt.plot(t_T, abs(C_e_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_e(t).png')
    plt.clf()
    
    plt.figure(14)
    plt.plot(t_T, np.angle(C_g_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_g(t)arg.png')
    plt.clf()
    plt.figure(15)
    plt.plot(t_T, np.angle(C_s_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_s(t)arg.png')
    plt.clf()
    plt.figure(16)
    plt.plot(t_T, np.angle(C_e_series))
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/C_e(t)arg.png')
    plt.clf()
    
    plt.figure(9)
    plt.plot(t_T, abs(C_g_series)**2+abs(C_s_series)**2)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/n_internal(t).png')
    plt.clf()
    plt.figure(10)
    plt.plot(t_T, abs(C_gomegap_series)**2+abs(C_somegap_series)**2)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/n_external(t).png')
    plt.clf()
    plt.figure(11)
    plt.plot(t_T, abs(C_somegap_series)**2+abs(C_gomegap_series)**2+abs(C_g_series)**2+abs(C_s_series)**2)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/n_total(t).png') 
    plt.clf()
    
    fp=open(f'./Duan_AnotherMethod/image{tag}/LossAnalysis.txt', 'w')
    Loss=1-(np.linalg.norm(C_somegap)**2+np.linalg.norm(C_gomegap)**2+abs(C_s)**2+abs(C_g)**2)
    WaveDot=np.dot(C_gomegap, C_gomegap0) - np.dot(C_somegap, C_somegap0)
    Fidelity=abs(WaveDot)**2
    fp.write(f"WaveDot={WaveDot}\n")
    fp.write(f"Loss={Loss}\n")
    fp.write(f"Fidelity={Fidelity}\n")
    fp.close()
    print(Fidelity)

    plt.figure(18)
    fgoutt=fftshift(fft(C_gomegap))/np.sqrt(domega)
    plt.plot(t[int(493*N/1000):int(507*N/1000)], (abs(fgoutt))[int(493*N/1000):int(507*N/1000)], linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/fgout(t).png')
    plt.figure(19)
    fsoutt=fftshift(fft(C_somegap))/np.sqrt(domega)
    plt.plot(t[int(493*N/1000):int(507*N/1000)], (abs(fsoutt))[int(493*N/1000):int(507*N/1000)], linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/image{tag}/fsout(t).png')

def Lf(kappa_ex, kappa_in, g, gammae, Delta, qubitx):
    if(qubitx==0):
        return (-kappa_ex+kappa_in+1j*Delta)/(kappa_ex+kappa_in+1j*Delta)
    else:
        return (-kappa_ex+kappa_in+1j*Delta + g*g/(gammae+1j*Delta) )/(kappa_ex+kappa_in+1j*Delta + g*g/(gammae+1j*Delta) )


def Thr(N=int(2**15), Nt=2**20, T=float(10), Omega=0, g=1, kappa=1, Delta=0 , omega_b=10, fun=f, gammae=0.001, gamma=0, somega=0.1, qubitx=1):
    tag=int(g)
    
    fp=open(f'./Duan_AnotherMethod/theoretical{tag}/Condition_LossAnalysis.txt', 'w')
    fp.write(f"g={g}\n")
    fp.write(f"kappa={kappa}\n")
    fp.write(f"gammae={gammae}\n")
    fp.write(f"gamma={gamma}\n")
    fp.write(f"somega={somega}\n")
    fp.write(f"qubitx={qubitx}\n")

    dT=T/Nt
    t_T=np.linspace(0,T,Nt)
    domega=omega_b*2/N
    omega=np.linspace(-omega_b+domega/2, omega_b-domega/2, N)
    dt=(2*np.pi/Nt)/domega
    t=np.linspace(dt/2,(N-1/2)*dt, N)

    fomega_in=(1+0j)*Normalized2(omega,fun(omega,somega),N)
    plt.figure(17)
    ftin=fftshift(fft(fomega_in))
    plt.plot(t[int(493*N/1000):int(507*N/1000)], abs(ftin[int(493*N/1000):int(507*N/1000)]), linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/fin(t).png')
    plt.clf()
    plt.figure(1)
    plt.plot(omega, abs(fomega_in))
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/abs(f_in(omega)).png')
    plt.clf()
    plt.figure(2)
    plt.plot(omega, np.angle(fomega_in))
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/arg(f_in(omega)).png')
    plt.clf()
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

    C_gomegap=Lf(kappa, 0, g, gammae, omega, 1)*C_gomegap0
    C_somegap=Lf(kappa, 0, g, gammae, omega, 0)*C_somegap0

    fgomegap_out=C_gomegap/np.sqrt(domega)
    plt.figure(3)
    plt.plot(omega, abs(fgomegap_out),linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/abs(fg_out(omega)).png', dpi=600)
    plt.clf()
    plt.figure(4)
    plt.plot(omega, np.angle(fgomegap_out),linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/arg(fg_out(omega)).png', dpi=600)
    plt.clf()
    fsomegap_out=C_somegap/np.sqrt(domega)
    plt.figure(3)
    plt.plot(omega, abs(fsomegap_out),linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/abs(fs_out(omega)).png', dpi=600)
    plt.clf()
    plt.figure(4)
    plt.plot(omega, np.angle(fsomegap_out),linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/arg(fs_out(omega)).png', dpi=600)
    plt.clf()
    plt.figure(18)
    fgoutt=fftshift(fft(C_gomegap))/np.sqrt(domega)
    plt.plot(t[int(493*N/1000):int(507*N/1000)], (abs(fgoutt))[int(493*N/1000):int(507*N/1000)], linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/fgout(t).png')
    plt.clf()
    plt.figure(19)
    fsoutt=fftshift(fft(C_somegap))/np.sqrt(domega)
    plt.plot(t[int(493*N/1000):int(507*N/1000)], (abs(fsoutt))[int(493*N/1000):int(507*N/1000)], linewidth=0.3)
    plt.savefig(f'./Duan_AnotherMethod/theoretical{tag}/fsout(t).png')
    plt.clf()


    Loss=1-(np.linalg.norm(C_somegap)**2+np.linalg.norm(C_gomegap)**2+abs(C_s)**2+abs(C_g)**2)
    WaveDot=np.dot(C_gomegap, C_gomegap0) - np.dot(C_somegap, C_somegap0)
    Fidelity=abs(WaveDot)**2
    fp.write(f"WaveDot={WaveDot}\n")
    fp.write(f"Loss={Loss}\n")
    fp.write(f"Fidelity={Fidelity}\n")
    fp.close()
    print(Fidelity, Loss)
    

#main(g=10, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**15)

#for i in range(10):
#    main(g=10+i, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**17)

#Thr(g=10, T=120, kappa=10, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**15)
for i in range(10):
    Thr(g=1+i, T=120, kappa=100, fun=f1, Omega=0, gammae=1, gamma=0, somega=1, omega_b=10, qubitx=1/2, N=2**15, Nt=2**17)
