import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import exp,log,sqrt
from scipy.stats import norm,gamma

r=[10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 3, 22, 15, 32, 3]
n=[39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 4, 12, 41, 30, 51, 7]
x1=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
x2=[0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
I=21

def calcule_prob(alpha0,alpha1,alpha2,alpha12,tau,d1=x1,d2=x2):
    y=[]
    b=norm.rvs(loc=0, scale=1/tau, size=21)
    for i in range(21):
        z=alpha0 + alpha1 * x1[i] + alpha2 * x2[i] +alpha12 * x1[i] * x2[i] + b[i]
        y.append(exp(z)/(1+exp(z)))
    return y
  

def log_pos(alpha0,alpha1,alpha2,alpha12,tau,a=r,b=n):
    y=0
    prob= calcule_prob(alpha0,alpha1,alpha2,alpha12,tau,x1,x2)
    for i in range(21):
        y=y+r[i]*log(prob[i])+(n[i]-r[i])*log(1.0-prob[i])
    y=y+ norm.logpdf(alpha0, 0 ,1)+norm.logpdf(alpha1, 0 ,  1)+norm.logpdf(alpha2, 0 , 1)+norm.logpdf(alpha12, 0 , 1)+\
    gamma.logpdf(tau, 1, loc=0, scale=1)
    return y
  
def Metropolis_Hastings(N,alpha0,alpha1,alpha2,alpha12,tau):
    Y=np.zeros((N+1,5))
    Y[0,0]=alpha0
    Y[0,1]=alpha1
    Y[0,2]=alpha2
    Y[0,3]=alpha12
    Y[0,4]=tau
    for i in range(N):
        prop_alpha0=Y[i,0]+np.random.randn()*0.45
        top=log_pos(prop_alpha0,Y[i,1],Y[i,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        bottom=log_pos(Y[i,0],Y[i,1],Y[i,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        alpha=min(0,top-bottom)
        u=np.random.rand()
        if log(u)<=alpha:
            Y[i+1,0]=prop_alpha0
        else:
            Y[i+1,0]=Y[i,0]
        prop_alpha1=Y[i,1]+np.random.randn()*0.45
        top=log_pos(Y[i+1,0],prop_alpha1,Y[i,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        bottom=log_pos(Y[i+1,0],Y[i,1],Y[i,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        alpha=min(0,top-bottom)
        u=np.random.rand()
        if log(u)<=alpha:
            Y[i+1,1]=prop_alpha1
        else:
            Y[i+1,1]=Y[i,1]
        prop_alpha2=Y[i,2]+np.random.randn()*0.45
        top=log_pos(Y[i+1,0],Y[i+1,1],prop_alpha2,Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        bottom=log_pos(Y[i+1,0],Y[i+1,1],Y[i,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        alpha=min(0,top-bottom)
        u=np.random.rand()
        if log(u)<=alpha:
            Y[i+1,2]=prop_alpha2
        else:
            Y[i+1,2]=Y[i,2]
        prop_alpha12=Y[i,3]+np.random.randn()*0.45
        top=log_pos(Y[i+1,0],Y[i+1,1],Y[i+1,2], prop_alpha12,Y[i,4],r,n)+log(Y[i,4])
        bottom=log_pos(Y[i+1,0],Y[i+1,1],Y[i+1,2],Y[i,3],Y[i,4],r,n)+log(Y[i,4])
        alpha=min(0,top-bottom)
        u=np.random.rand()
        if log(u)<=alpha:
            Y[i+1,3]=prop_alpha12
        else:
            Y[i+1,3]=Y[i,3]
        prop_tau=exp(log(Y[i,4])+np.random.randn()*0.1)
        top=log_pos(Y[i+1,0],Y[i+1,1],Y[i+1,2], Y[i+1,3], prop_tau,r,n)+log(prop_tau)
        bottom=log_pos(Y[i+1,0],Y[i+1,1],Y[i+1,2],Y[i+1,3],Y[i,4],r,n)+log(Y[i,4])
        alpha=min(0,top-bottom)
        u=np.random.rand()
        if log(u)<=alpha:
            Y[i+1,4]=prop_tau
        else:
            Y[i+1,4]=Y[i,4]
    return Y
  
Y=Metropolis_Hastings(10000,0,0,0,0,10)

y0=Y[:,0]
y0=y0[1000:]
print ("alpha0",np.mean(y0))
plt.plot(np.arange(len(Y)),Y[:,0])

y1=Y[:,1]
y1=y1[1000:]
print ("alpha1",np.mean(y1))
plt.plot(np.arange(len(Y)),Y[:,1])

y2=Y[:,2]
y2=y2[1000:]
print ("alpha2",np.mean(y2))
plt.plot(np.arange(len(Y)),Y[:,2])

y3=Y[:,3]
y3=y3[1000:]
print ("alpha12",np.mean(y3))
plt.plot(np.arange(len(Y)),Y[:,3])

y4=Y[:,4]
y4=y4[1000:]
print ("sigma",1/(np.mean(y4)))
plt.plot(np.arange(len(Y)),Y[:,4])

