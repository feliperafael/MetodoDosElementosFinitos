import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

inicio = 0
fim = 1

def exata(x):
    return np.sin(2*np.pi*x)

def plot(x):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlabel('x')
    plt.ylabel('u')
    ax.plot(x,exata(x),'r')
    plt.show()

def solucao_aproximada(x):
    print("a")

# t = ponto de gauss
def peso_gauss(n,t):
    if n == 2:
        return 1.0
    elif n == 3:
        if t == 0:
            return 8.0/9.0
        else:
            return 5.0/9.0
    print "ERROR: invalid value of n"
    return -1

def pontos_gauss(n):
    pontos = np.zeros(n)
    if n == 2:
        pontos[0] = -(np.sqrt(3.0)/3.0)
        pontos[1] = -pontos[0]
    if n == 3:
        pontos[0] = 0
        pontos[1] = -np.sqrt(3.0/5.0)
        pontos[2] = -pontos[1]    
    return pontos

#def converte_espaco(x1,x2):
#    t =  np.linspace(-1, 1, num=num_elementos)
#    
#    F = (t*((x2-x1)/2.0) + (x2+x1)/2.0)* (x2+x1/2.0)
#    return F
    
def quadratura_de_gauss(n):
    print("quadratura de gauss")

# onde grau eh o grau do polinomio
# i e qual funcao phi
# t e o valor para avaliar
def phi_function(grau,i,t):
    if grau == 1: #linear
        if i == 1 :
            return 0.5*(1-t)
        if i == 2 :
            return 0.5*(1+t)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    if grau == 2: #quadratico
        if i == 1 :
            return 0.5*t*(t-1)
        if i == 2 :
            return -1*(t-1)*(t+1)
        if i == 3:
            return 0.5*t*(t+1)
        else:
            print("PHI FUNCTION RETURN ZERO")
            return 0
    print("ERROR PHI FUNCTION")
       
# k = ordem_polinomio
# n = numero de nos da integracao gauss
# num_elementos = num_elementos da malha de integracao
def calc_matriz_local(k,n,num_elementos):
    dimensao_matriz_global = k*num_elementos + 1
    ke = np.zeros(((k+1),(k+1)))
    K = np.zeros(((dimensao_matriz_global),(dimensao_matriz_global)))
    pontos_de_gauss = pontos_gauss(n)
    h = 2/num_elementos # (1 - (-1))/num_elementos
    t =  np.linspace(-1, 1, num=num_elementos)
    numero_pontos_de_gauss = 2
    if k == 1:
        for i in range(k+1):
            for j in range(k+1):
                for phi_de_gauss in range(k+1):
                    t_pontos_de_gauss = pontos_gauss(numero_pontos_de_gauss)
                    #print t_pontos_de_gauss
                    #print phi_function(k,phi_de_gauss+1,t_pontos_de_gauss)
                    ke[i][j] += phi_function(k,phi_de_gauss+1,t_pontos_de_gauss)
                    print ke
    elif k == 2:
        print("k=2 nao implementado")
    print ke
if __name__ == "__main__":
    num_elementos = 5
    k = 1 #ordem_polinomio = 1

    x = np.linspace(inicio, fim, num=num_elementos)
    i = 0
    calc_matriz_local(k,k+1,num_elementos)
