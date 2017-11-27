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

# k = ordem_polinomio
# n = numero de nos da integracao gauss
# num_elementos = num_elementos da malha de integracao
def calc_matriz_local(k,n,num_elementos):
    ke = np.zeros(((k+1),(k+1)))
    pontos_de_gauss = pontos_gauss(n)
    h = 2/num_elementos # (1 - (-1))/num_elementos
    t =  np.linspace(-1, 1, num=num_elementos)
    
    if k == 1:
        for i in range(k+1):
            for j in range(k+1):
                for p in range(k+1): 
                    ke[i][j] = 0
    elif k == 2:
        print("k=2 nao implementado")
        
if __name__ == "__main__":
    num_elementos = 5
    ordem_polinomio = 1

    x = np.linspace(inicio, fim, num=num_elementos)
    i = 0
    #print converte_espaco(x[i],x[i+1])
