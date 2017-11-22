import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from math import sqrt, log

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

if __name__ == "__main__":

    num_elementos = 512
    
    
    x = np.linspace(inicio, fim, num=num_elementos)
    
    plot(x)
