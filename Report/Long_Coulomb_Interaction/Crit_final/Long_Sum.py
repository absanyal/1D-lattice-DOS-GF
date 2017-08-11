# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 17:28:36 2017

@author: AB Sanyal
"""

#imports
import numpy as np
import matplotlib.pyplot as plt

#Hopping term
t = 1

#Two interaction terms
U1 = 7.8 * t
U2 = 0.0 * t

#Epsilon
epsilon = 0.1

#Cutoff m_c; G_(m_c + 1) = 0
m_c = 100

#f = 2 t cos ( k / 2)
def f(k):
    return ( 2 * t * np.cos(k / 2) )

#Kronecker-Delta function
def kron_delta(i, j):
    if ( int(i) == int(j) ):
        return 1
    else:
        return 0

#Potential

alpha = 0.0
beta = 1.0

def U(m):
    return ( U1 * np.exp( - ( alpha ) * m * m ) / ( pow( m, beta ) ) )

#z
def z(w):
    return ( complex( w, epsilon ) )

#A matrix
def A(m):
    if ( m == m_c + 1 ):
        return ( -f(k) / ( (z(w) - U(m)) * 1/(z(w) - U(m)) * f(k) * f(k) ) )
    else:
        return ( -f(k) / ( (z(w) - U(m) ) + f(k) * A(m + 1) ) )

#Arrays for storage
w_array = []
ImG1_array = []

#Segment for looping over k values

k = 0 #Initial value of k
k_stop = np.pi #Stopping value of k
c = 0 #Artificial offset

A_sum = []
A_sum = np.array(A_sum)

while (k <= k_stop): #Loop over k values
    w_array = []
    ImG1_array = []

    w = -5 #Start value of w

    while (w <= 15): #Loop over w
        G1 = 1 / ( z(w) - U(1) + f(k) * A(2) ) #Actual calculations are all happening here
        w_array.append(w)
        ImG1_array.append( (-1 / np.pi) * G1.imag)
        w += 0.01 #Increment on w

    if ( round(k, 4) == 0 ): #Colors the the Brilloin zone edges
        col = 'blue'
    elif ( round(k, 4) == round(k_stop, 4) ):
        col = 'red'
    else:
        col = 'black'

    if (k == 0):
        A_sum = ImG1_array[:]
    else:
        A_sum += np.array(ImG1_array)

    #plt.plot( w_array, ImG1_array, col ) #Plot for a particular k
    #plt.title( "U1 = " + str(U1) + " t, U2 = " + str(U2) + " t, alpha = " + str(alpha) + ", beta = " + str(beta) )

    print( "Finished for k =", k ) #Just to keep track of where in the loop we are

    #Incrementing counters
    k += ( np.pi / 20 )
    c += 0.05

A_sum = list(A_sum / 20)

plt.title( r"$U =\,$" + str(U1) + r"$\,t\,$" + "," + r"$\,\alpha\,=\,$" + str(alpha) + "," + r"$\,\beta\,=\,$" + str(beta) )
plt.plot(w_array, A_sum)
plt.ylim([0, 0.2])

plt.ylabel(r'$\rho \left( \omega \right )$')
plt.xlabel(r'$\omega$')
savename = "alpha_" + str(int(alpha*1000)) +"_beta_" + str(int(beta*1000)) + "_U_" + str(int(U1*1000)) + ".pdf"
#plt.show()
plt.savefig(savename, bbox_inches = 'tight')