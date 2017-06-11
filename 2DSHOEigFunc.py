#Program that plots the eigen functions of a 2D simple harmonic oscillator, with angular quantum number set to zero.
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

N=60000 # Number of iterations of the Numerov Algorithm
h = 0.0001 # Step size in the algorithm
h2 = pow(h,2) 


def main():

  epsilon = getEigenState() #Prompts the user for the energy eigenstate to be plotted

  x_out, z_out = NumerovAlg(N, h2, epsilon) 

  p = np.linspace(0, 2*np.pi, 100)  #1D array of angles in the polar coordinate system, between 0 and 2*pi

  # print (x_out)

  R, P = np.meshgrid(x_out, p) 
  #Creates 2D arrays R and P that hold the polar coordinates r and $\theta$ respectively. 

  X, Y = R*np.cos(P), R*np.sin(P) 

  Z = np.array([z_out for i in range(len(P))])
  # Z = np.outer(z_out, z_out)
  print(Z)
  print(X.shape,Y.shape, Z.shape)

  surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)

  # Customize the z axis.
  ax.text2D(0.05, 0.95, "Eigenfunctions with angular quantum number zero for radial quantum number " + str(int(epsilon - 0.5)), transform=ax.transAxes)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.set_xlabel('X axis')
  ax.set_ylabel('Y axis')
  ax.set_zlabel('$\Psi$')
  # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

  # Add a color bar which maps values to colors.
  fig.colorbar(surf, shrink=0.5, aspect=5)

  plt.show()


# Returns list x_out that is the radial axis, and z_out that is the value of the wavefunction on each axis
def NumerovAlg(N, h2, epsilon):
  z = 0.0
  k = 0.0
  x = -1*(N-2)*h

  k_minus_2 = epsilon + x-2*h # k_0
  k_minus_1 = epsilon + x-h # k_1
  a = 0.1
  z_minus_2 = 0 # z_0
  z_minus_1 = a # z_1

  x_out = []
  z_out = []

  n=-1*N+2
  while n < N-2:
    n+=1
    x += h;
    k = 2*epsilon - pow(x, 2)
    b = h2/12
    z = ( 2*(1-5*b*k_minus_1) * z_minus_1 - (1+b*k_minus_2) * z_minus_2 ) / (1 + b * k)

    # Save for plotting. Add new values of x and z to the lists x_out and z_out
    if (n % 10 == 0):
      x_out.append(x)
      z_out.append(z)

    # Shift for next iteration
    z_minus_2 = z_minus_1
    z_minus_1 = z
    k_minus_2 = k_minus_1
    k_minus_1 = k

  return x_out, z_out #The x_out corresponds to the radial axis since the potential well is radially symmetrical

def getEigenState():
  epsilon = 0
  while True:
    inp = input("Enter the radial quantum number of the eigenstate you want to plot(even integers only): ")
    
    if inp.isdigit() and int(inp) % 2 == 0:
      epsilon = (0.5 + float(inp))
      break
    else:
      print("That's not a positive integer!")

  return epsilon

main()
