#!/usr/bin/env python

from Likelihood import *

import matplotlib.pyplot as plt
#import pylab

def main():

    # Test the pois
    x = np.arange(-5, 5, .1)
    y = [gauss(point, 0, 1) for point in x]
    plt.figure()
    plt.plot(x, y)
    plt.savefig("gauss.pdf")
    
    # Create a likelihood model
    model = Likelihood()
    model.n = 100
    model.mu = 1.0
    model.alpha = 0
    model.delta = 20

    data = [110]
   
    # Test the minimization
    model.minimize(data, params=['mu','alpha'])
                    
    # Plot the likelihood as a function of mu
    x = scipy.linspace(0, 2, num=100)
    y = [model.nll(data, mu=p) for p in x]
    #pylab.plot(x, y)
    plt.figure()
    plt.plot(x,y)
    #plt.axis([0, 3, 0, 20])
    plt.savefig("nll.pdf")    
    return



    #x = scipy.linspace(0,10,11)
    #pylab.plot(x, pois(x,5))

if __name__=="__main__":
    main()
