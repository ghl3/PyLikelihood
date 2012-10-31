#!/usr/bin/env python


from ConfidenceInterval import *


def main():

    # Annoying boilerplate
    def pois(x, mu):
        return poisson.pmf(x, mu)
    def gauss(x, mu, sigma):
        return norm.pdf(x, mu, sigma)

    d = variable("d", 0, 10, 100)
    mu = variable("mu", 0, 10, 100)
    mu0 = variable("mu0", 0, 10, 100)

    def pdf(d, mu=5.0, mu0=5.0, sigma=2.0):
        return gauss(d, mu, 2.0)*gauss(mu0, mu, sigma)


    #func_spec = inspect.getargspec(pdf)
    #(all_arguments, all_defaults) = (func_spec.args, func_spec.defaults)
    #print "Spec: ", func_spec
    #print "Arguments: ", all_arguments
    #print "Defaults: ", all_defaults

    # Create the likelihood
    model = likelihood(pdf, data=d, params=[mu, mu0])

    print model

    from pprint import pprint
    pprint (vars(model))
    print model.mu
    print model.mu0

    model.mu = 5.0
    model.mu0 = 5.0
    #model.sigma = 1.0

    print model.norm
    model.eval(5)
    print model.norm

    return



    for point in range(0, 10):
        print point, ": ", model(point), ": ", model.eval(point)

    print model.get_interval(0.68)

    model.make_plot(interval=0.68)
    plt.savefig("plot.pdf")

    # clear the current figure and
    # plot the neyman interval
    plt.clf()
    neyman = model.get_neyman(0.68, "mu")
    for pair in neyman:
        #print pair
        (mu, x0, x1) = pair[0], pair[1][0], pair[1][1]
        plt.hlines(mu, x0, x1)
    plt.xlabel('x')
    plt.ylabel('mu')
    plt.savefig("neyman.pdf")

    print "inverted Neyman 4.5: ", model.invert_neyman(4.5, neyman)
    print "inverted Neyman 5.0: ", model.invert_neyman(5.0, neyman)
    print "inverted Neyman 5.5: ", model.invert_neyman(5.5, neyman)
    print "inverted Neyman 6.0: ", model.invert_neyman(6.0, neyman)
    print "inverted Neyman 8.0: ", model.invert_neyman(8.0, neyman)


if __name__=="__main__":
    main()
