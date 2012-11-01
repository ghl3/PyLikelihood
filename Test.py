#!/usr/bin/env python


from ConfidenceInterval import *

import math

def create_model():

    # Annoying boilerplate
    def pois(x, mu):
        return poisson.pmf(math.floor(x), mu)
    def gauss(x, mu, sigma):
        return norm.pdf(x, mu, sigma)

    d = variable("d", 0, 20, 100)
    s = variable("s", 0, 20, 100)
    b = variable("b", 0, 20, 100)

    def pdf(d, s, b, b0=5.0, sigma=1.0):
        #return gauss(d, mu, 2.0)*gauss(mu0, mu, sigma)
        return pois(d, s+b)*gauss(b0, b, sigma)

    # Create the likelihood
    model = likelihood(pdf, data=d, params=[s, b])
    model.logging.setLevel(logging.DEBUG)

    return model


def print_model(model):

    # Print the model
    print model
    from pprint import pprint
    pprint (vars(model))

    # Test the setting of parameters
    model.s = 1.0
    model.b = 5.0

    # Testing the normalization on evaluation
    print model.norm
    model.eval(d=5)
    print model.norm


def test_minimization(model, obs_data):

    # Test the minimization
    print "Fitting to data=10: "
    #obs_data = 10
    model.fitTo(obs_data, params=["s", "b"])
    nll_min = model.nll(obs_data)


def test_profile(model, obs_data):

    plt.clf()

    # Testing the profile
    model.fitTo(obs_data, params=["s", "b"])
    nll_min = model.nll(obs_data)
    profile_min = model.profile(obs_data, "s", nuisance=["b"])
    print "nll_min: ", nll_min
    print " profile min: ", profile_min

    # create the x values
    x = model.data.linspace()

    # Draw Nll
    model.fitTo(obs_data, params=["s", "b"])
    y = [model.nll(obs_data, s=d) - nll_min for d in x]
    plt.plot(x,y)

    # Draw profile nll
    model.fitTo(obs_data, params=["s", "b"])
    z = [model.profile(obs_data, poi="s", s=d, nuisance=["b"], ) for d in x]
    plt.plot(x,z)

    (ymin, ymax) = plt.ylim()
    (xmin, xmax) = plt.xlim()
    plt.ylim([0, ymax]) 
    plt.hlines(0.5, xmin, xmax)

    plt.xlabel('mu')
    plt.ylabel('profile likelihood(x)')
    plt.savefig("profile.pdf")

    plt.clf()


def test_interval(model, obs_data):

    # Test the interval functionality (over data)
    plt.clf()
    model.fitTo(obs_data, params=["s", "b"])
    model.make_plot(interval=0.68)
    plt.xlabel('data')
    plt.savefig("plot.pdf")

    # clear the current figure and
    # plot the neyman interval
    plt.clf()
    neyman = model.get_neyman(0.68, "s")
    for pair in neyman:
        (s, x0, x1) = pair[0], pair[1][0], pair[1][1]
        plt.hlines(s, x0, x1)
    plt.xlabel('x')
    plt.ylabel('s')
    plt.savefig("neyman.pdf")

    #print "inverted Neyman 4: ", model.invert_neyman(4, neyman)
    print "inverted Neyman 8: ",  model.invert_neyman(8, neyman)
    print "inverted Neyman 10: ", model.invert_neyman(10, neyman)
    print "inverted Neyman 12: ", model.invert_neyman(12, neyman)
    print "inverted Neyman 14: ", model.invert_neyman(14, neyman)


def test_mc(model, obs_data):

    model.fitTo(params=["s", "b"], d=obs_data)
    samples = model.sample_mc(['d'], 1000)
    values = [point['d'] for point in samples]
    print values
    n, bins, patches = plt.hist(values, 20, normed=1, facecolor='g')
    plt.xlabel('data')
    plt.ylabel('Probability')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    plt.savefig("sampled_data.pdf")

if __name__=="__main__":
    model = create_model()
    print_model(model)

    obs_data = 7

    #test_minimization(model, obs_data)
    #test_profile(model, obs_data)
    #test_interval(model, obs_data)
    test_mc(model, obs_data)
