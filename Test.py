#!/usr/bin/env python


from ConfidenceInterval import *

import math


def create_model():

    print "Creating Model"

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

    # Test the setting of parameters
    model.s = 1.0
    model.b = 5.0

    return model


def print_model(model):

    print "Printing Model"

    # Print the model
    print model
    from pprint import pprint
    pprint (vars(model))

    # Testing the normalization on evaluation
    print model.norm
    model.eval(d=5)
    print model.norm


def test_minimization(model, obs_data):

    print "Test Minimization"

    # Test the minimization
    print "Fitting to data=%s: " % obs_data
    #obs_data = 10
    model.fitTo(obs_data, params=["s", "b"])
    nll_min = model.nll(d=obs_data)
    print "Fitted State: ", model.total_state()
    print "Fitted Nll: ", nll_min


def test_profile(model, obs_data):

    print "Test Profile"

    plt.clf()

    # Testing the profile
    model.fitTo(obs_data, params=["s", "b"])
    nll_min = model.nll(d=obs_data)
    profile_min = model.profile(d=obs_data, poi="s", nuisance=["b"])
    print "nll_min: ", nll_min
    print " profile min: ", profile_min

    # create the x values
    x = model.data.linspace()

    # Draw Nll
    model.fitTo(obs_data, params=["s", "b"])
    y = [model.nll(d=obs_data, s=sig) - nll_min for sig in x]
    plt.plot(x,y)

    # Draw profile nll
    model.fitTo(obs_data, params=["s", "b"])
    z = [model.profile(d=obs_data, poi="s", s=sig, nuisance=["b"], ) for sig in x]
    plt.plot(x,z)

    (ymin, ymax) = plt.ylim()
    (xmin, xmax) = plt.xlim()
    plt.ylim([0, ymax]) 
    plt.hlines(0.5, xmin, xmax)

    plt.xlabel('mu')
    plt.ylabel('profile likelihood(x)')
    plt.savefig("profile.pdf")
    plt.clf()


def test_data_plot(model):

    print "Test Data Plot"

    # Test the interval functionality (over data)
    plt.clf()
    model.fitTo(obs_data, params=["s", "b"])
    model.make_plot(interval=0.68)
    plt.xlabel('data')
    plt.savefig("plot.pdf")
    plt.clf()


def test_likelihood_plot(model, obs_data):

    print "Test Likelihood Plot"

    # Test the interval functionality (over data)
    plt.clf()
    model.d = obs_data
    model.fitTo(obs_data, params=["s", "b"])

    print "Likelihood Plot, Model State: ", str(model.total_state())
    x = model.s_var.linspace()
    z = [model.eval(s=point) for point in x]
    plt.plot(x, z)
    plt.xlabel('s')
    plt.ylabel('likelihood(s)')

    #x1,x2,y1,y2 = plt.axis()
    #model.make_plot(interval=0.68)
    #plt.xlabel('data')
    plt.savefig("likelihood.pdf")
    plt.clf()


def test_neyman(model, obs_data):

    print "Test Neyman"

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
    plt.clf()

    #print "inverted Neyman 4: ", model.invert_neyman(4, neyman)
    print "inverted Neyman data=6, bkg=", model.b, " s within:", model.invert_neyman(6, neyman)
    print "inverted Neyman data=8, bkg=", model.b, " s within:", model.invert_neyman(8, neyman)
    print "inverted Neyman data=10, bkg=", model.b, " s within:", model.invert_neyman(10, neyman)
    print "inverted Neyman data=12, bkg=", model.b, " s within:", model.invert_neyman(12, neyman)
    print "inverted Neyman data=14, bkg=", model.b, " s within:", model.invert_neyman(14, neyman)

    #print "inverted Neyman data=10: ", model.invert_neyman(10, neyman)
    #print "inverted Neyman data=12: ", model.invert_neyman(12, neyman)
    #print "inverted Neyman data=14: ", model.invert_neyman(14, neyman)


def test_data(model):

    print "Test data"

    model.set_data(10)
    print "Set data to 10.  Has data: ", model.get_data()
    print "Model Data Variable: ", model.data.name, model.data
    print "Model Data by direct access: ", model.d


def test_mc(model, obs_data):

    print "Test mc"

    model.fitTo(obs_data, params=["s", "b"])
    #samples = model.sample_mc(['d'], 2000)
    samples = model.sample_mc(['d'], 2000)
    values = [point['d'] for point in samples]

    # Plot the sampled values
    plt.clf()
    n, bins, patches = plt.hist(values, bins=20, range=[0,20], normed=1, facecolor='g')
    plt.xlabel('data'),
    plt.ylabel('Probability')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    plt.savefig("sampled_data_mc.pdf")
    plt.clf()


def test_mcmc(model, obs_data):

    print "Test mcmc"

    model.fitTo(obs_data, params=["s", "b"])
    #samples = model.sample_mc(['d'], 2000)
    samples = model.sample_mcmc(['d'], 2000)
    values = [point['d'] for point in samples]

    # Plot the sampled values
    plt.clf()
    n, bins, patches = plt.hist(values, bins=20, range=[0,20], normed=1, facecolor='g')
    plt.xlabel('data'),
    plt.ylabel('Probability')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    plt.savefig("sampled_data_mcmc.pdf")
    plt.clf()


if __name__ == "__main__":
    model = create_model()
    print_model(model)

    obs_data = 7


    #test_data(model)
    #test_data_plot(model)
    #test_likelihood_plot(model, obs_data)
    #test_minimization(model, obs_data)
    #test_neyman(model, obs_data)
    #test_profile(model, obs_data)
    test_mc(model, obs_data)
    test_mcmc(model, obs_data)
