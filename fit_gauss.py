
class PlotPoints:
    def __init__(self, fig):
        self.fig = fig
        self.xcoords = np.array([])
        self.ycoords = np.array([])
        self.flag = False

    def connect(self):
        self.cidpress = self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        self.cidclick = self.fig.canvas.mpl_connect('button_release_event', self.on_click)
        print("  Hit 'q' to skip these data.\n  Click another button, or the mouse, on initial guesses for ingress and egress:")

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cidpress)
        self.fig.canvas.mpl_disconnect(self.cidclick)

    def on_press(self, event):
        if 'q' in event.key:
            self.flag = True
            self.disconnect()
            closeplot()
            print("  ")
            return
        print('  added point at %.1f, %.1f' % (event.xdata, event.ydata))
        self.xcoords = np.append(self.xcoords, event.xdata)
        self.ycoords = np.append(self.ycoords, event.ydata)
        if self.xcoords.size == 2:
            self.disconnect()
            closeplot()
            print("  ")

    def on_click(self, event):
        print('  added point at %.1f, %.1f' % (event.xdata, event.ydata))
        self.xcoords = np.append(self.xcoords, event.xdata)
        self.ycoords = np.append(self.ycoords, event.ydata)
        if self.xcoords.size == 2:
            self.disconnect()
            closeplot()
            print("  ")

    def gaussPars(self):
        sep = np.fabs(self.xcoords[0] - self.xcoords[1])
        peak = np.fabs(self.ycoords).mean()
        t0 = self.xcoords.mean()
        sigma = sep/20
        return dict(t0=t0, peak=peak, sep=sep, log_sigma2=np.log(sigma**2))


class TwoGaussians(Model):
    parameter_names = ('t0', 'sep', 'peak', 'log_sigma2')

    def get_value(self, x):
        mu1, mu2 = self.t0 - self.sep/2, self.t0 + self.sep/2
        g1 = np.exp(-0.5 * (x-mu1)**2 * np.exp(-self.log_sigma2))
        g2 = np.exp(-0.5 * (x-mu2)**2 * np.exp(-self.log_sigma2))
        return -self.peak * g1 + self.peak * g2



# Define a cost function for MCMC
def log_like(params, y, gp):
    # print(params)
    gp.set_parameter_vector(params)
    return gp.log_likelihood(y)

# Define a cost function for scipy.minimize
def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)
def grad_neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.grad_log_likelihood(y)[1]

def fit_gauss(x, y, yerr):
    fig, ax = plt.subplots()
    plt.plot(x, y)
    gauss = PlotPoints(fig)
    gauss.connect()
    plt.show()

    if gauss.flag:
        print("  No eclipse taken from these data.")
        continue

    kwargs = gauss.gaussPars()
    # hold values close to initial guesses
    bounds = dict(
        t0=(kwargs['t0']-kwargs['sep']/8, kwargs['t0']+kwargs['sep']/8),
        sep=(0.9*kwargs['sep'], 1.1*kwargs['sep']),
        log_sigma2=(np.log(kwargs['sep']**2/10000), np.log(kwargs['sep']**2/25)),
        peak=(0.9*kwargs['peak'], 1.1*kwargs['peak'])
    )
    kwargs['bounds'] = bounds

    mean_model = TwoGaussians(**kwargs)

    mean, median, std = sigma_clipped_stats(y)
    delta_t = np.mean(np.diff(x))*5
    kernel = terms.RealTerm(log_a=np.log(std**2), log_c=-np.log(delta_t))
    gp = celerite.GP(kernel, mean=mean_model, fit_mean=True)
    gp.compute(x, yerr)
    # print("  Initial log-likelihood: {0}".format(gp.log_likelihood(y)))


    # Fit for the maximum likelihood parameters
    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()


    # Find a solution using Stu's method
    soln = minimize(neg_log_like, initial_params, jac=grad_neg_log_like,
                    method="L-BFGS-B", bounds=bounds, args=(y, gp))
    if not soln.success:
        print('  Warning: may not have converged')
        print(soln.message)

    gp.set_parameter_vector(soln.x)
    mean_model.set_parameter_vector(gp.get_parameter_vector()[2:])

    out = soln['x']
    t_ecl = out[2]

    print("  Using MCMC to characterise error at peak likelihood...")


    # Use an MCMC model, starting from the solution we found, to model the errors
    ndim     = 6
    nwalkers = 50

    # Initial positions. Scatter by 0.00001, as this is one above the order of magnitude of the error
    #  we expect on t_ecl.
    p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
    scatter = 0.0001/t_ecl
    p0 *= scatter
    p0 += 1. - scatter
    p0 = np.transpose(np.repeat(out, nwalkers).reshape((ndim, nwalkers))) *p0

    # Construct a sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_like, args=[y, gp], threads=1)


    width=30

    # Burn in
    print("")
    nsteps = 200
    start_time = time.time()
    for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
        n = int((width+1) * float(i) / nsteps)
        sys.stdout.write("\r  Burning in...    [{}{}]".format('#'*n, ' '*(width - n)))
    pos, prob, state = result

    # Data
    sampler.reset()
    nsteps = 300

    start_time = time.time()
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
        n = int((width+1) * float(i) / nsteps)
        sys.stdout.write("\r  Sampling data... [{}{}]".format('#'*n, ' '*(width - n))) 
    print("")

    # corner.corner(sampler.flatchain, labels=['???', '???', 't_ecl', '???', 'a', 'b'])
    # plt.show()

    t_ecl = np.mean(sampler.flatchain[:,2])
    err = np.std(sampler.flatchain[:,2])
    sep = np.mean(sampler.flatchain[:,3])

    print("    Got a solution: {:.7f}+/-{:.7f}\n".format(t_ecl, err))#

    return t_ecl, err, sep


## example
centre_of_ecl, error, seperation_of_eclipses = fit_gauss(x, y, yerr)