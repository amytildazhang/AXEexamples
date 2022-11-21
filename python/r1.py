from autograd import jacobian
import autograd.numpy as np
from approxcv.models.abstract_lgiidregr import AbstractiidLGR
from approxcv.acvij import IJ
from approxcv.acvns import NS
from approxcv.models.examples import Radon
import regression.helpers as hp
import os
import time

# possible values for method: NS, IJ, MCV
args = hp.get_args()
method = args.method
integrate = False if method == 'MCV' else args.integrate
is_ij = method == "IJ"

modelfolder = "../data-raw/pretrained/r1"
savefolder = "./results/r1"


dat = hp.load_data(modelfolder)
y = dat['log_radon'].to_numpy()
cvfolds = dat['county'].factorize()[0]
J = np.max(cvfolds)

prefix_int = "-c" if integrate else ""
with open(savefolder + '/%s_yhats.csv' % (method + prefix_int), 'w') as yhatfile:
    yhatfile.write("model,i,yhat,loop,loop_time\n")
    for mod in np.arange(1,4):
        print(mod)
        prefix = "model%s_" % mod
        (parnames, params) = hp.load_params(modelfolder, prefix)
        X = hp.load_x(modelfolder, prefix)
        X[:,0] = 1.0
        N,P = X.shape
        prior = hp.prior_info(modelfolder, prefix)
        # pm_full = hp.optim_fulldata(parnames, params, X, Radon, y, prior)

        # pm = hp.ParamManager(pm_full.parnames, pm_full.params, pm_full.X)
        pm = hp.ParamManager(parnames, params, X)
        if integrate:
            county_mask = ['county' in x for x in pm.parnames][0:P]
            pm = pm.remove_coeffs(county_mask).save_params()

        model = Radon(pm.X, y, prior, pm.parnames,  integrate = integrate)
        ij_start = time.time()
        acv = IJ(model, pm.params, N) if is_ij else NS(model, pm.params, N)
        ij_time = time.time() - ij_start
        for j in range(J + 1):
            o = [i for i in range(N) if cvfolds[i] == j]
            # check for any missing parameters
            if not is_ij:
                acv, model = pm.remove_missing(o).transfer_to_model(acv, model)
                pm = pm.reset_params()

            yhats = hp.compute_yhats(acv, model, o, optim = method == "MCV")
            timedif = yhats['time']
            timedif += ij_time if j == 0 and is_ij else 0

            yhatfile.write(hp.yhat_line('%s,' % mod, yhats['yhat'], o, j + 1, timedif))
