import autograd.numpy as np
from approxcv.models.abstract_lgiidregr import AbstractiidLGR
from approxcv.acvij import IJ
from approxcv.acvns import NS
from approxcv.models.examples import LOL
import regression.helpers as hp
import os
import time
from pandas import read_csv

args = hp.get_args()
method = args.method
integrate = False if method == 'MCV' else args.integrate
is_ij = method == "IJ"

modelfolder = "../pretrained/lol"
savefolder = "./results/lol"

dat = hp.load_data(modelfolder)
phi = read_csv("%s/phi.csv" % modelfolder).to_numpy()[:,1]
y = hp.approx_y_poisson(np.log(1/phi), dat['kills'].to_numpy())
cvfolds = dat['player'].factorize()[0]
J = np.max(cvfolds)


prefix = ''
# draws = hp.load_samps(modelfolder, prefix).to_numpy()
(parnames, params) = hp.load_params(modelfolder, prefix)
X = hp.load_x(modelfolder, prefix)
N,P = X.shape
prior = hp.prior_info(modelfolder, prefix)

pm = hp.ParamManager(parnames, params, X)
if integrate:
    reff_mask = ['player' in x for x in pm.parnames][0:P]
    pm = pm.remove_coeffs(reff_mask).save_params()

model = LOL(pm.X, y, prior, pm.parnames,  integrate = integrate)
ij_start = time.time()
acv = IJ(model, pm.params, N) if is_ij else NS(model, pm.params, N)
ij_time = time.time() - ij_start


prefix_int = "-c" if integrate else ""
with open(savefolder + '/%s_yhats.csv' % (method + prefix_int), 'w') as yhatfile:
    yhatfile.write("i,yhat,loop,loop_time\n")
    for j in range(J + 1):
        print(j)
        o = [i for i in range(N) if cvfolds[i] == j]
        # check for any missing parameters
        if not is_ij:
            acv, model = pm.remove_missing(o).transfer_to_model(acv, model)
            pm = pm.reset_params()
        yhats = hp.compute_yhats(acv, model, o, optim = method == "MCV")

        timedif = yhats['time']
        timedif += ij_time if j == 0 and is_ij else 0
        yhatfile.write(hp.yhat_line('', yhats['yhat'], o, j + 1, timedif))

