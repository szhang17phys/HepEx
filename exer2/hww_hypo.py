## \file
## \ingroup hypotest
## \notebook
##
## \macro_image
## \macro_output
## \macro_code
##
## \author Lailin XU


import ROOT as R
from math import log, sqrt, fabs

# HWW simplified model: https://arxiv.org/abs/1412.2641

N1 = 407.
s1 = 48.6
b1 = 335.
N2, s2, b2 = 143., 16., 129.

# Treat ee+mumu and emu as two separate measurements
# =================

hs1 = N1 - b1
hs2 = N2 - b2

# Delta log-likelihood
lam = 2*( - hs1 + N1 * log( (hs1+b1)/b1 ))
lam += 2*( - hs2 + N2 * log( (hs2+b2)/b2 ))

# P-value
p = R.Math.chisquared_cdf_c(lam, 2.)
# Convert p-value to signifiance
nsig = R.Math.normal_quantile_c(p, 1.)
# Use RooStats to convert p-value to signifiance
nsig_rs = R.RooStats.PValueToSignificance(p)
# Two-sided siginifance
nsig2 = R.Math.normal_quantile_c(p*0.5, 1.)

print("Lambda = {0}, p-value= {1}, significance= {2} {3}".format(lam, p, nsig, nsig2))

# One single POI
# =================
# Maximize the likelihood to get the best estimate of mu
a = (s1 + s2)*s1*s2
b = (s1+s2)*(b1*s2 + b2*s1) - (N1+N2)*s1*s2
c = b1*b2*(s1+s2) - N1*s1*b2 - N2*s2*b1

# MLE of mu
mu = (-b + sqrt( pow(b, 2) - 4.*a*c))/(2.*a)

print("mu = {0}".format(mu))

hs1 = mu*s1
hs2 = mu*s2

lam = 2*( - hs1 + N1 * log( (hs1+b1)/b1 ))
lam += 2*( - hs2 + N2 * log( (hs2+b2)/b2 ))

p = R.Math.chisquared_cdf_c(lam, 1.)
nsig = R.Math.normal_quantile_c(p, 1.)
nsig_rs = R.RooStats.PValueToSignificance(p)
nsig2 = R.Math.normal_quantile_c(p*0.5, 1.)

print("Lambda = {0}, p-value= {1}, significance= {2} {3}".format(lam, p, nsig, nsig2))

# Calculate NLL explicitly
x = 0
lnL0 = ( -x* s1 - b1 + N1*log(x*s1+b1) ) + (-x*s2 - b2 + N2*log(x*s2+b2) )
print("mu={0}, NLL= {1}".format(x, lnL0))
x = mu
lnL1 = ( -x* s1 - b1 + N1*log(x*s1+b1) ) + (-x*s2 - b2 + N2*log(x*s2+b2) )
print("mu={0}, NLL= {1}".format(x, lnL1))
dNLL = 2*(lnL1-lnL0)
print("dNLL={0}, Z0= {1}".format(dNLL, sqrt(dNLL)))


# Use RooFit/RooStat to do the hypotest
# =================
# Instantiate a workspace
w = R.RooWorkspace("w")

# Create pdf components
w.factory("expr::s_emu('mu*n_emu',mu[1,0,10],n_emu[48.6])")
w.factory("sum:nexp_emu(s_emu,b_emu[335])")
w.factory("expr::s_ll('mu*n_ll',mu,n_ll[16])")
w.factory("sum:nexp_ll(s_ll,b_ll[129])")
w.factory("Poisson:emu(nobs_emu[0,1000],nexp_emu)")
w.factory("Poisson:ll(nobs_ll[0,1000],nexp_ll)")

# Create the total model
w.factory("PROD:model(emu,ll)")

# Data
r_emu = w.var("nobs_emu")
r_ll = w.var("nobs_ll")
data = R.RooDataSet("ds", "ds", R.RooArgSet(r_emu, r_ll))

r_emu.setVal(407)
r_ll.setVal(143)
data.add(R.RooArgSet(r_emu, r_ll))

w.Import(data, R.RooFit.Rename("ObsData"))

model = w.pdf("model")
model.Print()

# Create the ModelConfig
mc=R.RooStats.ModelConfig("ModelConfig", w)
# Set up the Model
mc.SetPdf(w.pdf("model"))
mc.SetParametersOfInterest(R.RooArgSet(w.var("mu")))
mc.SetObservables(R.RooArgSet(w.var("nobs_emu"), w.var("nobs_ll")))

mc.Print()

w.Import(mc)

# Take a peek at the workspac
w.Print("t")

# Fitting
model.fitTo(data)

data.Print("v")

# Asymptotic calculator
# =================

# The S+B model (Alternative hypo)
sbModel = w.obj("ModelConfig")
pois = sbModel.GetParametersOfInterest()
pois[0].setVal(1)
poi = pois.first()
sbModel.SetSnapshot(R.RooArgSet(pois))

# PDF
pdf = sbModel.GetPdf()

# save snapshot before any fit has been done
params = pdf.getParameters(data)
snapshotName_init = "snapshot_paramsVals_initial"
w.saveSnapshot(snapshotName_init, params)

# The B model (Null hypo)
bModel = sbModel.Clone()
bModel.SetName("B_only_model")
pois[0].setVal(0)
bModel.SetSnapshot(R.RooArgSet(pois))
bModel.Print()

w.Print()

# Asymptotic calculator
ac = R.RooStats.AsymptoticCalculator(data, sbModel, bModel)
ac.SetOneSidedDiscovery(True)

# Get the hypo test result
asResult = ac.GetHypoTest()
asResult.Print()
pvalue_as = asResult.NullPValue()


# By hand calculation
# =======================
w.loadSnapshot(snapshotName_init)
sbModel = w.obj("ModelConfig")
pdf = sbModel.GetPdf()
# Get the nuisance parameters and global observables
constrainedParams = sbModel.GetNuisanceParameters()
glbObs = sbModel.GetGlobalObservables()
# Create the neg-log-likelihood
nll_sb = pdf.createNLL(data, R.RooFit.Constrain(constrainedParams), R.RooFit.GlobalObservables(glbObs),
                                       R.RooFit.NumCPU(2), R.RooFit.Optimize(2))
nllval = nll_sb.getVal()
print("Starting NLL value:", nllval)

# Do the minimization
minim = R.RooMinimizer(nll_sb)
strategy = R.Math.MinimizerOptions.DefaultStrategy()
minim.setStrategy(strategy)
minim.optimizeConst(2)
minimizer = R.Math.MinimizerOptions.DefaultMinimizerType()
algorithm = R.Math.MinimizerOptions.DefaultMinimizerAlgo()
print("\n =========== Unconditinal fit =========\n")
status = minim.minimize(minimizer, algorithm)

obs_nll_min = nll_sb.getVal()
reverse = (poi.getVal() < 0)

# Fix POI to 0 (B-only model) and do the minimization again
print("\n =========== Conditinal fit =========\n")
w.loadSnapshot(snapshotName_init)
poi.setVal(0)
poi.setConstant(1)

status = minim.minimize(minimizer, algorithm)

obs_nll_min_bkg = nll_sb.getVal()

# The asymptotic statistic: q0 = nll(b-only, mu=0) - nll(s+b)
# Significance: Z = sqrt( 2*q0 )
obs_q0 = 2*(obs_nll_min_bkg - obs_nll_min)
# Check the sign: excess or deficit? 
if reverse: obs_q0 = -obs_q0
sign = 0
if obs_q0!=0: sign = obs_q0 / fabs(obs_q0)
obs_sig = sign*sqrt(fabs(obs_q0));
print("\nUnconditional NLL value:", obs_nll_min)
print("Conditional NLL value:", obs_nll_min_bkg)
print("dNLL={0}, Z0= {1}".format(obs_q0, sqrt(obs_q0)))
print("==> Asymmptotic signficance: ", obs_sig)

