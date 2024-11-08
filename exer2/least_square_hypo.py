## \file
## \ingroup tutorial_pyroot
## \notebook
##
## \macro_image
## \macro_output
## \macro_code
##
## \author Lailin XU

# Import the ROOT libraries
import ROOT as R
from math import pow, sqrt, fabs
R.gROOT.SetStyle("ATLAS")
import array

# Input data
y1 = 0.9
e1 = 0.1
y2 = 1.4 
e2 = 0.2

# Hypo test: Wald test
# ================

# Define a test stat: y = y1 - y2, y ~ Normal(0, sigma)
y = y2 - y1
e = sqrt(pow(e1, 2) + pow(e2, 2))

# H0: y=0, Wald test
w = y/e
p = 2*R.Math.normal_cdf(-fabs(w))

print("test_stat: {0}, p-value: {1}".format(w, p))

# Hypo test: least square
# ================

# MLE y
y = (y1/pow(e1, 2) + y2/pow(e2, 2)) / (1/pow(e1, 2) + 1/pow(e2, 2))
ye = 1/sqrt( 1/pow(e1, 2) + 1/pow(e2, 2) )

chi2 = pow( (y1-y)/e1, 2) + pow( (y2-y)/e2, 2)

# P-value
p = R.Math.chisquared_cdf_c(chi2, 1.)
print("MLE y={2}, ye={3}, Chi2= {0}, p-value= {1}".format(chi2, p, y, ye))

# Plotting
myc = R.TCanvas("c", "c", 800, 600)
myc.SetFillColor(0)

myc.cd()

n = 2
ax = array.array("f", [1, 2])
axe = array.array("f", [0]*n)
ay = array.array("f", [y1, y2])
aye = array.array("f", [e1, e2])
gr = R.TGraphErrors(n, ax, ay, axe, aye)

gr_comb = R.TGraphErrors(n)
gr_comb.SetPoint(0, 1, y)
gr_comb.SetPointError(0, 0, ye)
gr_comb.SetPoint(1, 2, y)
gr_comb.SetPointError(1, 0, ye)
gr_comb.SetFillColor(R.kGreen)
gr_comb.SetFillStyle(3002)

gr.Draw("AP")
gr.GetXaxis().SetTitle("Measurements")
gr.GetYaxis().SetTitle("y")

gr_comb.Draw("same C 3")
myc.SaveAs("exer2_hypo.png")
