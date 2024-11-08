## \file
## \ingroup tutorial_tmva
## \notebook
## TMVA example, for classification
##  with following objectives:
##  * Apply a BDT with TMVA
##
## \macro_image
## \macro_output
## \macro_code
##
## \author Lailin XU

from ROOT import TMVA, TFile, TTree, TCut, TH1F, TCanvas, gROOT, TLegend
from subprocess import call
from os.path import isfile
from array import array
 
gROOT.SetStyle("ATLAS")

# Setup TMVA
TMVA.Tools.Instance()

# Reader. One reader for each application.
reader = TMVA.Reader("Color:!Silent")
reader_S = TMVA.Reader("Color:!Silent")
reader_B = TMVA.Reader("Color:!Silent")
 
# Inputs
# =============
# Load data
# An unknown sample
trfile = "real_data_float.root"
if not isfile(trfile):
  cmd = "python convert_tr.py"
  os.system(cmd)
data = TFile.Open(trfile)
tree = data.Get('data')

# Known signal
# trfile_S = "class_data.root"
trfile_S = "class_data_float.root"
data_S = TFile.Open(trfile_S)
tree_S = data_S.Get('signal;1')

# Known background
tree_B = data_S.Get('background;1')
 
# Set input variables. Do this for each reader
branches = {}
for branch in tree.GetListOfBranches():
    branchName = branch.GetName()
    branches[branchName] = array('f', [-999])
    tree.SetBranchAddress(branchName, branches[branchName])
    reader.AddVariable(branchName, branches[branchName])

branches_S = {}
for branch in tree_S.GetListOfBranches():
    branchName = branch.GetName()
    branches_S[branchName] = array('f', [-999])
    tree_S.SetBranchAddress(branchName, branches_S[branchName])
    reader_S.AddVariable(branchName, branches_S[branchName])

branches_B = {}
for branch in tree_B.GetListOfBranches():
    branchName = branch.GetName()
    branches_B[branchName] = array('f', [-999])
    tree_B.SetBranchAddress(branchName, branches_B[branchName])
    reader_B.AddVariable(branchName, branches_B[branchName])
 
# Book method(s)
# =============
# BDT
methodName1 = "BDT"
weightfile = 'dataset/weights/TMVAClassification_{0}.weights.xml'.format(methodName1)
reader.BookMVA( methodName1, weightfile )
reader_S.BookMVA( methodName1, weightfile )
reader_B.BookMVA( methodName1, weightfile )
# BDTG
methodName2 = "BDTG"
weightfile = 'dataset/weights/TMVAClassification_{0}.weights.xml'.format(methodName2)
reader.BookMVA( methodName2, weightfile )
reader_S.BookMVA( methodName2, weightfile )
reader_B.BookMVA( methodName2, weightfile )
# MLP
methodName3 = "MLP"
weightfile = 'dataset/weights/TMVAClassification_{0}.weights.xml'.format(methodName3)
reader.BookMVA( methodName3, weightfile )
reader_S.BookMVA( methodName3, weightfile )
reader_B.BookMVA( methodName3, weightfile )

# Loop events for evaluation
# ================

# Book histograms
nbins, xmin, xmax=20, -1, 1
# Signal
tag = "S"
hname="BDT_{0}".format(tag)
h1 = TH1F(hname, hname, nbins, xmin, xmax)
h1.Sumw2()
hname="BDTG_{0}".format(tag)
h2 = TH1F(hname, hname, nbins, xmin, xmax)
h2.Sumw2()
nbins, xmin, xmax=20, 0, 1
hname="MLP_{0}".format(tag)
h7 = TH1F(hname, hname, nbins, xmin, xmax)
h7.Sumw2()

nevents = tree_S.GetEntries()
for i in range(nevents):
  tree_S.GetEntry(i)

  BDT = reader_S.EvaluateMVA(methodName1)
  BDTG = reader_S.EvaluateMVA(methodName2)
  MLP = reader_S.EvaluateMVA(methodName3)
  h1.Fill(BDT)
  h2.Fill(BDTG)
  h7.Fill(MLP)

# Background
nbins, xmin, xmax=20, -1, 1
tag = "B"
hname="BDT_{0}".format(tag)
h3 = TH1F(hname, hname, nbins, xmin, xmax)
h3.Sumw2()
hname="BDTG_{0}".format(tag)
h4 = TH1F(hname, hname, nbins, xmin, xmax)
h4.Sumw2()
nbins, xmin, xmax=20, 0, 1
hname="MLP_{0}".format(tag)
h8 = TH1F(hname, hname, nbins, xmin, xmax)
h8.Sumw2()

nevents = tree_B.GetEntries()
for i in range(nevents):
  tree_B.GetEntry(i)

  BDT = reader_B.EvaluateMVA(methodName1)
  BDTG = reader_B.EvaluateMVA(methodName2)
  MLP = reader_B.EvaluateMVA(methodName3)
  h3.Fill(BDT)
  h4.Fill(BDTG)
  h8.Fill(MLP)

# New sample
nbins, xmin, xmax=20, -1, 1
tag = "N"
hname="BDT_{0}".format(tag)
h5 = TH1F(hname, hname, nbins, xmin, xmax)
h5.Sumw2()
hname="BDTG_{0}".format(tag)
h6 = TH1F(hname, hname, nbins, xmin, xmax)
h6.Sumw2()
nbins, xmin, xmax=20, 0, 1
hname="MLP_{0}".format(tag)
h9 = TH1F(hname, hname, nbins, xmin, xmax)
h9.Sumw2()

nevents = tree.GetEntries()
ntot = 0
l_ns = [0]*3
for i in range(nevents):
  tree.GetEntry(i)

  BDT = reader.EvaluateMVA(methodName1)
  BDTG = reader.EvaluateMVA(methodName2)
  MLP = reader.EvaluateMVA(methodName3)
  h5.Fill(BDT)
  h6.Fill(BDTG)
  h9.Fill(MLP)

  # apply a cut to classify the event as a signal
  ntot +=1
  if BDT >0: l_ns[0] += 1
  if BDTG >0: l_ns[1] += 1
  if MLP>0.5: l_ns[2] += 1

# Helper function to normalize hists
def norm_hists(h):

  h_new = h.Clone()
  hname = h.GetName() + "_normalized"
  h_new.SetName(hname)
  h_new.SetTitle(hname)
  ntot = h.Integral()
  if ntot!=0:
    h_new.Scale(1./ntot)

  return h_new

# Plotting
myc = TCanvas("c", "c", 800, 600)
myc.SetFillColor(0)
myc.cd()

# Compare the performance for BDT
nh1 = norm_hists(h1)
nh1.GetXaxis().SetTitle("BDT")
nh1.GetYaxis().SetTitle("A.U.")
nh1.Draw("hist")
nh3 = norm_hists(h3)
nh3.SetLineColor(2)
nh3.SetMarkerColor(2)
nh3.Draw("same hist")
nh5 = norm_hists(h5)
nh5.SetLineColor(4)
nh5.SetMarkerColor(4)
nh5.Draw("same")

ymin = 0
ymax = max(nh1.GetMaximum(), nh3.GetMaximum(), nh5.GetMaximum())
nh1.GetYaxis().SetRangeUser(ymin, ymax*1.5)

# Draw legends
lIy = 0.92
lg = TLegend(0.60, lIy-0.25, 0.85, lIy)
lg.SetBorderSize(0)
lg.SetFillStyle(0)
lg.SetTextFont(42)
lg.SetTextSize(0.04)
lg.AddEntry(nh1, "Signal", "l")
lg.AddEntry(nh3, "Background", "l")
lg.AddEntry(nh5, "Data", "l")
lg.Draw()

myc.Draw()
myc.SaveAs("BDT_cla_app_1.png")

# Compare the performance for BDTG
nh1 = norm_hists(h2)
nh1.GetXaxis().SetTitle("BDTG")
nh1.GetYaxis().SetTitle("A.U.")
nh1.Draw("hist")
nh3 = norm_hists(h4)
nh3.SetLineColor(2)
nh3.SetMarkerColor(2)
nh3.Draw("same hist")
nh5 = norm_hists(h6)
nh5.SetLineColor(4)
nh5.SetMarkerColor(4)
nh5.Draw("same")

ymin = 0
ymax = max(nh1.GetMaximum(), nh3.GetMaximum(), nh5.GetMaximum())
nh1.GetYaxis().SetRangeUser(ymin, ymax*1.5)

# Draw legends
lIy = 0.92
lg = TLegend(0.60, lIy-0.25, 0.85, lIy)
lg.SetBorderSize(0)
lg.SetFillStyle(0)
lg.SetTextFont(42)
lg.SetTextSize(0.04)
lg.AddEntry(nh1, "Signal", "l")
lg.AddEntry(nh3, "Background", "l")
lg.AddEntry(nh5, "Data", "l")
lg.Draw()

myc.Draw()
myc.SaveAs("BDT_cla_app_2.png")

# Compare the performance for MLP
nh1 = norm_hists(h7)
nh1.GetXaxis().SetTitle("MLP")
nh1.GetYaxis().SetTitle("A.U.")
nh1.Draw("hist")
nh3 = norm_hists(h8)
nh3.SetLineColor(2)
nh3.SetMarkerColor(2)
nh3.Draw("same hist")
nh5 = norm_hists(h9)
nh5.SetLineColor(4)
nh5.SetMarkerColor(4)
nh5.Draw("same")

ymin = 0
ymax = max(nh1.GetMaximum(), nh3.GetMaximum(), nh5.GetMaximum())
nh1.GetYaxis().SetRangeUser(ymin, ymax*1.5)

# Draw legends
lIy = 0.92
lg = TLegend(0.60, lIy-0.25, 0.85, lIy)
lg.SetBorderSize(0)
lg.SetFillStyle(0)
lg.SetTextFont(42)
lg.SetTextSize(0.04)
lg.AddEntry(nh1, "Signal", "l")
lg.AddEntry(nh3, "Background", "l")
lg.AddEntry(nh5, "Data", "l")
lg.Draw()

myc.Draw()
myc.SaveAs("BDT_cla_app_3.png")


# Classification results
# =======================
print("\nTest classification: real data ")
rate_s = [0]*3
for i, ns in enumerate(l_ns):
  rate_s[i] = ns/float(ntot)

print("Estimated signal events: BDT {0} / {1} = {2}".format(l_ns[0], ntot, rate_s[0]))
print("Estimated signal events: BDTG {0} / {1} = {2}".format(l_ns[1], ntot, rate_s[1]))
print("Estimated signal events: MLP {0} / {1} = {2}".format(l_ns[2], ntot, rate_s[2]))
