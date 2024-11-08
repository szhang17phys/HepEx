## \file
## \ingroup tutorial_tmva
## \notebook
## A simply script to convert tree branches
##
## \macro_image
## \macro_output
## \macro_code
##
## \author Lailin XU

from ROOT import TMVA, TFile, TTree, TCut, TH1F, TCanvas, gROOT, TLegend
from array import array


def convert(trfile, tr, tr_new):
  """
  convert tree branches Double to Float
  """

  data = TFile.Open(trfile)
  tree = data.Get(tr)
  
  trfile_new = trfile.replace(".root", "_float.root")
  data_new = TFile.Open(trfile_new, "UPDATE")
  tree_new = TTree(tr_new, tr_new)
  
  x1 = array('f', [0])
  x2 = array('f', [0])
  tree_new.Branch("X1", x1, "X1/F")
  tree_new.Branch("X2", x2, "X2/F")


  nevents = tree.GetEntries()
  for i in range(nevents):
    tree.GetEntry(i)
    x1[0] = tree.X1
    x2[0] = tree.X2
  
    tree_new.Fill()
  
  data_new.cd()
  tree_new.Write()
  data_new.Close()
  
  data.Close()
  
  
  
if __name__ == "__main__":

  trfile = "real_data.root"
  convert(trfile, "data;2", "data")

  trfile = "class_data.root"
  convert(trfile, "signal;1", "signal")

  trfile = "class_data.root"
  convert(trfile, "background;1", "background")
