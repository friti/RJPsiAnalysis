import ROOT
import numpy as np
import datetime
import os
production_tag = datetime.date.today().strftime('%Y%b%d')

print("Starting")
#data pu histo 
data_file = 'PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root'
input_data = ROOT.TFile.Open(data_file, 'read')
data_histo = input_data.Get("pileup")
data_histo_up = input_data.Get("pileup_plus")
data_histo_down = input_data.Get("pileup_minus")
print("Data histo loaded")

#mc nanoAOD input
#mc_file = 'RJPsi_mc_2020Nov25_89.root'
#mc_file = 'RJPsi_mc_10614.root'
mc_file = 'RJPsi_mc_'+production_tag+'.root'
input_mc = ROOT.TFile.Open(mc_file)
mc_tree = input_mc.Get("Events")
print("MC tree loaded")

#output
#newfile = ROOT.TFile("output_pu_reweight.root", "RECREATE")
newfile = ROOT.TFile('RJPsi_mc_pu_'+production_tag+'.root', "RECREATE")
newtree = mc_tree.CloneTree(0)
#newtree = mc_tree.CopyTree("","",-1,0)
pu_weights = np.array([0.])
pu_weights_up = np.array([0.])
pu_weights_down = np.array([0.])
puweight_branch = newtree.Branch("puWeight",pu_weights,'pu_weights/D')
puweight_up_branch = newtree.Branch("puWeight_up",pu_weights_up,'pu_weights_up/D')
puweight_down_branch = newtree.Branch("puWeight_down",pu_weights_down,'pu_weights_down/D')
print("New branches ready")

mc_histo = data_histo.Clone("autoPU")
mc_histo.Reset()
input_mc.Get("Events").Project("autoPU","Pileup_nTrueInt")

ret = mc_histo.Clone("hweights")
ret_up = mc_histo.Clone("hweights_up")
ret_down = mc_histo.Clone("hweights_down")


mc_vals = []
data_vals = []
data_vals_up = []
data_vals_down = []
#check if same number of bins
if(data_histo.GetNcells() != mc_histo.GetNcells()):
    print(data_histo.GetNcells, mc_histo.GetNcells)
    print("Numerator and Denominator have different number of bins")
    
#takes values of the bins of the two histograms
for i in range(data_histo.GetNcells()):
    mc_vals.append(mc_histo.GetBinContent(i))
    data_vals.append(data_histo.GetBinContent(i))
    data_vals_up.append(data_histo_up.GetBinContent(i))
    data_vals_down.append(data_histo_down.GetBinContent(i))
print("histo integral ",mc_histo.Integral())
    
#normalisation
scale_mc = 1.0/mc_histo.Integral()
mc_vals = [i * scale_mc for i in mc_vals]
scale_data = 1.0/data_histo.Integral()
data_vals = [i * scale_data for i in data_vals]
scale_data_up = 1.0/data_histo_up.Integral()
data_vals_up = [i * scale_data_up for i in data_vals_up]
scale_data_down = 1.0/data_histo_down.Integral()
data_vals_down = [i * scale_data_down for i in data_vals_down]
print("Values of the histograms normalized")

#weights = np.zeros(len(mc_vals),dtype = float)
#newfile = ROOT.TFile("pureweight.root", "RECREATE")
#newtree = mc_tree.CloneTree(0)
#puweight_branch = newtree.Branch("puWeight",weights,'weights/F')

weights = []
weights_up = []
weights_down = []

#does the ratio of the histos
for i in range(len(mc_vals)):
    if(mc_vals[i] !=0):
        weight = data_vals[i]/mc_vals[i]
        weight_up = data_vals_up[i]/mc_vals[i]
        weight_down = data_vals_down[i]/mc_vals[i]
    else:
        weight = 1
        weight_up = 1
        weight_down = 1
    #print(weight)
    weights.append(weight);
    weights_up.append(weight_up);
    weights_down.append(weight_down);
    #    puweight_branch.Fill()
print("Ratio of the histos done.")
#new mc_histo with weights
for i in range(len(weights)):
    ret.SetBinContent(i,weights[i])
    ret_up.SetBinContent(i,weights_up[i])
    ret_down.SetBinContent(i,weights_down[i])

print("weights histos made")
#for each event decides which weight to give

for i in range(mc_tree.GetEntries()):
    mc_tree.GetEntry(i)
    nvtx = newtree.Pileup_nTrueInt
    binx = max(1, min(ret.GetNbinsX(), ret.GetXaxis().FindBin(newtree.Pileup_nTrueInt)));
    binx_up = max(1, min(ret_up.GetNbinsX(), ret_up.GetXaxis().FindBin(newtree.Pileup_nTrueInt)));
    binx_down = max(1, min(ret_down.GetNbinsX(), ret_down.GetXaxis().FindBin(newtree.Pileup_nTrueInt)));
    if (nvtx < mc_histo.GetNbinsX()):
        weight = ret.GetBinContent(binx)
        weight_up = ret_up.GetBinContent(binx_up)
        weight_down = ret_down.GetBinContent(binx_down)
    else:
        weight = 1
        weight_up = 1
        weight_down = 1
    pu_weights[0] = weight
    pu_weights_up[0] = weight_up
    pu_weights_down[0] = weight_down
    newtree.Fill()
#output file
newfile.Write()
newfile.Close()
input_data.Close()
#delete old file
#print("Removing the old file.")
#os.system('rm RJPsi_mc_'+production_tag+'.root')
#newfile.cd()
#newtree.Write() 
#newfile.Close()

    



