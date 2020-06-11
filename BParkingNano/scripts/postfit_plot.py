import ROOT

ROOT.gROOT.SetBatch()
f=ROOT.TFile("/work/friti/new/CMSSW_10_2_15/src/HiggsAnalysis/fit22Maytau2/fitDiagnostics.root","r")

histo_m=f.Get("shapes_fit_s/rjpsi/mc_mu")
histo_m.SetFillColor(ROOT.kRed)
histo_m.SetLineColor(ROOT.kRed)

histo_t=f.Get("shapes_fit_s/rjpsi/mc_tau")
histo_t.SetFillColor(ROOT.kGreen)
histo_t.SetLineColor(ROOT.kGreen)

histo_x=f.Get("shapes_fit_s/rjpsi/mis_id")
histo_x.SetFillColor(ROOT.kBlue)
histo_x.SetLineColor(ROOT.kBlue)

histo_d=f.Get("shapes_fit_s/rjpsi/data")
histo_d.SetLineColor(ROOT.kBlack)

legend = ROOT.TLegend(0.55,0.66,0.85,0.88)

stack= ROOT.THStack("","")
stack.Add(histo_m)
legend.AddEntry(histo_m, "MC mu", "f")
stack.Add(histo_t)
legend.AddEntry(histo_t, "MC tau", "f")
stack.Add(histo_x)
legend.AddEntry(histo_x, "MisID", "f")
legend.AddEntry(histo_d, "Data", "lp")

maximum = max(histo_d.GetMaximum(), stack.GetMaximum())
stack.SetMaximum(maximum*1.3)


c=ROOT.TCanvas()
#stack.SetTitle("m_{miss}^{2} post fit;m_{miss}^{2} [GeV^{2}];Counts")
#stack.SetTitle("Q^{2} post fit;Q^{2} [GeV^{2}];Counts")
stack.SetTitle("tau")
stack.Draw("hist")
histo_d.Draw("E same")
#ROOT.gPad.SetLogy()
'''
for i in range(1,7):
    stack.GetXaxis().SetBinLabel(i,str(i+2))
    histo_d.GetXaxis().SetBinLabel(i,str(i+2))
'''

legend.Draw("SAME")
c.Update()
c.SaveAs("fit.png")
raw_input("ciaociao")
f.Close()

