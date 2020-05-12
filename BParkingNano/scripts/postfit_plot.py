import ROOT

f=ROOT.TFile("/work/friti/new/CMSSW_10_2_15/src/HiggsAnalysis/fit12May5/fitDiagnostics.root","r")

histo_m=f.Get("shapes_fit_s/rjpsi/mc_mu")
histo_m.SetFillColor(ROOT.kOrange+7)
histo_m.SetLineColor(ROOT.kOrange+7)

histo_t=f.Get("shapes_fit_s/rjpsi/mc_tau")
histo_t.SetFillColor(ROOT.kGreen)
histo_t.SetLineColor(ROOT.kGreen)

histo_x=f.Get("shapes_fit_s/rjpsi/mis_id")
histo_x.SetFillColor(ROOT.kBlue)
histo_x.SetLineColor(ROOT.kBlue)

histo_d=f.Get("shapes_fit_s/rjpsi/data")
histo_d.SetLineColor(ROOT.kBlack)

legend = ROOT.TLegend(0.7,0.66,0.94,0.88)

stack= ROOT.THStack("E_star","E_star")
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
stack.Draw("hist")
histo_d.Draw("E same")
legend.Draw("SAME")
f.Close()

