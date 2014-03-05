import ROOT as rt
import rootlogon

rootlogon.style()

file = rt.TFile("outfile.root")

h2 = file.Get("h20x25")
h4_2 = file.Get("h40x25")
h4_5 = file.Get("h40x50")
h8 = file.Get("h80x50")

#data_file = rt.TFile("hlt_tree_run208390.root","READ")

h2.SetLineColor(rt.kRed)
h4_2.SetLineColor(rt.kBlue)
h4_5.SetLineColor(rt.kGreen)
h8.SetLineColor(rt.kCyan)

hists = [h2, h4_5, h4_2, h8]

for h in hists: h.SetLineWidth(2)

canvas = rt.TCanvas()
canvas.SetGridy(1)
hists[0].Draw()
hists[0].GetXaxis().SetTitle("r_{9} Cut")
hists[0].GetYaxis().SetTitle("Sig. Eff. HLT_26_18_mass70")
for h in hists[1:]:
    h.Draw("Same")

canvas.BuildLegend().SetFillColor(0)


raw_input("RAW INPUT")
