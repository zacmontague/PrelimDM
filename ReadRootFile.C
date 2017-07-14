
#include <map>

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

void overlayNormalizedHistograms(std::vector<TH1D*> histograms, TString fileName) {
    TCanvas* canvasOverlay = new TCanvas("canvasOverlay", "canvasOverlay", 1500, 1500);
    TLegend* legendOverlay = new TLegend(0.6, 0.5, 0.8, 0.8);

    int entries = histograms[0] -> Integral();
    
    for (unsigned int i(0); i < histograms.size(); ++i) {
	histograms[i] -> Scale(entries / histograms[i] -> Integral());
        histograms[i] -> Draw("SAME");
        if (i == 4) {
            histograms[i] -> SetLineColor(histograms.size() + 1);
        }
        else {
            histograms[i] -> SetLineColor(i + 1);
        }
        histograms[i] -> SetLineWidth(4);
        legendOverlay -> AddEntry(histograms[i], histograms[i] -> GetName(), "l");
    }
    legendOverlay -> Draw();
    canvasOverlay -> SaveAs(fileName);
}

void ReadRootFile() {
    std::map<TString, TH1D*> map;
    std::vector<TString> blahstring = {"taggerISR_Eta"};/*, "taggerISR_MinDeltaEta", "taggerISR_MinMassRatio", "taggerISR_MinPTRatio"};*/
    std::vector<std::vector<TH1D*>> histogramsToOverlay;
    
    TFile* file = new TFile("~/Delphes/OutFile.root", "READ");
    for(unsigned int i(0); i < blahstring.size(); ++i) {
	std::vector<TH1D*> v1 = {(TH1D*) file->Get(blahstring[i]+"_Signal"), (TH1D*)file->Get(blahstring[i] + "_Other")};
	overlayNormalizedHistograms(v1, "~/" + blahstring[i] + ".png");
    }
}
