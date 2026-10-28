#include "Plotting.h"

#include <TCanvas.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TPaletteAxis.h>
#include <iostream>

#include <algorithm>
#include <cmath>
#include <memory>

namespace plotting {

// Text in the top right corner of every plot
const char* lumiLabel = "109 fb^{-1} (13.6 TeV)";

void setStyle()
{
    gStyle->SetOptStat(0);   // no statistics box
    gStyle->SetOptTitle(0);  // no histogram title above the plot
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetPadTickX(1);  // ticks on all four sides
    gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(0);
    TH1::AddDirectory(false);  // the histogram copies made below are not attached to any open file
}

// Move the underflow and overflow into the first and last bin.
void addOverflowToEdgeBins(TH1D& h)
{
    const int n = h.GetNbinsX();
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1));
    h.SetBinError(1, std::hypot(h.GetBinError(0), h.GetBinError(1)));
    h.SetBinContent(n, h.GetBinContent(n) + h.GetBinContent(n + 1));
    h.SetBinError(n, std::hypot(h.GetBinError(n), h.GetBinError(n + 1)));
    h.SetBinContent(0, 0.);
    h.SetBinContent(n + 1, 0.);
}

void drawDataMC(const TH1D* data, const std::vector<TH1D*>& mcStack, const std::string& fileName, bool logY)
{
    // ---- Copies of the histograms, so the originals stay unchanged ---------
    TH1D hData(*data);
    addOverflowToEdgeBins(hData);
    hData.SetMarkerStyle(20);
    hData.SetMarkerSize(0.8);
    hData.SetLineColor(kBlack);

    THStack stack;
    std::vector<std::unique_ptr<TH1D>> hMC;  // unique_ptr deletes the copies at the end of the function
    for (const TH1D* h : mcStack) {
        hMC.push_back(std::make_unique<TH1D>(*h));
        addOverflowToEdgeBins(*hMC.back());
        hMC.back()->SetLineColor(h->GetFillColor());
        stack.Add(hMC.back().get());
    }

    // Sum of all MC histograms, drawn as hatched statistical uncertainty band
    TH1D hMCTotal(*hMC.front());
    for (size_t i = 1; i < hMC.size(); ++i) hMCTotal.Add(hMC[i].get());
    hMCTotal.SetFillStyle(3354);
    hMCTotal.SetFillColor(kGray + 2);
    hMCTotal.SetLineWidth(0);
    hMCTotal.SetMarkerSize(0);

    // ---- Canvas with two pads: distributions on top, ratio below -----------
    TCanvas canvas("canvas", "", 800, 800);
    TPad upper("upper", "", 0., 0.30, 1., 1.);
    TPad lower("lower", "", 0., 0., 1., 0.30);
    upper.SetMargin(0.15, 0.04, 0.02, 0.10);  // left, right, bottom, top
    lower.SetMargin(0.15, 0.04, 0.35, 0.03);
    upper.Draw();
    lower.Draw();

    // ---- Upper pad: data and stacked MC -----------------------------------
    upper.cd();
    upper.SetLogy(logY);

    const double yMax = std::max(hData.GetMaximum(), hMCTotal.GetMaximum());
    stack.SetMaximum(logY ? yMax * 500. : yMax * 1.5);  // leave room for the legend
    stack.SetMinimum(logY ? 0.5 : 0.);
    stack.Draw("hist");
    stack.GetXaxis()->SetLabelSize(0.);  // the x-axis labels are on the lower pad
    stack.GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());
    stack.GetYaxis()->SetTitleSize(0.055);
    stack.GetYaxis()->SetTitleOffset(1.3);
    stack.GetYaxis()->SetLabelSize(0.045);
    hMCTotal.Draw("E2 same");
    hData.Draw("E1 X0 same");

    TLegend legend(0.62, 0.62, 0.93, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.045);
    legend.AddEntry(&hData, "Data", "pe");
    for (auto it = hMC.rbegin(); it != hMC.rend(); ++it) {  // top of the stack first
        legend.AddEntry(it->get(), (*it)->GetTitle(), "f");
    }
    legend.AddEntry(&hMCTotal, "MC stat. unc.", "f");
    legend.Draw();

    TLatex latex;
    latex.SetNDC();  // coordinates from 0 to 1 instead of axis values
    latex.SetTextFont(62);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.18, 0.83, "CMS");
    latex.SetTextFont(52);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.18, 0.78, "Private Work");
    latex.SetTextFont(42);
    latex.SetTextAlign(31);  // right-aligned
    latex.DrawLatex(0.96, 0.92, lumiLabel);

    // ---- Lower pad: data / MC ratio ----------------------------------------
    lower.cd();
    lower.SetGridy();

    TH1D hRatio(hData);     // points: data / MC, with the data statistical uncertainty
    TH1D hBand(hMCTotal);   // band around 1: the MC statistical uncertainty
    for (int i = 1; i <= hMCTotal.GetNbinsX(); ++i) {
        const double mc = hMCTotal.GetBinContent(i);
        hRatio.SetBinContent(i, mc > 0 ? hData.GetBinContent(i) / mc : 0.);
        hRatio.SetBinError(i, mc > 0 ? hData.GetBinError(i) / mc : 0.);
        hBand.SetBinContent(i, 1.);
        hBand.SetBinError(i, mc > 0 ? hMCTotal.GetBinError(i) / mc : 0.);
    }

    // The band is drawn first, so it defines the axes of the lower pad
    hBand.SetMinimum(0.5);
    hBand.SetMaximum(1.5);
    hBand.GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
    hBand.GetXaxis()->SetTitleSize(0.13);
    hBand.GetXaxis()->SetTitleOffset(1.1);
    hBand.GetXaxis()->SetLabelSize(0.11);
    hBand.GetYaxis()->SetTitle("Data / MC");
    hBand.GetYaxis()->SetTitleSize(0.11);
    hBand.GetYaxis()->SetTitleOffset(0.6);
    hBand.GetYaxis()->SetLabelSize(0.10);
    hBand.GetYaxis()->SetNdivisions(505);
    hBand.GetYaxis()->CenterTitle();
    hBand.Draw("E2");

    TLine one(hBand.GetXaxis()->GetXmin(), 1., hBand.GetXaxis()->GetXmax(), 1.);
    one.SetLineStyle(2);  // dashed
    one.Draw();
    hRatio.Draw("E1 X0 same");

    // ---- Save --------------------------------------------------------------
    canvas.SaveAs((fileName + ".pdf").c_str());
    canvas.SaveAs((fileName + ".png").c_str());
}

void drawMC(const std::vector<TH1D*>& mcStack, const std::string& fileName, bool logY)
{
    std::vector<std::shared_ptr<TH1D>> hMC;  // unique_ptr deletes the copies at the end of the function
    for (const TH1D* h : mcStack) {
        hMC.push_back(std::make_shared<TH1D>(*h));
        addOverflowToEdgeBins(*hMC.back());
        hMC.back()->SetLineColor(h->GetFillColor());
        hMC.back()->SetFillStyle(0);
    }

    // ---- Canvas with two pads: distributions on top, ratio below -----------
    TCanvas canvas("canvas", "", 800, 800);
    TPad upper("upper", "", 0., 0.30, 1., 1.);
    TPad lower("lower", "", 0., 0., 1., 0.30);
    upper.SetMargin(0.15, 0.04, 0.02, 0.10);  // left, right, bottom, top
    lower.SetMargin(0.15, 0.04, 0.35, 0.03);
    upper.Draw();
    lower.Draw();

    // ---- Upper pad: data and stacked MC -----------------------------------
    upper.cd();
    upper.SetLogy(logY);
    
    for (auto h : hMC) {
        h.get()->GetXaxis()->SetLabelSize(0.);  // the x-axis labels are on the lower pad
        h.get()->Draw("Same");
    }

    TLegend legend(0.62, 0.62, 0.93, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.045);
    for (auto it = hMC.rbegin(); it != hMC.rend(); ++it) {  // top of the stack first
        legend.AddEntry(it->get(), (*it)->GetTitle(), "f");
    }
    legend.Draw();

    TLatex latex;
    latex.SetNDC();  // coordinates from 0 to 1 instead of axis values
    latex.SetTextFont(62);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.18, 0.83, "CMS");
    latex.SetTextFont(52);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.18, 0.78, "Private Work");
    latex.SetTextFont(42);
    latex.SetTextAlign(31);  // right-aligned
    latex.DrawLatex(0.96, 0.92, lumiLabel);

    // ---- Lower pad: data / MC ratio ----------------------------------------
    lower.cd();
    lower.SetGridy();

    TH1D hRatio(*hMC.back());     // Zähler: b-Jets
    TH1D hBand(*hMC[0]);          // Nenner/Referenz: non-b-Jets
    for (int i = 1; i <= hMC[0]->GetNbinsX(); ++i) {
        const double non_b = hMC[0]->GetBinContent(i);
        const double b = hMC.back()->GetBinContent(i);
        
        hRatio.SetBinContent(i, non_b > 0 ? b / non_b : 0.);
        hRatio.SetBinError(i, 0.); 
        
        hBand.SetBinContent(i, 1.);
        hBand.SetBinError(i, 0.);
    }

    // The band is drawn first, so it defines the axes of the lower pad
    hBand.SetMinimum(0.0);
    hBand.SetMaximum(5.0);
    hBand.GetXaxis()->SetTitle(hMC[0]->GetXaxis()->GetTitle());
    hBand.GetXaxis()->SetTitleSize(0.13);
    hBand.GetXaxis()->SetTitleOffset(1.1);
    hBand.GetXaxis()->SetLabelSize(0.11);
    hBand.GetYaxis()->SetTitle("b / non-b");
    hBand.GetYaxis()->SetTitleSize(0.11);
    hBand.GetYaxis()->SetTitleOffset(0.6);
    hBand.GetYaxis()->SetLabelSize(0.10);
    hBand.GetYaxis()->SetNdivisions(505);
    hBand.GetYaxis()->CenterTitle();
    hBand.Draw("E2");

    TLine one(hBand.GetXaxis()->GetXmin(), 1., hBand.GetXaxis()->GetXmax(), 1.);
    one.SetLineStyle(2);  // dashed
    one.Draw();
    hRatio.Draw("E1 X0 same");
    // ---- Save --------------------------------------------------------------
    canvas.SaveAs((fileName + ".pdf").c_str());
    canvas.SaveAs((fileName + ".png").c_str());
}

void drawMC2D(TH2D* h, const std::string& fileName)
{
    auto c = new TCanvas("", "", 0, 0, 1100, 550);

    c->SetLogz(1);

    std::cout << h->GetCorrelationFactor() << std::endl;

    h->Draw("colz");

    // ---- Save --------------------------------------------------------------
    c->SaveAs((fileName + ".pdf").c_str());
    c->SaveAs((fileName + ".png").c_str());
}

}  // namespace plotting
