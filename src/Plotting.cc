#include "Plotting.h"

#include "Config.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>

namespace plotting {

void setStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetEndErrorSize(0);
}

// Return a copy of the histogram with underflow/overflow moved into the first/last bin.
static TH1D* copyWithOverflow(const TH1D* h, const std::string& name)
{
    TH1D* copy = static_cast<TH1D*>(h->Clone(name.c_str()));
    copy->SetDirectory(nullptr);
    const int n = copy->GetNbinsX();
    copy->SetBinContent(1, h->GetBinContent(0) + h->GetBinContent(1));
    copy->SetBinError(1, std::hypot(h->GetBinError(0), h->GetBinError(1)));
    copy->SetBinContent(n, h->GetBinContent(n) + h->GetBinContent(n + 1));
    copy->SetBinError(n, std::hypot(h->GetBinError(n), h->GetBinError(n + 1)));
    copy->SetBinContent(0, 0.);
    copy->SetBinContent(n + 1, 0.);
    return copy;
}

void drawDataMC(const TH1D* data, const std::vector<TH1D*>& mcByFlavour, const std::string& xTitle, bool logY,
                const std::string& fileName)
{
    // We draw copies, so the original histograms stay unchanged.
    TH1D* hData = copyWithOverflow(data, "hData");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.8);
    hData->SetLineColor(kBlack);

    THStack* stack = new THStack("stack", "");
    std::vector<TH1D*> hMC;
    TH1D* hMCTotal = nullptr;
    for (size_t i = 0; i < mcByFlavour.size(); ++i) {
        TH1D* h = copyWithOverflow(mcByFlavour[i], "hMC_" + config::flavours[i].name);
        h->SetFillColor(config::flavours[i].color);
        h->SetLineColor(config::flavours[i].color);
        stack->Add(h);
        hMC.push_back(h);

        if (hMCTotal == nullptr) {
            hMCTotal = static_cast<TH1D*>(h->Clone("hMCTotal"));
        } else {
            hMCTotal->Add(h);
        }
    }

    TCanvas* canvas = new TCanvas("canvas", "", 800, 800);

    // ---- upper pad: distributions ------------------------------------------
    TPad* upper = new TPad("upper", "", 0., 0.30, 1., 1.);
    upper->SetBottomMargin(0.02);
    upper->SetLeftMargin(0.15);
    upper->SetRightMargin(0.04);
    upper->SetLogy(logY);
    upper->Draw();
    upper->cd();

    const double yMax = std::max(hData->GetMaximum(), hMCTotal->GetMaximum());
    stack->SetMaximum(logY ? yMax * 500. : yMax * 1.5);
    stack->SetMinimum(logY ? 0.5 : 0.);
    stack->Draw("hist");
    stack->GetXaxis()->SetLabelSize(0.);
    stack->GetYaxis()->SetTitle("Jets");
    stack->GetYaxis()->SetTitleSize(0.055);
    stack->GetYaxis()->SetTitleOffset(1.3);
    stack->GetYaxis()->SetLabelSize(0.045);

    // MC statistical uncertainty as hatched band
    hMCTotal->SetFillStyle(3354);
    hMCTotal->SetFillColor(kGray + 2);
    hMCTotal->SetLineWidth(0);
    hMCTotal->SetMarkerSize(0);
    hMCTotal->Draw("E2 same");

    hData->Draw("E1 X0 same");

    TLegend* legend = new TLegend(0.62, 0.62, 0.93, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.045);
    legend->AddEntry(hData, "Data", "pe");
    // reverse order, so the legend matches the stack from top to bottom
    for (int i = (int)hMC.size() - 1; i >= 0; --i) {
        legend->AddEntry(hMC[i], ("t#bar{t}, " + config::flavours[i].label).c_str(), "f");
    }
    legend->AddEntry(hMCTotal, "MC stat. unc.", "f");
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.055);
    latex.SetTextFont(62);
    latex.DrawLatex(0.18, 0.83, "CMS");
    latex.SetTextFont(52);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.18, 0.78, "Private Work");
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.96, 0.92, "109 fb^{-1} (13.6 TeV)");

    // ---- lower pad: data / MC ratio ----------------------------------------
    canvas->cd();
    TPad* lower = new TPad("lower", "", 0., 0., 1., 0.30);
    lower->SetTopMargin(0.03);
    lower->SetBottomMargin(0.35);
    lower->SetLeftMargin(0.15);
    lower->SetRightMargin(0.04);
    lower->SetGridy();
    lower->Draw();
    lower->cd();

    // ratio = data / MC (data statistical uncertainty only)
    TH1D* hRatio = static_cast<TH1D*>(hData->Clone("hRatio"));
    // MC statistical uncertainty band around 1
    TH1D* hBand = static_cast<TH1D*>(hMCTotal->Clone("hBand"));
    for (int i = 1; i <= hMCTotal->GetNbinsX(); ++i) {
        const double mc = hMCTotal->GetBinContent(i);
        hRatio->SetBinContent(i, mc > 0 ? hData->GetBinContent(i) / mc : 0.);
        hRatio->SetBinError(i, mc > 0 ? hData->GetBinError(i) / mc : 0.);
        hBand->SetBinContent(i, 1.);
        hBand->SetBinError(i, mc > 0 ? hMCTotal->GetBinError(i) / mc : 0.);
    }

    hBand->SetMinimum(0.5);
    hBand->SetMaximum(1.5);
    hBand->GetXaxis()->SetTitle(xTitle.c_str());
    hBand->GetXaxis()->SetTitleSize(0.13);
    hBand->GetXaxis()->SetTitleOffset(1.1);
    hBand->GetXaxis()->SetLabelSize(0.11);
    hBand->GetYaxis()->SetTitle("Data / MC");
    hBand->GetYaxis()->SetTitleSize(0.11);
    hBand->GetYaxis()->SetTitleOffset(0.6);
    hBand->GetYaxis()->SetLabelSize(0.10);
    hBand->GetYaxis()->SetNdivisions(505);
    hBand->GetYaxis()->CenterTitle();
    hBand->Draw("E2");

    TLine* one = new TLine(hBand->GetXaxis()->GetXmin(), 1., hBand->GetXaxis()->GetXmax(), 1.);
    one->SetLineStyle(2);
    one->Draw();
    hRatio->Draw("E1 X0 same");

    canvas->SaveAs((fileName + ".pdf").c_str());
    canvas->SaveAs((fileName + ".png").c_str());

    // Deleting the canvas also deletes the pads; the histograms, stack, legend
    // and line we created with new/Clone must be deleted by us.
    delete canvas;
    delete stack;
    delete legend;
    delete one;
    delete hData;
    delete hMCTotal;
    delete hRatio;
    delete hBand;
    for (TH1D* h : hMC) delete h;
}

}  // namespace plotting
