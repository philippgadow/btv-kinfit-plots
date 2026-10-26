// ============================================================================
//  makePlots.cc  --  data/MC comparison of jet variables in ttbar -> e mu events
//
//  1. open the data and ttbar MC files with RDataFrame
//  2. select events and jets, and make one column per jet variable ("Jet_*")
//  3. fill histograms (MC split into b, c and light jets)
//  4. draw data/MC plots and save the histograms
//
//  Run with:  root-framework.run make run
// ============================================================================

#include "Config.h"
#include "Plotting.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
//  Which jets do we use?
//
//  One entry of the tree is one *event*. The per-jet quantities are arrays:
//  jetpt and jeteta have exactly njetpt entries, but the other per-jet arrays
//  (btagUParTAK4B, flavour, jetrank, close_*, far_*, ...) always have 8 (MC)
//  or 7 (data) entries. Only the first njetpt entries belong to the event!
//  The rest are padding or left-overs from other events.
//
//  This function returns the positions (indices) of the jets we want to use:
//  jets that belong to the event and pass the jet selection.
// ----------------------------------------------------------------------------
ROOT::RVecI goodJetIndices(int njet, const ROOT::RVecD& jetrank, const ROOT::RVecD& jetpt,
                           const ROOT::RVecD& jeteta)
{
    ROOT::RVecI indices;
    for (int i = 0; i < njet && i < (int)jetrank.size(); ++i) {
        if (jetrank[i] != i) continue;  // corrupted entry (happens in a few data events)
        if (jetpt[i] < 30. || std::abs(jeteta[i]) > 2.5) continue;  // jet selection
        indices.push_back(i);
    }
    return indices;
}

// ----------------------------------------------------------------------------
//  Event selection and per-jet columns. The same for data and MC.
//
//  After this, every Jet_* column is an array with one entry per good jet,
//  so e.g. Jet_pt[0] and Jet_btagUParTB[0] belong to the same jet.
// ----------------------------------------------------------------------------
ROOT::RDF::RNode selectEvents(ROOT::RDF::RNode df)
{
    return df.Define("goodJets", goodJetIndices, {"njetpt", "jetrank", "jetpt", "jeteta"})
        .Filter("goodJets.size() >= 2", "at least two good jets")
        .Define("Jet_pt", "Take(jetpt, goodJets)")
        .Define("Jet_eta", "Take(jeteta, goodJets)")
        .Define("Jet_flavour", "Take(flavour, goodJets)")
        .Define("Jet_btagUParTB", "Take(btagUParTAK4B, goodJets)");
    // TODO (task 1): add the close_*, far_* and j2ll_* variables in the same way,
    //                e.g. .Define("Jet_close_mlj", "Take(close_mlj, goodJets)")
}

// A variable to plot: a Jet_* column, the x-axis title and the binning
struct Variable {
    std::string column;
    std::string xTitle;
    int nBins;
    double xMin;
    double xMax;
    bool logY;
};

int main()
{
    // ---- Settings ----------------------------------------------------------------
    const std::vector<Variable> variables = {
        {"Jet_btagUParTB", "UParT AK4 b-tag discriminant", 50, 0.0, 1.0, false},
        {"Jet_pt", "Jet p_{T} [GeV]", 50, 0.0, 500.0, false},
        {"Jet_eta", "Jet #eta", 50, -2.5, 2.5, false},
        // TODO (task 1): add plots for the close_*, far_* and j2ll_* variables
    };

    // MC normalisation: scale every MC event by  cross section * luminosity / sum of generator weights
    const double luminosity = 109000.;               // 109 fb^-1 = 109000 pb^-1
    const double ttbarCrossSection = 98.03;          // pb: sigma(tt) = 923.6 pb times BR(tt -> 2l2nu) = (3 x 0.1086)^2
    const double ttbarSumGenWeights = 29473520581.6875;  // sum of genweight of ALL generated events
    const double mcScale = ttbarCrossSection * luminosity / ttbarSumGenWeights;

    const std::string outDir = "output";

    // ---- Setup -------------------------------------------------------------------
    gROOT->SetBatch(true);  // don't open graphics windows, just write files
    plotting::setStyle();
    gSystem->mkdir((outDir + "/plots").c_str(), true);
    ROOT::EnableImplicitMT();  // use all CPU cores

    // ---- Open the files --------------------------------------------------------------
    // Passing a list of files makes RDataFrame chain them together.
    ROOT::RDataFrame dataFrame("Events", config::dataFiles);
    ROOT::RDataFrame ttbarFrame("Events", config::ttbarFiles);
    ROOT::RDF::Experimental::AddProgressBar(ttbarFrame);  // the big MC file takes a moment

    ROOT::RDF::RNode data = selectEvents(dataFrame);
    ROOT::RDF::RNode ttbar = selectEvents(ttbarFrame);

    // Event weight for MC. 'weight' = generator weight * pileup weight * ttbar weight.
    ttbar = ttbar.Define("mcWeight", [mcScale](double w) { return w * mcScale; }, {"weight"});

    // For MC, make one column per variable and flavour, containing only the jets
    // of that flavour, e.g. "Jet_pt_b" = pT of all b jets in the event.
    for (const auto& var : variables) {
        for (const auto& flav : config::flavours) {
            ttbar = ttbar.Define(var.column + "_" + flav.name, var.column + "[" + flav.selection + "]");
        }
    }

    // ---- Book histograms -----------------------------------------------------------
    // Filling an array column puts every jet of the event into the histogram;
    // in MC each jet gets the event weight.
    // Nothing is read yet: RDataFrame only runs the event loop when a result is first
    // used, and then fills all booked histograms of that sample in one go.
    std::vector<ROOT::RDF::RResultPtr<TH1D>> dataHists;              // [variable]
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> mcHists;  // [variable][flavour]
    for (const auto& var : variables) {
        const std::string axes = ";" + var.xTitle + ";Jets";
        const std::string name = "data_" + var.column;
        dataHists.push_back(data.Histo1D({name.c_str(), axes.c_str(), var.nBins, var.xMin, var.xMax}, var.column));

        std::vector<ROOT::RDF::RResultPtr<TH1D>> byFlavour;
        for (const auto& flav : config::flavours) {
            const std::string column = var.column + "_" + flav.name;
            const std::string mcName = "ttbar_" + column;
            byFlavour.push_back(
                ttbar.Histo1D({mcName.c_str(), axes.c_str(), var.nBins, var.xMin, var.xMax}, column, "mcWeight"));
        }
        mcHists.push_back(byFlavour);
    }
    auto nDataEvents = data.Count();
    auto nMCEvents = ttbar.Sum("mcWeight");

    // ---- Results -----------------------------------------------------------------
    // Using a result with * starts the event loop of that sample (if not done yet).
    const double nData = *nDataEvents;
    const double nMC = *nMCEvents;
    std::cout << "\nSelected events:  data " << std::llround(nData) << "   ttbar MC " << std::llround(nMC)
              << "   data/MC " << std::setprecision(3) << nData / nMC << "\n\n";

    TFile outFile((outDir + "/histograms.root").c_str(), "RECREATE");
    for (size_t i = 0; i < variables.size(); ++i) {
        dataHists[i]->Write();
        std::vector<TH1D*> mcByFlavour;
        for (auto& h : mcHists[i]) {
            h->Write();
            mcByFlavour.push_back(h.GetPtr());
        }

        const std::string plotFile = outDir + "/plots/" + variables[i].column;
        plotting::drawDataMC(dataHists[i].GetPtr(), mcByFlavour, variables[i].xTitle, variables[i].logY, plotFile);
    }

    std::cout << "\nHistograms written to " << outFile.GetName() << ", plots to " << outDir << "/plots/\n";
    return 0;
}
