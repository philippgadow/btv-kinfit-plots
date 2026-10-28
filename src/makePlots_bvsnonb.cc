// ============================================================================
//  makePlots.cc  --  data/MC comparison of jet variables in ttbar -> e mu events
//
//    Settings                        input files, variables to plot, MC normalisation
//    Step 1  defineJets()            OBJECT DEFINITION: which jets of an event do we use?
//    Step 2  selectEvents()          EVENT SELECTION:   which events do we keep?
//    Step 3  bookHistograms()        HISTOGRAMMING:     which distributions do we fill?
//    Step 4  main()                  RUN AND OUTPUT:    run the event loop, write the
//                                                       histograms and draw the plots
//
//  Run with:  root-framework.run make run
// ============================================================================

#include "Plotting.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <Rtypes.h>  // ROOT colour constants (kRed, kAzure, ...)
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ============================================================================
//  Settings
// ============================================================================

// Input files
// All files in one list are chained together and processed as one sample

// ttbar MC
const std::vector<std::string> ttbarFiles = {
    "data/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_None_BDT.root",
};

// A variable to plot: the name of a per-jet column (see Step 1), the x-axis title and the binning
struct Variable {
    std::string column;
    std::string xTitle;
    int nBins;
    double xMin;
    double xMax;
    bool logY;
};

const std::vector<Variable> variables = {
    {"Jet_btagUParTB", "UParT AK4 b-tag discriminant", 50, 0.0, 1.0, false},
    {"Jet_pt", "Jet p_{T} [GeV]", 50, 0.0, 500.0, false},
    {"Jet_eta", "Jet #eta", 50, -2.5, 2.5, false},
    {"Jet_close_mlj", "close m_{lj} [GeV]", 50, 0.0, 500, false},
    {"Jet_close_deta", "close d_{eta} ", 50, 0.0, 5, false},
    {"Jet_close_dphi", "close d_{phi} ", 50, 0, 3.14, false},
    {"Jet_close_ptrel", "close pt_{rel} [GeV]", 50, 0.0, 500, false},
    {"Jet_close_lj2ll_deta", "close lj2ll d_{eta} ", 50, 0.0, 5, false},
    {"Jet_close_lj2ll_dphi", "close lj2ll d_{phi} ", 50, 0, 3.14, false},
    {"Jet_far_mlj", "far far m_{lj} [GeV] ", 50, 0.0, 500, false},
    {"Jet_far_deta", "far far d_{eta} ", 50, 0.0, 5, false},
    {"Jet_far_dphi", "far far d_{phi} ", 50, 0, 3.14, false},
    {"Jet_far_ptrel", "far far pt_{rel} [GeV]", 50, 0.0, 500, false},
    {"Jet_far_lj2ll_deta", "far lj2ll d_{eta} ", 50, 0.0, 5, false},
    {"Jet_far_lj2ll_dphi", "far lj2ll d_{phi} ", 50, 0, 3.14, false},
    {"Jet_j2ll_deta", "j2ll d_{eta} ", 50, 0.0, 6, false},
    {"Jet_j2ll_dphi", "j2ll d_{phi} ", 50, 0, 3.14, false},
    {"Jet_kindisc", "kindisc", 50, -1, 1, false}
};

// MC jet flavour categories, using the hadron flavour (5 = b, 4 = c, 0 = light).
// Each category gets its own histogram. They are stacked in the plots in this
// order (first = bottom of the stack).
struct FlavourCategory {
    std::string name;       // used in histogram names
    std::string label;      // legend entry
    std::string selection;  // condition on the per-jet column Jet_flavour
    int color;
};

const std::vector<FlavourCategory> flavours = {
    {"non_b", "t#bar{t}, non-b jets", "Jet_flavour != 5", kAzure + 1},
    {"b",     "t#bar{t}, b jets",     "Jet_flavour == 5", kRed + 1}
};

const std::string outDir = "output";

// ============================================================================
//  Step 1: Object definition -- which jets do we use?
//
//  One entry of the tree is one *event*. The per-jet quantities are arrays:
//  jetpt and jeteta have exactly njetpt entries, but the other per-jet arrays
//  (btagUParTAK4B, flavour, jetrank, close_*, far_*, ...) always have a fixed
// number of entries. Only the first njetpt entries belong to the event!
//  The rest are padding or left-overs from other events.
//
//  So we first find the positions (indices) of the "good jets" of the event:
//  jets that belong to the event and pass the jet selection. Then we build
//  new columns "Jet_*" that contain only these jets.
// ============================================================================

// Returns the indices of the good jets of one event.
ROOT::RVecI goodJetIndices(int njet, const ROOT::RVecD& jetrank, const ROOT::RVecD& jetpt,
                           const ROOT::RVecD& jeteta)
{
    ROOT::RVecI indices;
    for (int i = 0; i < njet && i < (int)jetrank.size(); ++i) {
        if (jetrank[i] != i) continue;  // corrupted entry (happens in a few data events)
        if (jetpt[i] < 30. || std::abs(jeteta[i]) > 2.5) continue;  // jet selection: pT > 30 GeV, |eta| < 2.5
        indices.push_back(i);
    }
    return indices;
}

// Adds the jet columns to the dataframe. The same for data and MC.
//
// After this, every Jet_* column is an array with one entry per good jet,
// so e.g. Jet_pt[0] and Jet_btagUParTB[0] belong to the same jet.
// "Take(jetpt, goodJets)" picks the entries of jetpt at the positions in goodJets.
ROOT::RDF::RNode defineJets(ROOT::RDF::RNode df)
{
    return df.Define("goodJets", goodJetIndices, {"njetpt", "jetrank", "jetpt", "jeteta"})
        .Define("nGoodJets", "goodJets.size()")
        .Define("Jet_pt", "Take(jetpt, goodJets)")
        .Define("Jet_eta", "Take(jeteta, goodJets)")
        .Define("Jet_flavour", "Take(flavour, goodJets)")
        .Define("Jet_btagUParTB", "Take(btagUParTAK4B, goodJets)")
        .Define("Jet_close_mlj", "Take(close_mlj, goodJets)")
        .Define("Jet_close_deta", "Take(close_deta, goodJets)")
        .Define("Jet_close_dphi", "Take(close_dphi, goodJets)")
        .Define("Jet_close_ptrel", "Take(close_ptrel, goodJets)")
        .Define("Jet_close_lj2ll_deta", "Take(close_lj2ll_deta, goodJets)")
        .Define("Jet_close_lj2ll_dphi", "Take(close_lj2ll_dphi, goodJets)")
        .Define("Jet_far_mlj", "Take(far_mlj, goodJets)")
        .Define("Jet_far_deta", "Take(far_deta, goodJets)")
        .Define("Jet_far_dphi", "Take(far_dphi, goodJets)")
        .Define("Jet_far_ptrel", "Take(far_ptrel, goodJets)")
        .Define("Jet_far_lj2ll_deta", "Take(far_lj2ll_deta, goodJets)")
        .Define("Jet_far_lj2ll_dphi", "Take(far_lj2ll_dphi, goodJets)")
        .Define("Jet_j2ll_deta", "Take(j2ll_deta, goodJets)")
        .Define("Jet_j2ll_dphi", "Take(j2ll_dphi, goodJets)")
        .Define("Jet_kindisc", "Take(kindisc, goodJets)");
}

// ============================================================================
//  Step 2: Event selection -- which events do we keep?
//
//  Filter() keeps only the events for which the condition is true. The name
//  of the cut (second argument) is printed in the cut flow report at the end.
//  More cuts can be added by chaining more .Filter(...) calls.
// ============================================================================
ROOT::RDF::RNode selectEvents(ROOT::RDF::RNode df)
{
    return df.Filter("nGoodJets >= 2", "at least two good jets");
}

// ============================================================================
//  Step 3: Histogramming -- which distributions do we fill?
//
//  For every variable we fill one data histogram and one MC histogram per
//  flavour category. Filling an *array* column (e.g. Jet_pt) puts every jet
//  of the event into the histogram; in MC every jet gets the event weight.
//
//  Histo1D does not fill anything yet! RDataFrame only runs the event loop
//  when a result is first used (in Step 4), and then fills all booked
//  histograms of a sample in one go. That is why Histo1D returns a
//  RResultPtr<TH1D> and not a TH1D.
// ============================================================================

// All histograms of one variable
struct Histograms {
    ROOT::RDF::RResultPtr<TH1D> data;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> mcByFlavour;  // same order as 'flavours'
};

Histograms bookHistograms(ROOT::RDF::RNode ttbar, const Variable& var)
{
    Histograms hists;
    const std::string axisTitles = ";" + var.xTitle + ";Fraction of jets";

    // MC: one histogram per flavour, filled only with the jets of that flavour.
    for (const auto& flav : flavours) {
        const std::string column = var.column + "_" + flav.name;  // e.g. Jet_pt_b
        const std::string mcName = "ttbar_" + column;
        auto jetsOfThisFlavour = ttbar.Define(column, var.column + "[" + flav.selection + "]");
        hists.mcByFlavour.push_back(jetsOfThisFlavour.Histo1D(
            {mcName.c_str(), axisTitles.c_str(), var.nBins, var.xMin, var.xMax}, column, "mcWeight"));
    }
    return hists;
}

// ============================================================================
//  Step 4: Run everything and write the output
// ============================================================================
int main()
{
    // ---- Setup -------------------------------------------------------------------
    gROOT->SetBatch(true);  // don't open graphics windows, just write files
    plotting::setStyle();
    gSystem->mkdir((outDir + "/plots").c_str(), true);
    ROOT::EnableImplicitMT();  // use all CPU cores

    // ---- Open the files ------------------------------------------------------------
    ROOT::RDataFrame ttbarFrame("Events", ttbarFiles);
    ROOT::RDF::Experimental::AddProgressBar(ttbarFrame);  // the big MC file takes a moment

    // ---- Step 1 + 2: define the jets, then select the events ---------------------------
    ROOT::RDF::RNode ttbar = selectEvents(defineJets(ttbarFrame));

    // MC only: the event weight. 'weight' = generator weight * pileup weight * ttbar weight,
    // and we scale it to the data luminosity.
    const double mcScale = 1;
    ttbar = ttbar.Define("mcWeight", [mcScale](double weight) { return weight * mcScale; }, {"weight"});

    // ---- Step 3: book the histograms -----------------------------------------------------
    std::vector<Histograms> allHists;  // one entry per variable
    for (const auto& var : variables) {
        allHists.push_back(bookHistograms(ttbar, var));
    }
    auto nMCEvents = ttbar.Sum("mcWeight");
    auto ttbarReport = ttbar.Report();

    // ---- Step 4: run the event loop and write everything out ------------------------------
    // Using a result (with * or ->) starts the event loop of that sample, once.
    std::string targetVar = "Jet_kindisc";

    TFile outFile((outDir + "/histograms.root").c_str(), "RECREATE");
    for (size_t i = 0; i < variables.size(); ++i) {
        const Variable& var = variables[i];
        Histograms& hists = allHists[i];

        std::vector<TH1D*> mcStack;

        for (size_t j = 0; j < flavours.size(); ++j) {
            TH1D* h = hists.mcByFlavour[j].GetPtr();

            h->Scale(1.0 / h->Integral()); 

            h->Write();
            h->SetFillColor(flavours[j].color);
            h->SetTitle(flavours[j].label.c_str());
            mcStack.push_back(h);
        }

        plotting::drawMC(mcStack, outDir + "/plots/" + var.column, var.logY);
    }
    
    std::cout << "\nHistograms written to " << outFile.GetName() << ", plots to " << outDir << "/plots/\n";
    return 0;
}