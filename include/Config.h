#pragma once
// ============================================================================
//  Config.h  --  input files and MC jet flavour categories
// ============================================================================

#include <Rtypes.h>  // ROOT colour constants (kRed, kAzure, ...)

#include <string>
#include <vector>

namespace config {

// All files in one list are chained together and processed as one sample.
inline const std::vector<std::string> dataFiles = {
    "data/MuonEGRun2024C_BDT.root",
    "data/MuonEGRun2024D_BDT.root",
    "data/MuonEGRun2024E_BDT.root",
    "data/MuonEGRun2024F_BDT.root",
    "data/MuonEGRun2024G_BDT.root",
    "data/MuonEGRun2024H_BDT.root",
    "data/MuonEGRun2024I_BDT.root",
};

inline const std::vector<std::string> ttbarFiles = {
    "data/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_None_BDT.root",
};

// MC jet flavour categories, using the hadron flavour (5 = b, 4 = c, 0 = light).
// They are stacked in the plots in this order (first = bottom of the stack).
struct FlavourCategory {
    std::string name;       // used in column and histogram names
    std::string label;      // legend entry
    std::string selection;  // condition on the per-jet column Jet_flavour
    int color;
};

inline const std::vector<FlavourCategory> flavours = {
    {"light", "light jets", "Jet_flavour != 5 && Jet_flavour != 4", kAzure + 1},
    {"c",     "c jets",     "Jet_flavour == 4",                     kGreen + 2},
    {"b",     "b jets",     "Jet_flavour == 5",                     kRed + 1},
};

}  // namespace config
