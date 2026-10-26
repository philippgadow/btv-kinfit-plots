#pragma once
// ============================================================================
//  Plotting.h  --  drawing data/MC comparison plots
// ============================================================================

#include <TH1D.h>

#include <string>
#include <vector>

namespace plotting {

// Global plot style (fonts, no statistics box, ...). Call once at the start.
void setStyle();

// Draw data (black points) on top of the stacked MC histograms (one per
// flavour category, same order as config::flavours), with a data/MC ratio
// panel below. Underflow and overflow are added to the first and last bin.
// Saves <fileName>.pdf and <fileName>.png
void drawDataMC(const TH1D* data, const std::vector<TH1D*>& mcByFlavour, const std::string& xTitle, bool logY,
                const std::string& fileName);

}  // namespace plotting
