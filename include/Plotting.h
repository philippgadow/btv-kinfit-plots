#pragma once
// ============================================================================
//  Plotting.h  --  drawing data/MC comparison plots
// ============================================================================

#include <TH1D.h>
#include <TH2D.h>

#include <string>
#include <vector>

namespace plotting {

// Global plot style (fonts, no statistics box, ...). Call once at the start.
void setStyle();

// Draw data (black points) on top of the stacked MC histograms, with a
// data/MC ratio panel below, and save it as <fileName>.pdf and <fileName>.png.
//
//  - mcStack: the MC histograms, bottom of the stack first. Each one is drawn
//    with its fill colour and appears in the legend with its title.
//  - The x-axis title is taken from the data histogram.
//  - Underflow and overflow are added to the first and last bin.
//
// The histograms you pass in are not modified.
void drawDataMC(const TH1D* data, const std::vector<TH1D*>& mcStack, const std::string& fileName,
                bool logY = false);

void drawMC(const std::vector<TH1D*>& mcStack, const std::string& fileName,
                bool logY = false);

void drawMC2D(TH2D* h, const std::string& fileName);

}  // namespace plotting
