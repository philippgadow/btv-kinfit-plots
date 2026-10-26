# UParT b-tag data/MC comparison in ttbar → eμ

This code makes data/MC comparison plots of jet variables: the UParT b-tag
discriminant, jet pT and jet η. It uses ROOT's
[RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) and C++.
The MC is normalised to 109 fb⁻¹, split into b, c and light jets, and drawn as
a stacked histogram.

## Quick start

ROOT on this machine is installed as a *snap* package. The snap comes with its
own compiler, so every command has to run inside the snap environment, which
`root-framework.run` gives you:

```bash
root-framework.run make          # compile -> bin/makePlots
root-framework.run make run      # compile if needed, then run (~30 s)
root-framework.run make clean    # delete compiled files
```

Output:
- `output/plots/<variable>.pdf/.png`: the data/MC plots
- `output/histograms.root`: all histograms

> Snap pitfall: snap ROOT can't read files under `/tmp`. Keep all your files in
> your home directory.

## Code layout

| File | What it does |
|---|---|
| `src/makePlots.cc` | **Start here.** The whole analysis: event and jet selection, variables, MC normalisation, histograms |
| `include/Config.h` | Lists of input files and the MC flavour categories (b, c, light) |
| `include/Plotting.h`, `src/Plotting.cc` | Draws the stacked data/MC plot with the ratio panel. You don't need to read this to get started |
| `Makefile` | Compiles every `src/*.cc` file and links them into `bin/makePlots` |

## How the input files are organised (read this!)

- There is one tree, `Events`. **One entry is one event**, not one jet.
- Per-jet quantities are stored as arrays inside each event:
  - `jetpt` and `jeteta` have exactly `njetpt` entries.
  - The other per-jet branches (`btagUParTAK4B`, `flavour`, `jetrank`,
    `close_*`, `far_*`, `j2ll_*`, ...) always have 8 entries in MC and 7 in data.
- ⚠️ **Only the first `njetpt` entries of these arrays belong to the event.**
  The entries after them are `NaN` or values left over from other events.
  `goodJetIndices()` in `makePlots.cc` returns the positions of the jets to use.
  To get a per-jet variable, always write `Take(<branch>, goodJets)`.
- `flavour` is the hadron flavour (5 = b, 4 = c, 0 = light). It is 0 for all
  jets in data.
- MC `weight` = generator weight × pileup weight × ttbar weight. The code
  multiplies it by σ(tt̄) × L / (sum of generator weights).

To look at a file yourself:

```bash
rootls -t data/MuonEGRun2024C_BDT.root     # list all branches
root -l data/MuonEGRun2024C_BDT.root
root [1] ROOT::RDataFrame df("Events", "data/MuonEGRun2024C_BDT.root");
root [2] df.Display({"njetpt", "jetpt", "btagUParTAK4B", "close_mlj"}, 3)->Print()
```

## Task 1: add the close_*, far_* and j2ll_* variables

Each jet has these extra variables:

- `close_mlj`, `close_deta`, `close_dphi`, `close_ptrel`, `close_lj2ll_deta`, `close_lj2ll_dphi`
- `far_mlj`, `far_deta`, `far_dphi`, `far_ptrel`, `far_lj2ll_deta`, `far_lj2ll_dphi`
- `j2ll_deta`, `j2ll_dphi`

Make a data/MC plot for each of them. In `src/makePlots.cc` there are two
`TODO (task 1)` comments:

1. In `selectEvents()`, make a per-jet column, the same way as `Jet_btagUParTB`.
2. Add an entry to the `variables` list with a sensible axis title and binning.
   To choose the binning, first look at the values with `Display` (see above),
   or quickly draw them in the ROOT prompt.

Before you plot, find out what each variable means. What are "close" and
"far", and what is `ptrel`?

## Ideas for later

- Add the other UParT outputs (`btagUParTAK4CvL`, `btagUParTAK4CvB`, ...).
- Plot only the leading jet (`Jet_pt[0]`), or each `jetrank` separately.
- Add event-level plots (number of good jets, `npvn`).
- Systematic uncertainty bands from the `weight_syst_*` weights.
- Other backgrounds (single top, DY, ...) as additional MC samples.
