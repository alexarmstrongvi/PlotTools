PlotTools
=============

Classes and tools for making ROOT plots from flat ntuples

# Contents

* [Introduction](#introduction)
* [Setup](#setup)
* [How the sausage is made](#how-the-sausage-is-made)

# Introduction
Making plots with ROOT can be a mess.
The way the actual histogram content of TH1 (e.g. bin values and edges) is integrated with a TH1 plot (e.g. line color, label size, fill style, axis labels, etc...) is not optimal.
Moreover, this is the opposite for TGraphs which do separate the content from the presentation.

In the process of making plotting scripts, one often has a lot of duplicate code.
For example, one might have code for making and formatting a 1D plot duplicated in one script for making stacked histograms and another for making a ratio plot.
This package attempts to modularize the process of making plots as much as possible (still a work in progress) so that duplicate code is minimized and creating new types of plots from simpler building blocks is quick and easy.

# Setup
The package can be setup anywhere and then added to PYTHONPATH to make the classes accessible for python scripts.
I prefer to add it to the home directory 
```bash
cd ${DIR_OF_CHOICE}
git clone git@github.com:alexarmstrongvi/PlotTools.git
source ./PlotTools/bash/add_to_python_path.sh
```

Make sure to source [`add_to_python_path.sh`](bash/add_to_python_path.sh) whenever starting in a new environment.

# How the sausage is made
Everything centers around running TTree::Draw on flatntuples to make plots that get combined into any number of predefined plot types (e.g. MC background stack overlaid with data in the top panel and a Data/MC ratio plot in the bottom panel)

All the classes and tools facilitate writing looper scripts that generally follow the procedure below:
1) Load in specified flat ntuples
2) Group the ntuples into TChains corresponding to specified samples (e.g. data, diboson, fake background, etc...)
3) Loop over specified regions and build TCuts for each sample
4) For each region, loop over specified plots that carry out the following:
    1) Use TTree::Draw to build all the inital histogram objects from the samples, applying the intended TCut
    2) Group histograms as needed for the specified plot
    3) Create and format the TCanvas
    4) Add all needed plots onto the TCanvas and save as an image

While this is the most common type of script to write, it is not the only one.
Below are examples of scripts that use these classes to make plots, yield tables, fake factors, etc...:
- TODO

TODO: Describe use of loopers and config files
