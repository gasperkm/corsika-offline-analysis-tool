_# corsika-offline-analysis-tool
Analysis tool for CORSIKA + Offline simulations.

# information on installation
Before installing, the Offline software must be installed and their directory paths updated in the offline_setup.sh.in bash script and save it as offline_setup.sh.
Before running for the first time, create folders input and results in the base directory.

# Program usage
For program usage, run "make usage". Each of the subprograms has a separate function:
- analysis_tool: Has a general statistics code that plots a histogram of any observable or a scatter plot of two observables (-f and -t options). In addition, it can also perform a multivariate analysis (MVA), with the use of the ROOT package TMVA (-m and -mg options).
- histogram_replot: Any data from histograms or scatter plots is written in .dat files, which can be used to replot them with this subprogram.
- hadd: Adds two or more ROOT files with the same structure to a single file. Useful for gathering ADST analysis files from Offline.
- tmvagui: Runs the GUI subprogram for viewing results of MVA training and testing on any file (default is tmva_output.root).

# Skipping standard input values
To make things quicker, it is possible to use a standard input file instead of writing into the terminal every time. However, all needed stdin values need to be supplied, otherwise the program will most likely fail or go into an endless loop.

#--- Example 1 -------------------------------------------------------------------
Example: ./analysis_tool -m < infile
Input file (infile):
xmax,shfoot,shwsize,risetime
1
1
MLPBNN
0.52

The above values are prompted by:
1. Select the observables to include in the MVA (comma separated):
2. Select one to analyze (1-2):
3. Take all events as background (0) or only inverse events (1)?
4. Select one of the above MVA methods for analysis (comma separate multiple methods):
5. Select the cut to be performed on the MVA variable:
#---------------------------------------------------------------------------------

#--- Example 2 -------------------------------------------------------------------
Example: ./analysis_tool -m ADST_inputfile1.root ADST_inputfile2.root < infile
Input file (infile):
0
xmax,shfoot,shwsize,risetime
1
1
MLPBNN
0.52

The above values are prompted by:
1. Write out observables to tmva.root (0) or just create a MVA analysis on the existing tmva.root (1)?
2. Select the observables to include in the MVA (comma separated):
3. Select one to analyze (1-2):
4. Take all events as background (0) or only inverse events (1)?
5. Select one of the above MVA methods for analysis (comma separate multiple methods):
6. Select the cut to be performed on the MVA variable:
#---------------------------------------------------------------------------------

#--- Example 3 -------------------------------------------------------------------
Example: ./analysis_tool -f ADST_inputfile.root < infile
Input file (infile):
0
risetime:xmaxqual
0
1
1500

The above values are prompted by:
1. The opened file has 100 events. Read all of them (0) or select a specific event (1)?
2. Please enter the plot directive:
3. Enter maximal Xmax error to be used (set to 0, to select all showers):
4. Use vertical (0) or slant depth (1):
5. Enter maximal distance of SD tank from shower core (set to 0, to select all tanks):
#---------------------------------------------------------------------------------
