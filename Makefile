ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

SRCDIR=./src
IDIR=./include
MVADIR=./root_mva

INC=-I$(ROOTINC) -I$(IDIR) -I$(AUGEROFFLINEROOT)/include/adst
LIB=-L. -L$(AUGEROFFLINEROOT)/lib
LIBSO=-lRecEventKG -lAnalysisKG -lAugerEventIO -lAugerFramework -lAugerModules -lAugerTools -lAugerUtilities -lTMVA

#CFILES=$(shell find -maxdepth 1 -name "*.cc")
CFILES=$(SRCDIR)/analysis_tool.cc $(SRCDIR)/massanalyse.cc $(SRCDIR)/adstanalyse.cc $(SRCDIR)/adst_mva.cc $(SRCDIR)/combine_connection.cc

all: analysis_tool main_usage histogram_replot replot_usage rootadd rootadd_usage tmvacombine tmvacombine_usage tmvagui tmvagui_usage
debug: analysis_tool_dbg main_usage
usage: main_usage replot_usage rootadd_usage tmvacombine_usage tmvagui_usage

analysis_tool: $(CFILES) $(IDIR)/analysis_tool.h
	@echo "\nCompiling the main program."
	$(CXX) $(INC) $(LIB) -std=c++11 $(CFILES) -o analysis_tool $(ROOTLIB) $(LIBSO)

analysis_tool_debug: $(CFILES) $(IDIR)/analysis_tool.h
	@echo "\nCompiling the program."
	$(CXX) $(INC) $(LIB) -std=c++11 $(CFILES) -g -o analysis_tool $(ROOTLIB) $(LIBSO)

main_usage:
	@echo "\nUsage: ./analysis_tool [OPTIONS]"
	@echo "[OPTIONS]:"
	@echo "   -f [input ROOT analysis file]"
	@echo "   -t [.tar.gz of multiple ROOT analysis files] (when using this option, the files inside the tar-ball"
	@echo "      should not be in any folders)."
	@echo "   -m [all ADST ROOT files used in multivariate analysis]"
	@echo "   -mg [all ADST ROOT files used in multivariate analysis] (same as above, but with a GUI after traning MVA).\n"

histogram_replot: $(SRCDIR)/histogram_replot.cc
	@echo "\nCompiling Histogram replot program."
	$(CXX) $(INC) $(LIB) -std=c++11 $< -o histogram_replot $(ROOTLIB) $(LIBSO)

replot_usage:	
	@echo "\nUsage: ./histogram_replot [OPTIONS] [INPUT FILE NAMES]"
	@echo "[OPTIONS]:"
	@echo "   -f filenames following this option will all be replotted in the same way"
	@echo "   -s in addition to the replot, make also statistical analysis\n"

rootadd: $(SRCDIR)/hadd.C
	@echo "\nCompiling program for combining separate ADST root files into one."
	$(CXX) $(INC) $(LIB) $< -o hadd $(ROOTLIB) $(LIBSO)

rootadd_usage:	
	@echo "\nUsage: ./hadd [OUTPUT FILE NAME] [INPUT FILE NAMES]\n"

tmvacombine: $(SRCDIR)/combine.cc
	@echo "\nCompiling program for combining ADST root files from different versions of Offline into one."
	$(CXX) $(INC) $(LIB) $< $(SRCDIR)/adst_mva.cc -o combine $(ROOTLIB) $(LIBSO)

tmvacombine_usage:	
	@echo "\nUsage: ./combine [OUTPUT FILE NAME] [INPUT FILE NAMES]\n"

tmvagui: $(MVADIR)/TMVAGui.C
	@echo "\nCompiling program for viewing MVA results."
	$(CXX) -I$(ROOTINC) -I. $< -o tmvagui $(ROOTLIB)

tmvagui_usage:	
	@echo "\nUsage: ./tmvagui [MVA ANALYSIS INPUT FILE NAME]\n"

clean:
	rm -f analysis_tool $(IDIR)/workstation.h $(IDIR)/OfflineInclude.h histogram_replot hadd tmvagui combine $(MVADIR)/*.so $(MVADIR)/*.d
