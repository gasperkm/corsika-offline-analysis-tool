#ifndef _ADSTMVA_H_
#define _ADSTMVA_H_

#include <string>
#include <vector>
#include <fstream>

#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMarker.h"
#include "TBox.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TColor.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "OfflineInclude.h"

#define SUBPLOT true
#define ALLEYES 5

// Class holding the root file structure
class Observables {
public:
   Observables();
   virtual ~Observables();

   double xmax;		// depth of shower maximum
   double x0;		// GH shape parameter X0
   double x1;		// GH shape parameter X1 (first interaction)
   double lambda;	// GH shape parameter lambda
   double shfoot;	// Shower foot
   double fdenergy;	// Total energy from FD
};

class AdstMva {
public:
   // Constructor and destructor
   AdstMva();
   virtual ~AdstMva();

   // Rewrite
   void RewriteObservables(int innr, Observables sig, Observables back);

   // Reconstructed event variables (from the ADST format)
   RecEvent *fRecEvent;
   RecEventFile *fFile;
   DetectorGeometry *fDetGeo;

   // All input filenames that will be used for the MVA
   std::vector<std::string> inname;
   // Output file created for the MVA
   TFile *outfile;
   TTree *sig_tree;
   TTree *back_tree;
   std::string outname;

   // Active eyes
   std::vector<FdRecShower> acteyes;
   int GetEyeLongestTrack();

   int GetShowerFoot(int longestEye, std::vector<FDEvent> fdevt);
   double shfootlimit;
   double shfoot;
};

#endif
