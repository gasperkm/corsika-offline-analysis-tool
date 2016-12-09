#include "analysis_tool.h"
#include "massanalyse.h"
#include "adstanalyse.h"
#include "adst_mva.h"
#include "workstation.h"
#include "colors.h"
#include "combine.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <vector>

using namespace std;

/*// Int to string conversion
string IntToStr(int nr)
{
   stringstream ss;
   ss << nr;
   return ss.str();
}

// Double to string conversion
string DblToStr(double nr)
{
   stringstream ss;
   ss << nr;
   return ss.str();
}*/

// Check if file is in format massanalysis (1) or ADST (2)
int CheckFormat(char *infile)
{
//   char ctemp[1024];

//   sprintf(ctemp, "%s/%s", BASEDIR, infile);
   tempfile = new TFile(infile,"READ");
   if(tempfile->IsOpen())
   {
      if( (tempfile->GetListOfKeys()->Contains("showsimdata")) || (tempfile->GetListOfKeys()->Contains("showsrecdata")) || (tempfile->GetListOfKeys()->Contains("showfrecdata")) )
      {
         cout << "File is in massanalysis format." << endl;
	 tempfile->Close();
         return 1;
      }
      else if( (tempfile->GetListOfKeys()->Contains("recEvent")) || (tempfile->GetListOfKeys()->Contains("eventInfo")) )
      {
         cout << "File is in ADST format." << endl;
	 tempfile->Close();
         return 2;
      }
      else
      {
	 tempfile->Close();
         return -1;
      }
   }
   else
   {
      tempfile->Close();
      return 0;
   }
   return 2;
}

void MethodList()
{
   Color::Modifier yellow(Color::FG_YELLOW);

   cout << endl << "Possible MVA methods to use:" << endl;
   cout << "- Cut optimisation:                                     Cuts, CutsD, CutsPCA, CutsGA, CutsSA" << endl
        << "- 1-dimensional likelihood:                             Likelihood, LikelihoodD, LikelihoodPCA, LikelihoodKDE, LikelihoodMIX" << endl
        << "- Multidimensional likelihood and nearest neighbours:   PDERS, PDERSD, PDERSPCA, PDEFoam, PDEFoamBoost, KNN" << endl
        << "- Linear discriminant:                                  LD, Fisher, FisherG, BoostedFisher, HMatrix" << endl
        << "- Functional discriminant:                              FDA_GA, FDA_SA, FDA_MC, FDA_MT, FDA_GAMT, FDA_MCMT" << endl
        << "- Neural networks:                                      MLP, MLPBFGS, MLPBNN, CFMlpANN, TMlpANN" << endl
        << "- Support vector machine:                               SVM" << endl
        << "- Boosted decision trees:                               BDT, BDTG, BDTB, BDTD, BDTF" << endl
        << "- Friedman's rulefit:                                   RuleFit" << endl
        << "- Default collection of methods:                        Default" << endl;

   cout << yellow << "Select one of the above MVA methods for analysis (comma separate multiple methods): ";
}

vector<string> AddVariables(TMVA::Factory *factory, vector<string> obs)
{
   if(obs[0] == "default")
   {
      vector<string> tempobs;

      if(obs.size() > 1)
      {
         for(int i = 1; i < obs.size(); i++)
         {
            if((obs[i] != "xmax") && (obs[i] != "shfoot") && (obs[i] != "shwsize") && (obs[i] != "nrmu") && (obs[i] != "risetime"))
               tempobs.push_back(obs[i]);
         }
      }

      obs.erase(obs.begin(), obs.end());
      obs.push_back("xmax");
      obs.push_back("shfoot");
      obs.push_back("shwsize");
      obs.push_back("nrmu");
      obs.push_back("risetime");

      if(tempobs.size() > 0)
      {
         for(int i = 0; i < tempobs.size(); i++)
	    obs.push_back(tempobs[i]);
      }
   }

   for(int i = 0; i < obs.size(); i++)
   {
      if(obs[i] == "xmax")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "x0")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "lambda")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "fdenergy")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "shfoot")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "shwsize")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "nrmu")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "curvature")
         factory->AddVariable(obs[i], 'F');

      if(obs[i] == "risetime")
         factory->AddVariable(obs[i], 'F');
   }

   return obs;
}

vector<string> BookTheMethod(TMVA::Factory *factory, vector<string> methods)
{
   if(methods[0] == "Default")
   {
      vector<string> tempmethods;

      if(methods.size() > 1)
      {
         for(int i = 1; i < methods.size(); i++)
         {
            if((methods[i] != "Cuts") && (methods[i] != "CutsD") && (methods[i] != "Likelihood") && (methods[i] != "LikelihoodPCA") && (methods[i] != "PDERS") && (methods[i] != "KNN") && (methods[i] != "LD") && (methods[i] != "FDA_GA") && (methods[i] != "MLPBNN") && (methods[i] != "SVN") && (methods[i] != "BDT") && (methods[i] != "RuleFit"))
               tempmethods.push_back(methods[i]);
         }
      }

      methods.erase(methods.begin(), methods.end());
      methods.push_back("Cuts");
      methods.push_back("CutsD");
      methods.push_back("Likelihood");
      methods.push_back("LikelihoodPCA");
      methods.push_back("PDERS");
      methods.push_back("KNN");
      methods.push_back("LD");
      methods.push_back("FDA_GA");
      methods.push_back("MLPBNN");
      methods.push_back("SVN");
      methods.push_back("BDT");
      methods.push_back("RuleFit");

      if(tempmethods.size() > 0)
      {
         for(int i = 0; i < tempmethods.size(); i++)
	    methods.push_back(tempmethods[i]);
      }
   }

   for(int i = 0; i < methods.size(); i++)
   {
      if(methods[i] == "Cuts")
         factory->BookMethod(TMVA::Types::kCuts, methods[i], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart");

      if(methods[i] == "CutsD")
         factory->BookMethod(TMVA::Types::kCuts, methods[i], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

      if(methods[i] == "CutsPCA")
         factory->BookMethod(TMVA::Types::kCuts, methods[i], "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

      if(methods[i] == "CutsGA")
         factory->BookMethod(TMVA::Types::kCuts, methods[i], "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

      if(methods[i] == "CutsSA")
         factory->BookMethod(TMVA::Types::kCuts, methods[i], "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

      if(methods[i] == "Likelihood")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[i], "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50");

      if(methods[i] == "LikelihoodD")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[i], "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate");

      if(methods[i] == "LikelihoodPCA")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[i], "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA"); 

      if(methods[i] == "LikelihoodKDE")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[i], "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50"); 

      if(methods[i] == "LikelihoodMIX")
         factory->BookMethod(TMVA::Types::kLikelihood, methods[i], "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50"); 

      if(methods[i] == "PDERS")
         factory->BookMethod(TMVA::Types::kPDERS, methods[i], "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600");

      if(methods[i] == "PDERSD")
         factory->BookMethod(TMVA::Types::kPDERS, methods[i], "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate");

      if(methods[i] == "PDERSPCA")
         factory->BookMethod(TMVA::Types::kPDERS, methods[i], "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA");

      if(methods[i] == "PDEFoam")
         factory->BookMethod(TMVA::Types::kPDEFoam, methods[i], "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");

      if(methods[i] == "PDEFoamBoost")
         factory->BookMethod(TMVA::Types::kPDEFoam, methods[i], "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T");

      if(methods[i] == "KNN")
         factory->BookMethod(TMVA::Types::kKNN, methods[i], "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim");

      if(methods[i] == "LD")
         factory->BookMethod(TMVA::Types::kLD, methods[i], "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

      if(methods[i] == "Fisher")
         factory->BookMethod(TMVA::Types::kFisher, methods[i], "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");

      if(methods[i] == "FisherG")
         factory->BookMethod(TMVA::Types::kFisher, methods[i], "H:!V:VarTransform=Gauss");

      if(methods[i] == "BoostedFisher")
         factory->BookMethod(TMVA::Types::kFisher, methods[i], "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring");

      if(methods[i] == "HMatrix")
         factory->BookMethod(TMVA::Types::kHMatrix, methods[i], "!H:!V:VarTransform=None");

      if(methods[i] == "FDA_GA")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1");

      if(methods[i] == "FDA_SA")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

      if(methods[i] == "FDA_MC")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1");

      if(methods[i] == "FDA_MT")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch");

      if(methods[i] == "FDA_GAMT")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim");

      if(methods[i] == "FDA_MCMT")
         factory->BookMethod(TMVA::Types::kFDA, methods[i], "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20");

      if(methods[i] == "MLP")
         factory->BookMethod(TMVA::Types::kMLP, methods[i], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");

      if(methods[i] == "MLPBFGS")
         factory->BookMethod(TMVA::Types::kMLP, methods[i], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator");

      if(methods[i] == "MLPBNN")
         factory->BookMethod(TMVA::Types::kMLP, methods[i], "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator");

      if(methods[i] == "CFMlpANN")
         factory->BookMethod(TMVA::Types::kCFMlpANN, methods[i], "!H:!V:NCycles=2000:HiddenLayers=N+1,N");

      if(methods[i] == "TMlpANN")
         factory->BookMethod(TMVA::Types::kTMlpANN, methods[i], "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3");

      if(methods[i] == "SVM")
         factory->BookMethod(TMVA::Types::kSVM, methods[i], "Gamma=0.25:Tol=0.001:VarTransform=Norm");

      if(methods[i] == "BDT")
         factory->BookMethod(TMVA::Types::kBDT, methods[i], "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

      if(methods[i] == "BDTG")
         factory->BookMethod(TMVA::Types::kBDT, methods[i], "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

      if(methods[i] == "BDTB")
         factory->BookMethod(TMVA::Types::kBDT, methods[i], "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

      if(methods[i] == "BDTD")
         factory->BookMethod(TMVA::Types::kBDT, methods[i], "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

      if(methods[i] == "BDTF")
         factory->BookMethod(TMVA::Types::kBDT, "BDTMitFisher", "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

      if(methods[i] == "RuleFit")
         factory->BookMethod(TMVA::Types::kRuleFit, methods[i], "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02");

   }

   return methods;
}

struct tokens: std::ctype<char> 
{
   tokens(): std::ctype<char>(get_table()) {}

   static std::ctype_base::mask const* get_table()
   {
      typedef std::ctype<char> cctype;
      static const cctype::mask *const_rc= cctype::classic_table();

      static cctype::mask rc[cctype::table_size];
      std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));

      rc[','] = std::ctype_base::space; 
      rc[' '] = std::ctype_base::space; 
      return &rc[0];
   }
};

// Main function
int main(int argc, char **argv)
{
   MassAnalyseTool *mantool = new MassAnalyseTool();
   AdstAnalyseTool *aantool = new AdstAnalyseTool();
   AdstMva *mvatool = new AdstMva();

   // if the argument is -f, only one file is supplied (argtype = 0): But this file can have one or more events
   // if the argument is -t, the files are saved into a tar-ball (argtype = 1)
   // if the argument is -m or -mg, the files are saved into the ADST format and will need MVA (argtype = 2)
   argtype = 0;

   int itemp;
   string stemp, stemp2;

//   printf("1\n");

   if(argc > 1)
   {
      for(int i = 1; i < argc; i++)
      {
         // Analysis of a single file
         if(strcmp("-f",argv[i]) == 0)
	 {
	    argtype = 0;

            if(i <= 2)
               fileformat = CheckFormat(argv[i+1]);
            if(fileformat == 1)
               mantool->SetAnalysisFilename(argv[i+1]);
	    else if(fileformat == 2)
               aantool->SetAnalysisFilename(argv[i+1]);
	    else
	    {
	       cout << "The selected files are not supported for this analysis. Please use only massanalysis or ADST files." << endl;
	       return -1;
	    }
	       
//	    cout << "Analysis file = " << mantool->GetAnalysisFilename() << endl;
	    i++;
	 }
	 // Analysis of a collection of files in a tar-ball
	 else if(strcmp("-t",argv[i]) == 0)
	 {
	    argtype = 1;
	    cout << "Using a tar-ball of shower reconstructions. Untarring them to directory " << BASEDIR << "/input." << endl;
	    // Deleting previous input files from the input directory
	    stemp = "rm -r " + string(BASEDIR) + "/input/*";
	    itemp = system(stemp.c_str());
	    // Saving the ROOT analysis file names from the tar-ball
	    stemp = "tar -ztf " + string(argv[i+1]) + " > " + string(BASEDIR) + "/input/tarball_names.txt";
	    itemp = system(stemp.c_str());
	    // Untarring the ROOT analysis files
	    stemp = "tar -zxf " + string(argv[i+1]) + " -C " + string(BASEDIR) + "/input";
	    itemp = system(stemp.c_str());
	    i++;
	 }
	 // Multivariate analysis with TMVA
	 else if( (strcmp("-m",argv[i]) == 0) || (strcmp("-mg",argv[i]) == 0) || (strcmp("-mgc",argv[i]) == 0) || (strcmp("-mc",argv[i]) == 0) )
	 {
	    if( (strcmp("-mg",argv[i]) == 0) || (strcmp("-mgc",argv[i]) == 0) )
	       mvatool->graphical = true;

	    if( (strcmp("-mgc",argv[i]) == 0) || (strcmp("-mc",argv[i]) == 0) )
	       argtype = 3;
	    else
	       argtype = 2;
	    cout << "Rewriting files and preparing for the multivariate analysis." << endl;

	    cout << "Nr. of input files: " << argc-2 << endl;

	    for(int i = 0; i < argc-2; i++)
	       mvatool->inname.push_back(string(argv[i+2]));
	 }
      }
   }
   else
   {
      cout << "No arguments were supplied." << endl;
      return 1;
   }

   if(fileformat == 2)
      aantool->argtype = argtype;

   // Situation, where we only have one input ROOT analysis file (that could have multiple events)
   if(argtype == 0)
   {
      // Format is massanalysis
      if(fileformat == 1)
      {
         delete aantool;
	 delete mvatool;

         // Open the analysis file and retrieve all default trees (those not tankvem and eyelong)
         mantool->infile = new TFile((mantool->GetAnalysisFilename()).c_str(),"READ");
         mantool->GetDefaultTrees();

         mantool->filesel = 0;
         mantool->GetActiveEyes();
         mantool->Hello();
         mantool->RunDirective(0);
         cout << "################ Finish RunDirective (" << (mantool->GetAnalysisFilename()) << ")" << endl;
         
         mantool->infile->Close();
         delete mantool;
      }
      // Format is ADST
      else if(fileformat == 2)
      {
         delete mantool;
	 delete mvatool;

         // Open the analysis file
	 aantool->fFile = new RecEventFile((aantool->GetAnalysisFilename()).c_str(), RecEventFile::eRead);
	 aantool->fFile->SetBuffers(&(aantool->fRecEvent));
	 aantool->fFile->ReadDetectorGeometry(*(aantool->fDetGeo));

	 aantool->ReadOption();  // Checks how many events are stored in the ADST file

	 if(aantool->readopt == 1)
	 {
            aantool->filesel = 0;
	    aantool->GetActiveEyes();
            aantool->Hello();
            aantool->RunDirective(aantool->readevent);
	 }
	 else if(aantool->readopt == 0)
         {
            aantool->filesel = 1;
            aantool->Hello();
            // Define all possible plotting types (will then delete the ones that are not used)
            aantool->gr = new TGraph(aantool->fFile->GetNEvents());
            aantool->grErr = new TGraphErrors(aantool->fFile->GetNEvents());
            aantool->grAsymmErr = new TGraphAsymmErrors(aantool->fFile->GetNEvents());
            aantool->histf = new TH1F("th1","",50,1,1);
            aantool->histf->SetBit(TH1::kCanRebin);
            aantool->hist2f = new TH2F("th2","",sqrt(aantool->fFile->GetNEvents()),0,4000,sqrt(aantool->fFile->GetNEvents()),0,4000);
            // Initiate the ranges for plots
            aantool->InitRanges();
	    aantool->evtcount = aantool->fFile->GetNEvents();

	    for(int i = 0; i < aantool->evtcount; i++)
	    {
	       cout << "Event number = " << i << endl;
               aantool->readevent = i;
               aantool->GetActiveTanks();
               aantool->GetActiveEyes();
               aantool->RunDirective(aantool->readevent);
	    }
/*            for(int i = 0; i < aantool->tarnames.size(); i++)
            {
               aantool->SetAnalysisFilename((char*)aantool->tarnames[i].c_str());

	       aantool->fFile = new RecEventFile((aantool->GetAnalysisFilename()).c_str(), RecEventFile::eRead);
	       aantool->fFile->SetBuffers(&(aantool->fRecEvent));
	       aantool->fFile->ReadDetectorGeometry(*(aantool->fDetGeo));
               aantool->GetActiveTanks();
               aantool->GetActiveEyes();
               
               aantool->RunDirective(i);
               cout << "################ Finish RunDirective (" << (aantool->GetAnalysisFilename()) << ")" << endl;

               aantool->fFile->Close();
            }*/
	 }
         cout << "################ Finish RunDirective (" << (aantool->GetAnalysisFilename()) << ")" << endl;
	 
	 aantool->fFile->Close();
        // delete aantool;
      }
   }
   // Situation, where we have multiple input ROOT analysis files saved in a tar-ball
   else if(argtype == 1)
   {
      itemp = 0;

      // Open the tarball_names.txt file and save the filenames
      ifstream infile;
      stemp = string(BASEDIR) + "/input/tarball_names.txt";
      infile.open(stemp.c_str(), ifstream::in);
      if(infile.is_open())
      {
         while(infile.good())
         {
	    if(infile.eof()) break;

            infile >> stemp;
            stemp = string(BASEDIR) + "/input/" + stemp;

            if(itemp == 0)
               fileformat = CheckFormat((char*)stemp.c_str());

            if(fileformat == 1)
               mantool->tarnames.push_back(stemp);
	    else if(fileformat == 2)
               aantool->tarnames.push_back(stemp);
	    else
	    {
	       cout << "The selected files are not supported for this analysis. Please use only massanalysis or ADST files." << endl;
	       return -1;
	    }
            
   	    infile.ignore(1,' ');
         }
      }
      infile.close();

      // Format is massanalysis
      if(fileformat == 1)
      {
         delete aantool;
	 delete mvatool;

         // Start the main loop execution over all files that were in the tar-ball
         mantool->filesel = 1;
         mantool->Hello();
         // Define all possible plotting types (will then delete the ones that are not used)
         mantool->gr = new TGraph(mantool->tarnames.size());
         mantool->grErr = new TGraphErrors(mantool->tarnames.size());
         mantool->grAsymmErr = new TGraphAsymmErrors(mantool->tarnames.size());
         mantool->histf = new TH1F("th1","",50,1,1);
         mantool->histf->SetBit(TH1::kCanRebin);
         mantool->hist2f = new TH2F("th2","",sqrt(mantool->tarnames.size()),0,4000,sqrt(mantool->tarnames.size()),0,4000);
         // Initiate the ranges for plots
         mantool->InitRanges();

         for(int i = 0; i < mantool->tarnames.size(); i++)
         {
            mantool->SetAnalysisFilename((char*)mantool->tarnames[i].c_str());

            mantool->infile = new TFile((mantool->GetAnalysisFilename()).c_str(),"READ");
            mantool->GetDefaultTrees();
            mantool->GetActiveTanks();
            mantool->GetActiveEyes();
            
            mantool->RunDirective(i);
            cout << "################ Finish RunDirective (" << (mantool->GetAnalysisFilename()) << ")" << endl;

            mantool->infile->Close();
         }
         delete mantool;
      }
      // Format is ADST
      else if(fileformat == 2)
      {
         delete mantool;
	 delete mvatool;

         // Start the main loop execution over all files that were in the tar-ball
         aantool->filesel = 1;
         aantool->Hello();
         // Define all possible plotting types (will then delete the ones that are not used)
         aantool->gr = new TGraph(aantool->tarnames.size());
         aantool->grErr = new TGraphErrors(aantool->tarnames.size());
         aantool->grAsymmErr = new TGraphAsymmErrors(aantool->tarnames.size());
         aantool->histf = new TH1F("th1","",50,1,1);
         aantool->histf->SetBit(TH1::kCanRebin);
         aantool->hist2f = new TH2F("th2","",sqrt(aantool->tarnames.size()),0,4000,sqrt(aantool->tarnames.size()),0,4000);
         // Initiate the ranges for plots
         aantool->InitRanges();
	 aantool->evtcount = aantool->tarnames.size();

         for(int i = 0; i < aantool->evtcount; i++)
         {
            aantool->SetAnalysisFilename((char*)aantool->tarnames[i].c_str());

	    aantool->fFile = new RecEventFile((aantool->GetAnalysisFilename()).c_str(), RecEventFile::eRead);
	    aantool->fFile->SetBuffers(&(aantool->fRecEvent));
	    aantool->fFile->ReadDetectorGeometry(*(aantool->fDetGeo));
            aantool->GetActiveTanks();
            aantool->GetActiveEyes();
            
            aantool->RunDirective(i);
            cout << "################ Finish RunDirective (" << (aantool->GetAnalysisFilename()) << ")" << endl;

            aantool->fFile->Close();
         }
         delete aantool;
      }
   }
   else if( (argtype == 2) || (argtype == 3) )
   {
      delete mantool;
      delete aantool;

      // Prepare colors for colored terminal output
      Color::Modifier red(Color::FG_RED);
      Color::Modifier cyan(Color::FG_CYAN);
      Color::Modifier yellow(Color::FG_YELLOW);
      Color::Modifier def(Color::FG_DEFAULT);

      // Prepare output file
      stemp = "analysis_out/tmva.root";

      // Determine if a writeout to tmva.root is needed, or if only a TMVA analysis should be performed
      int writeAnalysis = -1;
      if( (mvatool->inname.size() == 0) || ((argtype == 3) && (mvatool->inname.size() == 1)) )
      {
         if(mvatool->inname.size() == 0)
            cout << cyan << "No input files supplied. Only running a MVA analysis on the existing tmva_rewrite.root." << def << endl;
	 else
	    cout << cyan << "No input files supplied. Only running a MVA analysis on the existing " << mvatool->inname[0] << def << endl;
	 writeAnalysis = 1;
      }
      else
      {
         cout << yellow << "Write out observables to tmva_rewrite.root and perform analysis (0), only create a MVA analysis on the existing tmva_rewrite.root (1) or only write out observables (2)? " << def;
         cin >> writeAnalysis;
         while( (writeAnalysis != 0) && (writeAnalysis != 1) && (writeAnalysis != 2) )
         {
            cout << red << "Error: Please select a valid option. " << yellow << "Write out observables to tmva_rewrite.root and perform analysis (0), only create a MVA analysis on the existing tmva_rewrite.root (1) or only write out observables (2)? " << def;
            cin >> writeAnalysis;
         }
      }

      // TODO: Add writeouts for (mean - error) and (mean + error)

      mvatool->outname = string(BASEDIR) + "/" + stemp;
      if( (writeAnalysis == 0) || (writeAnalysis == 2) )
      {
         Observables obssig, obsback;
         TTree *back_tree;

         for(int j = 0; j < 3; j++)
	 {
	    if(j == 0)
               mvatool->outfile = TFile::Open((mvatool->outname).c_str(),"RECREATE");
	    else if(j == 1)
	    {
               mvatool->outname = string(BASEDIR) + "/analysis_out/tmva_negerror.root";
               mvatool->outfile = TFile::Open((mvatool->outname).c_str(),"RECREATE");
	    }
	    else if(j == 2)
	    {
               mvatool->outname = string(BASEDIR) + "/analysis_out/tmva_poserror.root";
               mvatool->outfile = TFile::Open((mvatool->outname).c_str(),"RECREATE");
	    }

            // Prepare observable values for signal and background
            back_tree = new TTree[mvatool->inname.size()];
#ifdef OFFLINEOLD
            for(int i = 0; i < mvatool->inname.size(); i++)
               back_tree[i].SetNameTitle(("TreeOldB" + IntToStr(i+1)).c_str(), ("Background tree without events from file " + mvatool->inname[i]).c_str());
            mvatool->all_tree = new TTree("TreeA", "Background tree with all events, including signal events.");
#elif defined OFFLINENEW
            mvatool->PrepareOtherTrees(mvatool->inname.size(), 1);
            for(int i = 0; i < mvatool->inname.size(); i++)
               back_tree[i].SetNameTitle(("TreeNewB" + IntToStr(i+1)).c_str(), ("Background tree without events from file " + mvatool->inname[i]).c_str());
            mvatool->all_tree = new TTree("TreeA", "Background tree with all events, including signal events.");
#endif
      
            // Start rewriting (for all input files)
            for(int i = 0; i < mvatool->inname.size(); i++)
               mvatool->RewriteObservables(i, obssig, obsback, back_tree, j);

#ifdef OFFLINEOLD
            mvatool->PrepareOtherTrees(mvatool->inname.size(), 1);
            mvatool->all_tree->Write();
            for(int i = 0; i < mvatool->inname.size(); i++)
               back_tree[i].Write();
            mvatool->PrepareOtherTrees(mvatool->inname.size(), 0);
#elif defined OFFLINENEW
            mvatool->all_tree->Write();
            mvatool->PrepareOtherTrees(mvatool->inname.size(), 0);
            for(int i = 0; i < mvatool->inname.size(); i++)
               back_tree[i].Write();
#endif

            // Close all open files
            mvatool->outfile->Close();
            mvatool->fFile->Close();

	    delete[] back_tree;
	 }
  
         // again setup the outname, because we have written out the mean and error values to separate files
         mvatool->outname = string(BASEDIR) + "/" + stemp;

         // If performing the MVA analysis straight away, rewrite the root file structure
         if(writeAnalysis == 0)
	 {
            stemp = string(BASEDIR) + "/analysis_out/tmva_rewrite.root";
            char *tmpfiles[3];
            for(int i = 0; i < 3; i++) tmpfiles[i] = new char[1024];
            sprintf(tmpfiles[0], "./combine");
            sprintf(tmpfiles[1], "%s", stemp.c_str());
            sprintf(tmpfiles[2], "%s", (mvatool->outname).c_str());
            char **filesptr = tmpfiles;
            hadd(3, filesptr);
            mvatool->outname = stemp;
	 }
	 else if(writeAnalysis == 2)
	 {
            delete mvatool;
	    return 0;
         }
      }
      else
      {
         if( (argtype == 3) && (mvatool->inname.size() == 1) )
            mvatool->outname = string(BASEDIR) + "/" + mvatool->inname[0];
	 else
            mvatool->outname = string(BASEDIR) + "/analysis_out/tmva_rewrite.root";

      }

      // Start performing the MVA analysis
      stemp = string(BASEDIR) + "/analysis_out/tmva_output.root";
      TFile *ofile = TFile::Open(stemp.c_str(),"RECREATE");
      // Factory usage:
      // - user-defined job name, reappearing in names of weight files for training results ("TMVAClassification")
      // - pointer to an output file (ofile)
      // - options
      // Factory has the following options:
      // 	V = verbose
      // 	Silent = batch mode
      // 	Color = colored screen output
      // 	DrawProgressBar = progress bar display during training and testing
      // 	Transformations = the transformations to make (identity, decorrelation, PCA, uniform, gaussian, gaussian decorrelation)
      // 	AnalysisType = setting the analysis type (Classification, Regression, Multiclass, Auto)
      // Default values = !V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Auto
      TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",ofile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

      // Add variables to be used (similar for spectators and targets):
      // - name of the variable as defined in the input file
      // - type of the variable (can be int = 'I' or float/double = 'F')
      // Additionally, can also have a title and the units for the variable: factory->AddVariable("name","title","unit",'F');
      string signalName = "default";
      cout << endl << cyan << "Observables to add in the multivariate analysis:" << endl
           << "   xmax, x0, lambda, fdenergy, shfoot, shwsize, nrmu, curvature, risetime" << endl
           << "   default (xmax,shfoot,shwsize,nrmu,risetime)" << endl << def;
      cout << endl << yellow << "Select the observables to include in the MVA (comma separated): " << def;
      cin >> signalName;
      stringstream ss1(signalName);
      ss1.imbue(locale(locale(), new tokens()));
      istream_iterator<string> begin1(ss1);
      istream_iterator<string> end1;
      vector<string> observables(begin1, end1);

      observables = AddVariables(factory, observables);

      // Open up the input file and ask which signal tree we want to check
      int whichAnalysis = -1;
      
      signalName = "TreeS1";

      TFile *ifile = TFile::Open((mvatool->outname).c_str());
      if((ifile->GetListOfKeys()->Contains("TreeA")) || (ifile->GetListOfKeys()->Contains("TreeB1")))
      {
         cout << endl << cyan << "There are " << ((ifile->GetNkeys()-1)/2) << " signal trees available:" << endl;
         for(int i = 1; i <= ((ifile->GetNkeys())/2); i++)
         {
            signalName = "TreeS" + IntToStr(i);
            cout << "- " << i << ": " << ifile->GetKey(signalName.c_str())->GetTitle() << endl;
         }
         cout << endl << yellow << "Select one to analyze (1-" << ((ifile->GetNkeys()-1)/2) << "): " << def;
         cin >> whichAnalysis;

         while( (whichAnalysis > ((ifile->GetNkeys())/2)) || (whichAnalysis < 1) )
         {
            cout << red << "Error: Wrong selection. " << yellow << "There are " << ((ifile->GetNkeys()-1)/2) << " signal trees available. Select one to analyze (1-" << ((ifile->GetNkeys()-1)/2) << "): " << def;
            cin >> whichAnalysis;
         }

         signalName = "TreeS" + IntToStr(whichAnalysis);
      }

      vector<int> backtrees;
      string backgroundName = "TreeS3";
      int nrbacktrees = (ifile->GetNkeys()-1)/2;

      if(nrbacktrees > 2)
      {
         cout << endl << cyan << "There are " << nrbacktrees-1 << " other signal trees available (" << signalName << " already taken as signal). How many trees should be taken as background (use 0, to skip and use simple background estimation instead)? " << def;
	 cin >> nrbacktrees;

	 if(nrbacktrees > 0)
	 {
	    cout << endl << cyan << "The available trees are:" << endl;
	    for(int i = 1; i <= (ifile->GetNkeys()-1)/2; i++)
	    {
	       if(i != whichAnalysis)
	       {
                  backgroundName = "TreeS" + IntToStr(i);
                  cout << "- " << i << ": " << ifile->GetKey(backgroundName.c_str())->GetTitle() << endl;
	       }
	    }
	    for(int i = 1; i <= nrbacktrees; i++)
	    {
               cout << endl << yellow << "Select one of these to use as background: " << def;
               cin >> itemp;
	       backtrees.push_back(itemp);
	       cout << "Currently selected trees for background: ";
	       for(int j = 0; j < i; j++) cout << "TreeS" << backtrees[j] << ", ";
	    }
	    cout << endl;
	 }
	 else
	    cout << "Swithcing to normal background estimation." << endl;
      }
      else
         nrbacktrees = 0;

      if(nrbacktrees < 1) nrbacktrees = 0;

      // Getting signal and background trees for training and testing (can supply multiple trees)
      TTree *signal = (TTree*)ifile->Get(signalName.c_str());
      cout << cyan << signalName << " selected as signal tree." << endl << def;
      TTree *background[nrbacktrees];

      if(nrbacktrees < 1)
      {
         cout << endl << yellow << "Take all events as background (0) or only inverse events (1)? " << def;
         cin >> itemp;
         while( (itemp != 0) && (itemp != 1) )
         {
            cout << red << "Error: Wrong selection. " << yellow << "Take all events as background (0) or only inverse events (1)? " << def;
            cin >> itemp;
         }

         if(itemp == 0)
         {
            background[0] = (TTree*)ifile->Get("TreeA");
            cout << cyan << "TreeA selected as background tree." << endl << def;
         }
         else if(itemp == 1)
         {
            signalName = "TreeB" + IntToStr(whichAnalysis);
            background[0] = (TTree*)ifile->Get(signalName.c_str());
            cout << cyan << signalName << " selected as background tree." << endl << def;
         }
      }
      else
      {
         for(int i = 0; i < nrbacktrees; i++)
	 {
	    backgroundName = "TreeS" + IntToStr(backtrees[i]);
            background[i] = (TTree*)ifile->Get(backgroundName.c_str());
            cout << cyan << backgroundName << " selected as background tree." << endl << def;
	 }
      }

      // Add all the trees to the factory with overall weights for the complete sample
      factory->AddSignalTree(signal, 1.0);
      if(nrbacktrees < 2)
         factory->AddBackgroundTree(background[0], 1.0);
      else
         for(int i = 0; i < nrbacktrees; i++)
	    factory->AddBackgroundTree(background[i], 1.0);

      // Preparing and training from the trees:
      // - preselection cuts make cuts on input variables, before even starting the MVA
      // - options
      // These are the possible options:
      // 	nTrain_Signal = number of training events of class Signal (0 takes all)
      // 	nTrain_Background = number of training events of class Background (0 takes all)
      // 	nTest_Signal = number of test events of class Signal (0 takes all)
      // 	nTest_Background = number of test events of class Background (0 takes all)
      // 	SplitMode = method of choosing training and testing events (Random, Alternate, Block)
      //	NormMode = renormalisation of event-by-event weights for training (NumEvents: average weight of 1 per event for signal and background, EqualNumEvents: average weight of 1 per event for signal and sum of weights for background equal to sum of weights for signal, None)
      //	V = verbose
      //	MixMode = method of mixing events of different classes into one dataset (SameAsSplitMode, Random, Alternate, Block)
      //	SplitSeed = seed for random event shuffling (default = 100)
      //	VerboseLevel = level of verbose (Debug, Verbose, Info)
      factory->PrepareTrainingAndTestTree("", "", "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

      // Choose the MVA methods
      cout << cyan;
      MethodList();
      cout << def;

      cin >> signalName;
      stringstream ss2(signalName);
      ss2.imbue(locale(locale(), new tokens()));
      istream_iterator<string> begin2(ss2);
      istream_iterator<string> end2;
      vector<string> methods(begin2, end2);

      // Booking MVA methods:
      // - type of MVA method to be used (all defined in src/Types.h)
      // - the unique name for the MVA method suplied by the user
      // - options
      // The possible options for each method are defined here: http://tmva.sourceforge.net/optionRef.html
      // For example:
      // 	H = print method-specific help message
      // 	V = verbose
      // 	NeuronType = neuron activation function type (default = sigmoid)
      // 	VarTransform = list of variable transformations to do before training (D_Background,P_Signal,G,N_AllClasses -> N = Normalization for all classes)
      // 	NCycles = number of training cycles
      // 	HiddenLayers = hidden layer architecture (default = N,N-1)
      // 	TestRate = test for overtraining at each #th epoch (default = 10)
      // 	TrainingMethod = train with back propagation (BP), BFGS algorithm (BFGS) or generic algorithm (GA)
      // 	UseRegulator = use regulator to avoid overtraining
      methods = BookTheMethod(factory, methods);
//    factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" );
      // Train the selected methods and save them to the weights folder
      factory->TrainAllMethods();
      // Test the selected methods by applying the trained data to the test data set -> outputs saved to TestTree output file and then to the output ROOT file
      factory->TestAllMethods();
      // Evaluation of methods printed to stdout
      factory->EvaluateAllMethods();

      ifile->Close();
      delete factory;
      ofile->Close();

      // Open the GUI to check for best MVA
      if(mvatool->graphical)
         itemp = system("./tmvagui analysis_out/tmva_output.root");
      cout << cyan << "Closing the GUI and continuing by applying the selected MVA to the data." << endl << def;
      
      // Perform classification application (only apply it, if one method has been selected!)
      if(methods.size() > 1)
      {
         cout << cyan << "More than one MVA method selected. To perform classification application, select only one." << endl << def;
         return 0;
      }

      string applymva;
      // Select the method to use
      if(methods.size() == 1)
         applymva = methods[0];
      else if(methods.size() > 1)
      {
         cout << endl << cyan << "The following methods were trained and tested:" << endl;
         for(int i = 0; i < methods.size(); i++)
            cout << "   " << methods[i] << endl;
         cout << yellow << "Select the method to be applied to the data (only one): " << def;
         cin >> applymva;
      }
     
      // Open a reader and add variables (must be the same as for the training)
      TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

      float obsvars[observables.size()];

      for(int i = 0; i < observables.size(); i++)
      {
         if(observables[i] == "xmax")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "x0")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "lambda")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "fdenergy")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "shfoot")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "shwsize")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "nrmu")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "curvature")
            reader->AddVariable(observables[i], &obsvars[i]);

         if(observables[i] == "risetime")
            reader->AddVariable(observables[i], &obsvars[i]);
      }

      // Book the MVA with the produced weights file
      signalName = string(BASEDIR) + "/weights/TMVAClassification_" + applymva + ".weights.xml";
      string mvamethod = (string)(applymva + " method");
      reader->BookMVA(mvamethod, signalName);

/*      // GKM
      FILE *fpsig = fopen((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal.txt").c_str(),"w");			// signal + wrong back after MVA cut
      FILE *fpback = fopen((string(BASEDIR) + "/root_mva/plots/gkm_simple_back.txt").c_str(),"w");			// back + wrong signal after MVA cut
      FILE *fpsigstart = fopen((string(BASEDIR) + "/root_mva/plots/gkm_simple_signal_start.txt").c_str(),"w");	// signal before MVA cut
      FILE *fpbackstart = fopen((string(BASEDIR) + "/root_mva/plots/gkm_simple_back_start.txt").c_str(),"w");		// back before MVA cut
      FILE *fpall = fopen((string(BASEDIR) + "/root_mva/plots/gkm_simple_all.txt").c_str(),"w");			// signal + back
      vector<int> backval, sigval;
      for(int i = 0; i < nrbacktrees; i++)
      {
         backval.push_back(0);
         sigval.push_back(0);
      }
      // GKM
*/
      // Open the input file and prepare the TTree
      ifile = TFile::Open((mvatool->outname).c_str());

// New
      int evalTree;
      cout << endl << cyan << "The available trees that can be evaluated by the analysis:" << endl;
      for(int i = 1; i <= (ifile->GetNkeys()-1)/2; i++)
      {
         signalName = "TreeS" + IntToStr(i);
         cout << "- " << i << ": " << ifile->GetKey(signalName.c_str())->GetTitle() << endl;
      }
      cout << endl << yellow << "Select the tree to evaluate (use -1 to evaluate all): " << def;
      cin >> evalTree;

      if(evalTree == -1)
      {
         double cutmva;
         cout << endl << yellow << "Select the cut to be performed on the MVA variable: " << def;
         cin >> cutmva;

         for(int j = 1; j <= (ifile->GetNkeys()-1)/2; j++)
	 {
            signalName = "TreeS" + IntToStr(j);
            cout << endl << cyan << signalName << " has been selected for evaluation." << endl;

            TTree *signalapp = (TTree*)ifile->Get(signalName.c_str());

            for(int i = 0; i < observables.size(); i++)
               signalapp->SetBranchAddress((observables[i]).c_str(), &obsvars[i]);

	    cout << def;
            mvatool->CreateMVAPlots(cutmva, observables, signalapp, reader, mvamethod, obsvars, signalName);
	    cout << def;
	 }
      }
      else
      {
         signalName = "TreeS" + IntToStr(evalTree);
         cout << endl << cyan << signalName << " has been selected for evaluation." << endl;

         TTree *signalapp = (TTree*)ifile->Get(signalName.c_str());

         for(int i = 0; i < observables.size(); i++)
            signalapp->SetBranchAddress((observables[i]).c_str(), &obsvars[i]);

         double cutmva;
         cout << endl << yellow << "Select the cut to be performed on the MVA variable: " << def;
         cin >> cutmva;

         mvatool->CreateMVAPlots(cutmva, observables, signalapp, reader, mvamethod, obsvars, signalName);
      }
// New

/*      for(int ievt=0; ievt < signalapp->GetEntries(); ievt++)
      {
         signalapp->GetEntry(ievt);

         for(int i = 0; i < observables.size(); i++)
	 {
//	    cout << obsvars[i] << "\t";
	    fprintf(fpall, "%lf\t", obsvars[i]);
	    fprintf(fpsigstart, "%lf\t", obsvars[i]);

            if(reader->EvaluateMVA(mvamethod) >= cutmva)
	       fprintf(fpsig, "%lf\t", obsvars[i]);
            else
	       fprintf(fpback, "%lf\t", obsvars[i]);
	 }
//	 cout << reader->EvaluateMVA(mvamethod) << "\n";
	 fprintf(fpall, "%lf\n", reader->EvaluateMVA(mvamethod));
	 fprintf(fpsigstart, "%lf\n", reader->EvaluateMVA(mvamethod));

         if(reader->EvaluateMVA(mvamethod) >= cutmva)
	    fprintf(fpsig, "%lf\n", reader->EvaluateMVA(mvamethod));
	 else
	    fprintf(fpback, "%lf\n", reader->EvaluateMVA(mvamethod));

         if(reader->EvaluateMVA(mvamethod) >= cutmva)
	    sigval[0]++;
	 else
	    backval[0]++;
      } 

      ifile->Close();

      // GKM
      cout << cyan << "Signal TTree: There were " << sigval[0] << " signal events and " << backval[0] << " background events" << endl << def;
      // GKM

      ifile = TFile::Open((mvatool->outname).c_str());
      if(nrbacktrees < 2)
      {
         signalName = "TreeB" + IntToStr(whichAnalysis);

         TTree *backgroundapp = (TTree*)ifile->Get(signalName.c_str());

         for(int i = 0; i < observables.size(); i++)
            backgroundapp->SetBranchAddress((observables[i]).c_str(), &obsvars[i]);

         backval[0] = 0; sigval[0] = 0;

         for(int ievt=0; ievt < backgroundapp->GetEntries(); ievt++)
         {
            backgroundapp->GetEntry(ievt);

            for(int i = 0; i < observables.size(); i++)
            {
               fprintf(fpall, "%lf\t", obsvars[i]);
               fprintf(fpbackstart, "%lf\t", obsvars[i]);

               if(reader->EvaluateMVA(mvamethod) >= cutmva)
                  fprintf(fpsig, "%lf\t", obsvars[i]);
               else
                  fprintf(fpback, "%lf\t", obsvars[i]);
            }
            fprintf(fpall, "%lf\n", reader->EvaluateMVA(mvamethod));
            fprintf(fpbackstart, "%lf\n", reader->EvaluateMVA(mvamethod));

            if(reader->EvaluateMVA(mvamethod) >= cutmva)
               fprintf(fpsig, "%lf\n", reader->EvaluateMVA(mvamethod));
            else
               fprintf(fpback, "%lf\n", reader->EvaluateMVA(mvamethod));

            if(reader->EvaluateMVA(mvamethod) >= cutmva)
               sigval[0]++;
            else
               backval[0]++;
         } 
      }
      else
      {
         for(int i = 0; i < nrbacktrees; i++)
	 {
            backval[i] = 0; sigval[i] = 0;
            backgroundName = "TreeS" + IntToStr(backtrees[i]);

            TTree *backgroundapp = (TTree*)ifile->Get(backgroundName.c_str());

            for(int i = 0; i < observables.size(); i++)
               backgroundapp->SetBranchAddress((observables[i]).c_str(), &obsvars[i]);

            for(int ievt=0; ievt < backgroundapp->GetEntries(); ievt++)
            {
               backgroundapp->GetEntry(ievt);

               for(int i = 0; i < observables.size(); i++)
               {
                  fprintf(fpall, "%lf\t", obsvars[i]);
                  fprintf(fpbackstart, "%lf\t", obsvars[i]);

                  if(reader->EvaluateMVA(mvamethod) >= cutmva)
                     fprintf(fpsig, "%lf\t", obsvars[i]);
                  else
                     fprintf(fpback, "%lf\t", obsvars[i]);
               }
               fprintf(fpall, "%lf\n", reader->EvaluateMVA(mvamethod));
               fprintf(fpbackstart, "%lf\n", reader->EvaluateMVA(mvamethod));

               if(reader->EvaluateMVA(mvamethod) >= cutmva)
                  fprintf(fpsig, "%lf\n", reader->EvaluateMVA(mvamethod));
               else
                  fprintf(fpback, "%lf\n", reader->EvaluateMVA(mvamethod));

               if(reader->EvaluateMVA(mvamethod) >= cutmva)
                  sigval[i]++;
               else
                  backval[i]++;
            } 
	 }
      }
*/
      ifile->Close();
/*
      // GKM
      for(int i = 0; i < nrbacktrees; i++)
         cout << cyan << "Background TTree: There were " << sigval[i] << " signal events and " << backval[i] << " background events" << endl << def;
      fclose(fpsig);
      fclose(fpback);
      fclose(fpsigstart);
      fclose(fpbackstart);
      fclose(fpall);
      // GKM
*/
      cout << def << endl;
//      mvatool->CreateMVAPlots(cutmva, observables);

      delete reader;
      delete mvatool;
   }

   return 0;
}
