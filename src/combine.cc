//macro to add root files

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include "adst_mva.h"
#include "combine.h"

/*TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames );
void RewriteRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames );
void CheckKeys( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts );
void hadd(int nrfiles, char **files);*/

// Int to string conversion
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
}

string RemoveFromFirst(string input, string toremove, string replace)
{
   int size = toremove.size();
   int pos = input.find(toremove);

   if(pos >= 0)
   {
      return input.replace(pos, size, replace);
   }
   else
      return input.replace(0, size, replace);
}

int main(int argc, char **argv)
{
   system("rm -v input/hsimple*.root");
   char confirm = 'n';

   if(argc == 3)
   {
      printf("Only one input file selected. Rewriting it into normal tree format (y/n)? ");
      scanf("%c", &confirm);

      if( (confirm == 'y') || (confirm == 'Y') )
      {
         printf("Rewriting file %s\n", argv[2]);
         hadd(argc, argv);
      }
      else
         return 1;
   }
   else if(argc < 4)
   {
      // Stop program if you forget to enter input files
      printf("No input and output data files used (1 output file and at least 2 input files needed).\n");
      return 1;
   }
   else
   {
      // Add all input files to a single output file
      hadd(argc, argv);
   }

   return 0;
}

void hadd(int nrfiles, char **files)
{
   TList *FileList;
   TFile *Target;

   char ctemp[1024];
   int total;
   vector<int> nrkeys;
   vector<int> nrevts;
   vector<string> filenames;

   // Prepare the files to be merged
   if(gSystem->AccessPathName("input/hsimple1.root"))
   {
     for(int i = 2; i < nrfiles; i++)
     {
       sprintf(ctemp, "input/hsimple%d.root", i-1);
       gSystem->CopyFile(files[i], ctemp);
       filenames.push_back(files[i]);
     }
     total = nrfiles-2;
   }

   Target = TFile::Open(files[1], "RECREATE");

   FileList = new TList();
   for(int i = 0; i < total; i++)
   {
     sprintf(ctemp, "input/hsimple%d.root", i+1);
     FileList->Add(TFile::Open(ctemp));   
   }

   CheckKeys(Target, FileList, nrkeys, nrevts);
   for(int i = 0; i < nrkeys.size(); i++)
      printf("Keys in file %d = %d\n", i, nrkeys[i]);
   for(int i = 0; i < nrevts.size(); i++)
      printf("Events in key %d = %d\n", i, nrevts[i]);

   if(total > 1)
      MergeRootfile(Target, FileList, nrkeys, nrevts, filenames);
   else
      RewriteRootfile(Target, FileList, nrkeys, nrevts, filenames);

   for(int i = 0; i < total; i++)
   {
     sprintf(ctemp, "rm -v input/hsimple%d.root", i+1);
     system(ctemp);
   }
}

void CheckKeys( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts )
{
   int nrfiles = sourcelist->GetSize();
   int fileCount = 0;
   string strOldS, strNewS, strOldB, strNewB;
   TTree *treeOldS, *treeNewS, *treeOldB, *treeNewB, *treeAll;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   while(nextsource)
   {
      nrkeys.push_back(nextsource->GetNkeys());

      treeAll = (TTree*)nextsource->Get("TreeA");
      nrevts.push_back(treeAll->GetEntries());
      for(int j = 1; j <= (nrkeys[fileCount]-1)/4; j++)
      {
         strOldS = "TreeOldS" + IntToStr(j);
         strNewS = "TreeNewS" + IntToStr(j);
         strOldB = "TreeOldB" + IntToStr(j);
         strNewB = "TreeNewB" + IntToStr(j);

         treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
         treeNewS = (TTree*)nextsource->Get(strNewS.c_str());
         treeOldB = (TTree*)nextsource->Get(strOldB.c_str());
         treeNewB = (TTree*)nextsource->Get(strNewB.c_str());

         nrevts.push_back(treeOldS->GetEntries());
         nrevts.push_back(treeNewS->GetEntries());
         nrevts.push_back(treeOldB->GetEntries());
         nrevts.push_back(treeNewB->GetEntries());
      }

      nextsource = (TFile*)sourcelist->After(nextsource);
      fileCount++;
   }
}

void RewriteRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames )
{
   int sigKeys = 0;
   for(int i = 0; i < nrkeys.size(); i++)
      sigKeys += (nrkeys[i]-1)/4;
   printf("Number of signal keys in all files = %d\n", sigKeys);

   if(sigKeys == 1)
   {
      printf("Only one key found in root file. Aborting.\n");
      return;
   }

   int nrfiles = sourcelist->GetSize();
   string strOldS, strNewS, strOldB, strNewB;
   TTree *treeOldS, *treeNewS, *treeOldB, *treeNewB, *treeAll;
   TTree *outSig[sigKeys], *outBack[sigKeys], *outAll;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   printf("Get number of files: %d\n", nrfiles);

   TList *listAll, *listSig[sigKeys], *listBack[sigKeys];
   target->cd();
   listAll = new TList;
   for(int j = 0; j < sigKeys; j++)
   {
      listSig[j] = new TList;
      listBack[j] = new TList;
   }

   printf("Number of keys in file %s: %d (%d)\n", nextsource->GetName(), nrkeys[0], (nrkeys[0]-1)/4);
   
   for(int j = 1; j <= (nrkeys[0]-1)/4; j++)
   {
      strOldS = "TreeOldS" + IntToStr(j);
      strNewS = "TreeNewS" + IntToStr(j);
      strOldB = "TreeOldB" + IntToStr(j);
      strNewB = "TreeNewB" + IntToStr(j);

      nextsource->cd(path);

      treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
      treeNewS = (TTree*)nextsource->Get(strNewS.c_str());
      treeAll = (TTree*)nextsource->Get("TreeA");
      treeOldB = (TTree*)nextsource->Get(strOldB.c_str());
      treeNewB = (TTree*)nextsource->Get(strNewB.c_str());

      // Rewrite the signal trees
      printf("#- %d a ---------------------------------------------\n", j);

      if(treeOldS->GetEntries() > 0)
      {
         printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
         listSig[j-1]->Add(treeOldS);
      }

      if(treeNewS->GetEntries() > 0)
      {
         printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
         listSig[j-1]->Add(treeNewS);
      }

      // Rewrite the background trees
      printf("#- %d b ---------------------------------------------\n", j);

      if(treeOldB->GetEntries() > 0)
      {
         printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeOldB->GetEntries(), strOldB.c_str(), nextsource->GetKey(strOldB.c_str())->GetTitle());
         listBack[j-1]->Add(treeOldB);
      }

      if(treeNewB->GetEntries() > 0)
      {
         printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeNewB->GetEntries(), strNewB.c_str(), nextsource->GetKey(strNewB.c_str())->GetTitle());
         listBack[j-1]->Add(treeNewB);
      }

   }

   // Merge all events
   printf("#- c ---------------------------------------------\n");
   if(treeAll->GetEntries() > 0)
   {
      printf("Found %d entries in tree \"TreeA\" and title \"%s\"\n", (int)treeAll->GetEntries(), nextsource->GetKey("TreeA")->GetTitle());
      listAll->Add(treeAll);
   }

   printf("#- e saving the file --------------------------------------\n");

//   string treeTitle[sigKeys];

   target->cd();
   printf("Saving signal trees\n");
   for(int i = 0; i < sigKeys; i++)
   {
      outSig[i] = TTree::MergeTrees(listSig[i]);
      strOldS = "TreeS" + IntToStr(i+1);
      outSig[i]->SetName(strOldS.c_str());
//      treeTitle[i] = outSig[i]->GetTitle();
      outSig[i]->Write();
   }
   printf("Saving all events tree\n");
   outAll = TTree::MergeTrees(listAll);
   outAll->Write();
   printf("Saving background trees\n");
   for(int i = 0; i < sigKeys; i++)
   {
      outBack[i] = TTree::MergeTrees(listBack[i]);
      strOldS = "TreeB" + IntToStr(i+1);
      outBack[i]->SetName(strOldS.c_str());
//      outBack[i]->SetTitle((RemoveFromFirst(treeTitle[i], "Signal tree from", "Background tree without events from")).c_str());
      outBack[i]->Write();
   }

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}

void MergeRootfile( TDirectory *target, TList *sourcelist, vector<int> &nrkeys, vector<int> &nrevts, vector<string> &filenames )
{
   int sigKeys = 0;
   for(int i = 0; i < nrkeys.size(); i++)
      sigKeys += (nrkeys[i]-1)/4;
   printf("Number of signal keys in all files = %d\n", sigKeys);

   int nrfiles = sourcelist->GetSize();
   string strOldS, strNewS, strOldB, strNewB;
   TTree *treeOldS, *treeNewS, *treeOldB, *treeNewB, *treeAll;
   TTree *outSig[sigKeys], *outBack[sigKeys], *outAll;
   int selectCount[2];
   int fileCount = 0;

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove(0,2);

   TFile *nextsource = (TFile*)sourcelist->First();
   nextsource->cd(path);
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   printf("Get number of files: %d\n", nrfiles);

   fileCount = 0;
   selectCount[0] = 0;
   selectCount[1] = 0;

   TList *listAll, *listSig[sigKeys], *listBack[sigKeys];
   target->cd();
   listAll = new TList;
   for(int j = 0; j < sigKeys; j++)
   {
      listSig[j] = new TList;
      listBack[j] = new TList;
   }

   while(nextsource)
   {
      printf("Number of keys in file %s: %d (%d)\n", nextsource->GetName(), nrkeys[fileCount], (nrkeys[fileCount]-1)/4);
      
      for(int j = 1; j <= (nrkeys[fileCount]-1)/4; j++)
      {
         strOldS = "TreeOldS" + IntToStr(j);
         strNewS = "TreeNewS" + IntToStr(j);
         strOldB = "TreeOldB" + IntToStr(j);
         strNewB = "TreeNewB" + IntToStr(j);

         nextsource->cd(path);

         treeOldS = (TTree*)nextsource->Get(strOldS.c_str());
         treeNewS = (TTree*)nextsource->Get(strNewS.c_str());
         treeAll = (TTree*)nextsource->Get("TreeA");
         treeOldB = (TTree*)nextsource->Get(strOldB.c_str());
         treeNewB = (TTree*)nextsource->Get(strNewB.c_str());

         // Merge the signal trees
	 printf("#- %d %d a ---------------------------------------------\n", fileCount, j);

         if(treeOldS->GetEntries() > 0)
	 {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
	    listSig[selectCount[0]]->Add(treeOldS);
	    selectCount[0]++;
	 }

         if(treeNewS->GetEntries() > 0)
	 {
            printf("Found %d entries in tree \"%s\" and title \"%s\"\n", (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
	    listSig[selectCount[0]]->Add(treeNewS);
	    selectCount[0]++;
	 }

         // Merge the background trees
	 printf("#- %d %d b ---------------------------------------------\n", fileCount, j);

         if( (nrkeys[fileCount]-1)/4 == 1 )
	 {
	    printf("No background trees: selectCount[1] = %d\n", selectCount[1]);
	    if(fileCount == 0)
	    {
	       for(int i = 1; i < sigKeys; i++)
	       {
                  if(treeOldS->GetEntries() > 0)
	          {
                     printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
	             listBack[i]->Add(treeOldS);
	          }

                  if(treeNewS->GetEntries() > 0)
	          {
                     printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
	             listBack[i]->Add(treeNewS);
	          }
	       }
	    }
	    else
	    {
	       for(int i = 0; i < sigKeys-1; i++)
	       {
                  if(treeOldS->GetEntries() > 0)
	          {
                     printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldS->GetEntries(), strOldS.c_str(), nextsource->GetKey(strOldS.c_str())->GetTitle());
	             listBack[i]->Add(treeOldS);
	          }

                  if(treeNewS->GetEntries() > 0)
	          {
                     printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewS->GetEntries(), strNewS.c_str(), nextsource->GetKey(strNewS.c_str())->GetTitle());
	             listBack[i]->Add(treeNewS);
	          }
	       }
	    }
	 }
	 else
	 {
	    printf("Some background trees: selectCount[1] = %d\n", selectCount[1]);
	    if(fileCount == 0)
	    {
	       selectCount[1] = j-1;
	       printf("selectCount = %d\n", selectCount[1]);

	       for(int i = selectCount[1]; i < sigKeys-((nrkeys[fileCount+1]-1)/4); i++)
	       {
                  if( selectCount[1] == i )
		  {
                     if(treeOldB->GetEntries() > 0)
		     {
                        printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldB->GetEntries(), strOldB.c_str(), nextsource->GetKey(strOldB.c_str())->GetTitle());
	                listBack[i]->Add(treeOldB);
	             }

                     if(treeNewB->GetEntries() > 0)
		     {
                        printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewB->GetEntries(), strNewB.c_str(), nextsource->GetKey(strNewB.c_str())->GetTitle());
	                listBack[i]->Add(treeOldB);
		     }
		  }
		  else
	            selectCount[1]++;
	       }

	       for(int i = sigKeys-((nrkeys[fileCount+1]-1)/4); i < sigKeys; i++)
	       {
                  if(treeOldB->GetEntries() > 0)
		  {
                     printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldB->GetEntries(), strOldB.c_str(), nextsource->GetKey(strOldB.c_str())->GetTitle());
	             listBack[i]->Add(treeOldB);
		  }

                  if(treeNewB->GetEntries() > 0)
		  {
                     printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewB->GetEntries(), strNewB.c_str(), nextsource->GetKey(strNewB.c_str())->GetTitle());
	             listBack[i]->Add(treeNewB);
		  }
	       }
	    }
	    else
	    {
	       selectCount[1] = ((nrkeys[fileCount-1]-1)/4)+j-1;

	       for(int i = 0; i < ((nrkeys[fileCount-1]-1)/4); i++)
	       {
                  if(treeOldB->GetEntries() > 0)
		  {
                     printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldB->GetEntries(), strOldB.c_str(), nextsource->GetKey(strOldB.c_str())->GetTitle());
	             listBack[i]->Add(treeOldB);
		  }

                  if(treeNewB->GetEntries() > 0)
		  {
                     printf("Early, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewB->GetEntries(), strNewB.c_str(), nextsource->GetKey(strNewB.c_str())->GetTitle());
	             listBack[i]->Add(treeNewB);
		  }
	       }

	       for(int i = ((nrkeys[fileCount-1]-1)/4); i < sigKeys; i++)
	       {
                  if( selectCount[1] == i )
		  {
                     if(treeOldB->GetEntries() > 0)
		     {
                        printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeOldB->GetEntries(), strOldB.c_str(), nextsource->GetKey(strOldB.c_str())->GetTitle());
	                listBack[i]->Add(treeOldB);
		     }

                     if(treeNewB->GetEntries() > 0)
		     {
                        printf("Late, i = %d, selectCount = %d: Found %d entries in tree \"%s\" and title \"%s\"\n", i, selectCount[1], (int)treeNewB->GetEntries(), strNewB.c_str(), nextsource->GetKey(strNewB.c_str())->GetTitle());
	                listBack[i]->Add(treeOldB);
	   	     }
		  }
	       }
	    }
	 }

      }

      // Merge all events
      printf("#- %d c ---------------------------------------------\n", fileCount);
      if(treeAll->GetEntries() > 0)
      {
         printf("Found %d entries in tree \"TreeA\" and title \"%s\"\n", (int)treeAll->GetEntries(), nextsource->GetKey("TreeA")->GetTitle());
	 listAll->Add(treeAll);
      }

      // Move to the next file
      printf("#- %d d -----------------------------------------------\n", fileCount);
      fileCount++;
      nextsource = (TFile*)sourcelist->After(nextsource);
   }

   printf("#- e saving the file --------------------------------------\n");

   string treeTitle[sigKeys];

   target->cd();
   printf("Saving signal trees\n");
   for(int i = 0; i < sigKeys; i++)
   {
      outSig[i] = TTree::MergeTrees(listSig[i]);
      strOldS = "TreeS" + IntToStr(i+1);
      outSig[i]->SetName(strOldS.c_str());
      treeTitle[i] = outSig[i]->GetTitle();
      outSig[i]->Write();
   }
   printf("Saving all events tree\n");
   outAll = TTree::MergeTrees(listAll);
   outAll->Write();
   printf("Saving background trees\n");
   for(int i = 0; i < sigKeys; i++)
   {
      outBack[i] = TTree::MergeTrees(listBack[i]);
      strOldS = "TreeB" + IntToStr(i+1);
      outBack[i]->SetName(strOldS.c_str());
      outBack[i]->SetTitle((RemoveFromFirst(treeTitle[i], "Signal tree from", "Background tree without events from")).c_str());
      outBack[i]->Write();
   }


   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}
