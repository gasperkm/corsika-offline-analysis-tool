#ifndef _COMBINE_H_
#define _COMBINE_H_

#include <vector>
#include <string.h>

void MergeRootfile( TDirectory *target, TList *sourcelist, std::vector<int> &nrkeys, std::vector<int> &nrevts, std::vector<std::string> &filenames );
void RewriteRootfile( TDirectory *target, TList *sourcelist, std::vector<int> &nrkeys, std::vector<int> &nrevts, std::vector<std::string> &filenames );
void CheckKeys( TDirectory *target, TList *sourcelist, std::vector<int> &nrkeys, std::vector<int> &nrevts );
void hadd(int nrfiles, char **files);

#endif
