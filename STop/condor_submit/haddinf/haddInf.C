/******************************************************************************
 * This is a evolved from the macro $ROOTSYS/tutorial/io/hadd.C 
 * which is provided as a tutorial.
 * I have added several very useful features.
 * 1. Merge/Add histograms/trees with weights.
 *		For this, a weight should be specified for each of the input root files.
 * 2. Command line execution. No need to load root every time.
 *
 * @Author: Samantha Hewamanage <samantha@fnal.gov>
 * @Employer: Florida International University
 *****************************************************************************/

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "TProfile.h"
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

TFile *Target;
typedef vector<TFile*> tfileVec;
typedef vector<TFile*>::const_iterator tfileVec_it;

void MergeRootfile( TDirectory *target, const tfileVec& vFileList);

void CopyVec(const vector<string>& src, tfileVec& des, const int& x, const int&y)
{
	if (x>=0 && x< (int) src.size() && y>=0 && y<(int) src.size() && x<=y)
	{
		for (int i=x; i<=y; ++i) des.push_back(TFile::Open(src.at(i).c_str()));
	} else {
		cout << __FUNCTION__ << ":precondition fail! x,y, n=" << x << "/" << y << "/" << src.size() << endl;
		assert(false);
	}
}


// entry point for commandline executable
int main(int argc, char* argv[])
{

	if (argc<2)
	{
		cout << "\n======================== USAGE INFO ===============================================" << endl;
		cout << "       ./haddinf file1.root file2.root w1 w2\n" << endl;
		cout << "       w1, w2: are floating point numbers." << endl;
		cout << "               Each input file should have a corresponding weight (>0)." << endl;
		cout << "               A weight of 1 will leave those histograms unchanged." << endl;
		cout << "       Files names and weights can be given in any order." << endl;
		cout << "          But the first weight found will be assigned to the first file listed," << endl;
		cout << "          second weight found will be assigned to the second file listed and so on." << endl;
		cout << "       Output file will be named results.root and will be recreated with subseqent" << endl; 
		cout << "          executions." << endl;
		cout << "======================== END USAGE INFO =============================================\n" << endl;
		return 0;
	}

	tfileVec inputList;
	vector<string> inputFileNames;
	vector<float> inputWeights; 

	try {
		for (int i=1; i< argc; ++i)
		{
			//const char *argi = argv[i];
			string sargi(argv[i]);

			//check if this a number
			if (sargi.find_first_not_of("1234567890.") == string::npos) //not a number
			{
				const float w = atof(argv[i]);
				inputWeights.push_back(w);
			} else {
				if (sargi == "Merged.root")
				{
					cout << "A file name with Merged.root alrady exists! Please rename or it will be overidden!" << endl;
					return 1;
				}
				inputFileNames.push_back(sargi);
			}
		}  
	} catch (exception& e)
	{
		cout << e.what() << endl;
	}

	//do some basic sanity checks
	//if no weights are given, just add them. else use the weights.
	if (inputWeights.size() > 0) 	
	{
		if (inputFileNames.size() != inputWeights.size())
		{
			cout << "Every input root file must have a corresponding weight! please check!" << endl;
			return 0;
		}
	}
	//check if all files exists and readable
	for (vector<string>::const_iterator it = inputFileNames.begin(); it != inputFileNames.end(); ++it)
	{
		TFile *f = new TFile (it->c_str());
		if (f->IsZombie()) {
			cout << "File: " << (*it) << " not found or readable!" << endl;
			return 0;
		}
		delete f;
	}

	tfileVec vFileList;

	const int perMerge = 500; 
	cout << "NinputFiles = " << inputFileNames.size() << endl;
	int d =  inputFileNames.size() / perMerge;
	int r =  inputFileNames.size() % perMerge;
	int nTempFiles = d + 1;
	cout << "d/r/nTempFile = " << d << "/" << r << "/" << nTempFiles << endl;
	tfileVec tempFiles;  

	for (int t = 1 ; t<=nTempFiles; t++)
	{
		tfileVec vTempFiles;
		const int begin = (t-1)*perMerge;
		int upTo  =  perMerge-1;
		if (t==nTempFiles) upTo = r-1;
		const int end   = begin + upTo;
		cout << __LINE__ <<":t/begin/end(upTo)=" << t << "/" << begin << "/" << end  << " ("<< upTo <<")" << endl;
		CopyVec(inputFileNames, vTempFiles, begin, end );
		stringstream outfile;
		outfile << "Temp_" << t << ".root";
		
   	if (vTempFiles.size()>1) 
		{
			Target = new TFile(outfile.str().c_str(), "RECREATE");
			cout << "Created temp file " << Target->GetName()  << " and adding these:" << endl;
			for (unsigned c=0;c<vTempFiles.size();c++) cout << "\t " << vTempFiles.at(c)->GetName() <<endl;
			MergeRootfile(Target , vTempFiles );
		} else 
		{
			stringstream cmd;
			cmd << "cp -v " << vTempFiles.at(0)->GetName() << " " << outfile.str();
			system(cmd.str().c_str());
			Target = new TFile(outfile.str().c_str());
			cout << "Created temp file " << Target->GetName()  << " and adding these:" << endl;
			for (unsigned c=0;c<vTempFiles.size();c++) cout << "\t " << vTempFiles.at(c)->GetName() <<endl;
		}
		for (unsigned c=0;c<vTempFiles.size();c++) delete vTempFiles.at(c);
		tempFiles.push_back(Target);
	}

	cout << "tempFile = " << tempFiles.size() << endl;

	if (tempFiles.size()>1) 
	{
		Target = TFile::Open("Merged.root", "RECREATE");
		cout << "Created final file " << Target->GetName() << " and adding these: " << endl;
		//Since all intermediate temp files are "RECREATED", need to close them and open in "READ" mode 
		for (unsigned c=0;c<tempFiles.size();c++) 
		{
			//cout << "\t " << tempFiles.at(c)->GetName() <<endl;
			tempFiles.at(c)->Print();
			const string name(tempFiles.at(c)->GetName());
			tempFiles.at(c)->Close();
			tempFiles.at(c) = new TFile(name.c_str());
			//cout << "\t->"; tempFiles.at(c)->Print();
		}
		MergeRootfile(Target , tempFiles);
		//Target->Close();
		for (unsigned c=0;c<tempFiles.size();c++) 
		{
			const string name(tempFiles.at(c)->GetName());
			delete tempFiles.at(c);
			stringstream cmd;
			cmd << "rm -v " << name;
			system(cmd.str().c_str());
		}

	} else {
		system("mv -v Temp_1.root Merged.root"); 
	}

	//now cleanup. This makes valgrind happy :)
	for (tfileVec_it it = vFileList.begin(); it != vFileList.end(); ++it)
	{
		//delete it;
	}

	for (tfileVec_it it = tempFiles.begin(); it != tempFiles.end(); ++it)
	{
		//delete it;
	}

	cout << "Added " << inputFileNames.size() << " files to " << Target->GetName() << endl;
	Target->Close();

	return 1;
}

void haddInf() {
   // in an interactive ROOT session, edit the file names
   // Target and FileList, then
   // root > .L haddInf.C
   // root > haddInf()

   Target = TFile::Open( "result.root", "RECREATE" );

	tfileVec vFileList;
	tfileVec_it it;

	const float xSec[] = {
		276000.0,  //pt 
		8426.0,  //
		204.0  //
	};

	const float nEvents[] = {
		26966344, // #numbers from DBS, PREP page numbers are approximate
		30529026, // 
		13823375 // 
	};
	//const float scaleTo = 10000.0; //pb-1
	const float scaleTo = 19455.518; //A+B+C+D pb-1
	const int PTbins = 3;
	//double w[PTbins];

	for (int i=0; i<PTbins; ++i)
	{
		//w[i] = scaleTo * xSec[i]/nEvents[i];
	}

	cout << "Scaling to Lumi of " << scaleTo << endl;
//	vFileList.push_back(make_pair(TFile::Open("MG_QCD_HT_250To500.root"),w[0]));
//	vFileList.push_back(make_pair(TFile::Open("MG_QCD_HT_500To1000.root"),w[1]));
//	vFileList.push_back(make_pair(TFile::Open("MG_QCD_HT_1000ToInf.root"),w[2]));
	
	for (it = vFileList.begin(); it != vFileList.end(); ++it)
	{
	//	cout << "File/weight = " << it->first->GetName() << "/" << it->second << endl;
	}

   //MergeRootfile( Target, vFileList );

	//now cleanup. This makes valgrind happy :)
	for (it = vFileList.begin(); it != vFileList.end(); ++it)
	{
	//	delete (it->first);
	}

	Target->Close();

}

void MergeRootfile( TDirectory *target, const tfileVec& vFileList) {

   //  cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );
	
	tfileVec_it it = vFileList.begin();
   TFile *first_source = (*it);
   const float first_weight = 1;
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         //      cout << "Merging histogram " << obj->GetName() << endl;
         TH1 *h1 = (TH1*)obj;
			if (first_weight>0) {
				h1->Scale(first_weight);
			}

         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
			for (tfileVec_it nextsrc = vFileList.begin()+1; nextsrc != vFileList.end(); ++nextsrc ) {

            // make sure we are at the correct directory level by cd'ing to path
            (*nextsrc)->cd( path );
				const float next_weight = 1;
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
               TH1 *h2 = (TH1*)key2->ReadObj();
					if (next_weight>0) h2->Scale(next_weight);
               h1->Add( h2 );
               delete h2;
            }
         }

      } else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

         // loop over all source files create a chain of Trees "globChain"
         const char* obj_name= obj->GetName();

         globChain = new TChain(obj_name);
         globChain->Add(first_source->GetName());
			for (tfileVec_it nextsrc = vFileList.begin()+1; nextsrc != vFileList.end(); ++nextsrc ) {
            globChain->Add((*nextsrc)->GetName());
         }

      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory

         cout << "Found subdirectory " << obj->GetName() << endl;

         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

         // newdir is now the starting point of another round of merging
         // newdir still knows its depth within the target file via
         // GetPath(), so we can still figure out where we are in the recursion
         MergeRootfile( newdir, vFileList);

      } else {

         // object is of no type that we know or can handle
         cout << "Unknown object type, name: "
         << obj->GetName() << " title: " << obj->GetTitle() << endl;
      }

      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      if ( obj ) {
         target->cd();

         //!!if the object is a tree, it is stored in globChain...
         if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
         else
            obj->Write( key->GetName() );
      }

   } // while ( ( TKey *key = (TKey*)nextkey() ) )

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}
