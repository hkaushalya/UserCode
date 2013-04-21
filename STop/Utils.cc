#include "Utils.hh"
#include <iostream>
//#include "Tokenizer.h"
#include<vector>
#include<string>
#include <iomanip>
#include "assert.h"
#include <stdio.h>
/*hide boost header files from CINT */
/*#ifndef __CINT__
include <boost/foreach.hpp>
include <boost/tokenizer.hpp>
#endif
*/

using namespace std;
bool PtAComparator(const TLorentzVector& a,const TLorentzVector& b)
{
	return a.Pt() > b.Pt();
}
/*vector<double> NumStringToVec(const string s)
{
	cout << __FUNCTION__ << ":" << __LINE__ << endl;
	cout << "string got " << s << endl;
	vector<double> res;
	
    // instanciate Tokenizer class
    Tokenizer str;
    string token;
    int counter = 0;

    // set source string with default delimiters (space, tab, and newline char)
    str.set(s);

    // Tokenizer::next() returns a next available token from source string
    // If it reaches EOS, it returns zero-length string, "".
    while((token = str.next()) != "")
    {
        ++counter;
        cout << counter << ": " << token << endl;
		  //string part(token);
		  double num = atof(token.c_str());
		  res.push_back(num);
    }
    cout << endl;

	cout << __FUNCTION__ << ":" << __LINE__ << endl;
	return res;
}
*/
void DumpHist(const TH1* h)
{
	if (h == NULL)
	{
		cout << __FUNCTION__ << ": Nul pointer received! " << endl;
		assert(false);
	}
	if (h->GetDimension() != 1)
	{
		cout << __FUNCTION__ << ": Historam dimenstion is not 1! returning." << endl;
		return;
	}

	cout << __FUNCTION__ << ":: Bin contents for " << h->GetName()  << " with Nbins = " << h->GetNbinsX() << endl;

	cout << setw(25) << " bin " << setw(15) << "val" << setw(15) << "err" << endl;	
	for (int bin=0; bin <= h->GetNbinsX()+1; ++bin)
	{
		const double lox = h->GetBinLowEdge(bin);
		const double hix = h->GetXaxis()->GetBinUpEdge(bin);
		const double w   = h->GetBinWidth(bin);
		const double val = h->GetBinContent(bin);
		const double err = h->GetBinError(bin);
		cout << bin << "[" << lox << "-" << hix << "] = " << setw(15) << val << setw(15) << err << endl;  
	}

}
/* prints a nice progress bar */
void printProgBar( int percent ){
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}
