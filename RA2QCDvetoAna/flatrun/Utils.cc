#include "Utils.hh"
#include <iostream>
//#include "Tokenizer.h"
#include<vector>
#include<string>
/*hide boost header files from CINT */
/*#ifndef __CINT__
include <boost/foreach.hpp>
include <boost/tokenizer.hpp>
#endif
*/

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
