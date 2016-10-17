//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by Jim Leek <leek2@llnl.gov> and Adrian Wong <adrianskw@gmail.com>
// LLNL-CODE-704097
// All rights reserved.
// This file is part of C++ PolyFit. 
// For details, see https://github.com/llnl/CxxPolyFit.
// Please also CxxPolyFit/LICENSE.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include <cstdlib>
#include <iostream>
#include <string>
#include "Polynomial.hh"

using namespace std;

int main() {

  double a_coeff[]  = { 6, 7, 8, 4,
		        5, 7, 5,10,
		        6, 3 };/* ,32,12,
		       54,12,54, 1,
		       5, 8, 5, 3};*/
  vector<double> coeff(a_coeff,a_coeff+10); // find a way to assign the a_coeff automatically without specifying length

  Polynomial quad(coeff, 2,3);
  
  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.eval(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }

  cout << "# XX derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivXX(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }

  cout << "# YY derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivYY(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }

  cout << "# XY derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivXY(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }


  cout << "# YX derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivYX(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }


  return 0;
  
}
