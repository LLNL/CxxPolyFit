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

  double a_coeff[]  = {0,0,0,1,0,1};
  vector<double> coeff(a_coeff, a_coeff+6);

  Polynomial quad(coeff, 2, 2);
  
  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.eval(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }

  cout << "# X derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivX(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << 0 << " " << yy << endl;
  }

  cout << "# Y derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.derivY(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << 0 << " " << ii << " " << yy << endl;
  }

  return 0;
}
