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

  //  double a_powers[] = {0,1,2};
  double a_coeff[]  = {0,0,1};
  //  vector<double> powers(a_powers, a_powers+3); 
  vector<double> coeff(a_coeff, a_coeff+3);


  Polynomial quad(coeff, 1, 2);
  
  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    double yy = quad.eval(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << yy << endl;
  }

  cout << "#derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);    
    double yy = quad.grad1(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << yy << endl;
  }
  

  return 0;
}
