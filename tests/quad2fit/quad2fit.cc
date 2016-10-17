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

  double a_xpts[]  = {-2,0,
                      -1,0,
                      0, 0,
                      1, 0,
                      2, 0,
                      0,-2,
                      0,-1,
                      0,-0,
                      0, 1,
                      0, 2};
  double a_ypts[]  = {4,1,0,1,4, 
                      4,1,0,1,4};
  vector<double> xpts(a_xpts, a_xpts+20);
  vector<double> ypts(a_ypts, a_ypts+10);


  Polynomial quad(xpts, ypts, 2, 2);
  
  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.eval(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << ii << " " << yy << endl;
  }

  cout << "#derivative" << endl;

  for(double ii = -10; ii < 11; ++ii) {
    vector<double> xx;
    xx.push_back(ii);
    xx.push_back(ii);
    double yy = quad.grad1(xx);
    if(abs(yy) < .0001) { yy = 0; }  //Just avoiding a common thing that kills tests on different platforms
    cout << ii << " " << 0 << " " << yy << endl;
  }


  return 0;
}
