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

#include <Polynomial.hh>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

extern "C"
{
  /** All fortran arugments are passed by reference */
  extern void dgelss_(int* M, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, double* S, double* RCOND, int* RANK, double* WORK, int* LWORK, int* INFO);
}


/**
 * This constructor makes a Polynomial from a set of data points via regression.
 *  xs is a multidimensional array of the axes of the polynomial points.  
 *   ie, if we are fitting a 2D polynomial, xs will be of length numpts * 2.  
 *   The points are given as sets.  For a 2D polynomial, as pairs.
 *
 *  ys is a single dimensional array of the values at the xs points.  
 *   So it is of length numpts.
 *
 *  in_dims says how many dimensions we are fitting.  In our example above, 
 *   in_dims is 2.  1, 2, and 3 are the only allowed values.
 *
 *  in_order is the order of the polynomial.  5 is often a good number, as 
 *   the order goes higher the algorithm seems to become less stable.
 **/
Polynomial::Polynomial(const vector<double>& xs, const vector<double>& ys, int in_dims, int in_order) : dims(in_dims), order(in_order) {

  // construct least-squares matrix A
  getBasisPowers(dims, order, powers);
  numBasis = getNumBasis(dims, order); //AKA n
  int numPoints = ys.size();  //AKA m
  double* vals = new double[numBasis];
  double* A = new double[numBasis * numPoints];  //m x n n is number basis funcs, m number of points
  for(int ii = 0; ii < numPoints; ++ii) {
    vector<double> i_xs(xs.begin()+(dims*ii), xs.begin()+(dims*ii)+dims);  //Get the function values as calculated with existing basis functions for each point.
    vector<double> vals;
    basisEvals(i_xs,vals); //get all basis function values at each point
    
    //Since A is passed into LAPack, it must be in column major order, but it's the only thing that is.
    //So I just have this special copy here
    for(int jj = 0; jj < numBasis; ++jj) {
      A[jj*numPoints+ii] = vals[jj]; 
    }
  }
  // find least-squares solution of A c = f
  int nrhs  = 1;
  int work_size = 6 * max(numPoints,numBasis);  //allocate work space
  double* Awork = new double[numBasis * numPoints];
  double* Bwork = new double[nrhs * numPoints];
  double* S     = new double[numBasis];
  double* work  = new double[work_size];

  copy(A, A + (numBasis * numPoints), Awork); // don't overwrite A, Bwork is a single column of ys
  copy(ys.begin(), ys.end(), Bwork);

  double rcond = -1;
  int rank = 0;
  int info = 0;
  dgelss_(&numPoints,&numBasis,&nrhs,Awork,&numPoints,Bwork,&numPoints,S,&rcond,&rank,work,&work_size, &info);
  
  coefficients.resize(numBasis);
  copy(Bwork, Bwork+numBasis, coefficients.begin());   //You get the coefficients out of B[1:n] from column 1

}

/*
Polynomial::Polynomial(const vector<double>& in_powers, const vector<double>& in_coefficients, double in_dims, double in_order) : dims(in_dims), order(in_order) {
  numBasis = getNumBasis(dims, order);
  
  copy(in_powers.begin(), in_powers.end(), back_inserter(powers));
  copy(in_coefficients.begin(), in_coefficients.end(), back_inserter(coefficients));

}
*/

/**
 * Makes a polynomial by directly setting the coefficents, dimensions, and order. 
 * You have to know how the coefficents array is laid out to use this constructor.
 **/

Polynomial::Polynomial(const vector<double>& in_coefficients, double in_dims, double in_order) : dims(in_dims), order(in_order) {
  numBasis = getNumBasis(dims, order);
  getBasisPowers(dims, order, powers);
  copy(in_coefficients.begin(), in_coefficients.end(), back_inserter(coefficients));

}

Polynomial::Polynomial(const Polynomial& rhs) : dims(rhs.dims), order(rhs.order), numBasis(rhs.numBasis) {
  
  copy(rhs.powers.begin(), rhs.powers.end(), back_inserter(powers));
  copy(rhs.coefficients.begin(), rhs.coefficients.end(), back_inserter(coefficients));

}

//This is rather a silly algorithm, but it works and only happens rarely
int Polynomial::getNumBasis(int in_dims, int in_order) {

  double outNumBasis = 0;
  switch(in_dims) {

  case 1:
    return in_order+1;

  case 2:
    for(int ii = 0; ii < in_order+1; ++ii) {
      for(int j1 = in_order; j1 >= 0; --j1) {
	for(int j2 = in_order; j2 >= 0; --j2) {
	  if (j1+j2==ii) ++outNumBasis;
	}
      }
    }
    return outNumBasis;

  case 3:
    for(int ii = 0; ii < in_order+1; ++ii) {
      for(int j1 = in_order; j1 >= 0; --j1) {
	for(int j2 = in_order; j2 >= 0; --j2) {
	  for(int j3 = in_order; j3 >= 0; --j3) {
	    if (j1+j2+j3==ii) ++outNumBasis;
	  }
	}
      }
    }
    return outNumBasis;
  default:
    throw new invalid_argument("Invalid dimension passed to getNumBasis, only 1, 2, or 3 are valid values");
  }
  return 0;

}

//This is rather a silly algorithm, but it works and only happens rarely
void Polynomial::getBasisPowers(int in_dims, int in_order, vector<double>& out_powers) {

  double numBasis = getNumBasis(in_dims, in_order);
  double basisNum = 0;
  out_powers.resize(in_dims*numBasis);
  fill(out_powers.begin(), out_powers.end(), 0);

  switch(in_dims) {
  case 1:
    for(int ii = 0; ii < in_order+1; ++ii) {
      out_powers[ii] = ii;
    }
    return;
  case 2:
    for(int ii = 0; ii < in_order+1; ++ii) {
      for(int j1 = in_order; j1 >= 0; --j1) {
	for(int j2 = in_order; j2 >= 0; --j2) {	  
	  if (j1+j2==ii) {
	    out_powers[basisNum] = j1;             //1st dimension
	    out_powers[numBasis + basisNum] = j2;  //2nd dimension
	    ++basisNum;
	  }
	}
      }
    }
    return;
  case 3:
    for(int ii = 0; ii < in_order+1; ++ii) {
      for(int j1 = in_order; j1 >= 0; --j1) {
	for(int j2 = in_order; j2 >= 0; --j2) {
	  for(int j3 = in_order; j3 >= 0; --j3) {
	    if (j1+j2+j3==ii) {
	      out_powers[basisNum] = j1;                 //1st dimension
	      out_powers[numBasis + basisNum] = j2;      //2nd dimension
	      out_powers[(numBasis*2) + basisNum] = j3;  //3rd dimension
	      ++basisNum;
	    }
	    
	  }
	}
      }
    }
    return;

  default:
    throw new invalid_argument("Polynomial::getBasisPowers: Invalid dimension only 1, 2, or 3 are valid values");
  }


}

//Evaluate just the basis functions
void Polynomial::basisEvals(const vector<double>& xs, vector<double>& outVals) {
  outVals.resize(numBasis);
  fill(outVals.begin(), outVals.end(), 1);
  
  for(int bb = 0; bb < numBasis; ++bb) {
    for(int dd = 0; dd < dims; ++dd) {
      if(powers[numBasis*dd + bb] > 0) {  //-1 means the power doesn't have that term
	outVals[bb] *= pow(xs[dd],powers[numBasis*dd + bb]);
      } else { //if(powers[numBasis*dd + bb] == 0) {  a 0th power is always 1, so just leave it as 1
        ;//outVals[bb] = xs[dd];  
      }
    }
  }
}

//Does the coefficients multiplcations and the summation
double Polynomial::eval(const vector<double>& xs) {

  if(xs.size() != dims) {
    throw new invalid_argument("Polynomial::eval: xs must be the same length as polynomial dimensions.");
  }
  double retVal = 0;
  vector<double> basisVals;
  basisEvals(xs, basisVals);
  for(int bb = 0; bb < numBasis; ++bb) {
    retVal += coefficients[bb] * basisVals[bb];
  }
  return retVal;
}

void Polynomial::basisGradients(const vector<double>& xs, vector<double>& outGrads) {
  outGrads.resize(dims*numBasis);
  fill(outGrads.begin(), outGrads.end(), 1);

  for(int bb = 0; bb < numBasis; ++bb) {    // for each basis function ...
    for(int ii = 0; ii < dims; ++ii) {      // for each dimension ...
      for(int nn = 0; nn < dims; ++nn) {    // multiply each factor ...
	if (nn==ii)  {                      // multiply derivative factor
	  if (powers[numBasis*nn +bb] > 1) {
	    outGrads[numBasis*nn +bb] *= powers[numBasis*nn +bb] * pow(xs[nn], powers[numBasis*nn + bb]-1);
	  } else if (powers[numBasis*nn +bb] == 1) {
	    ; //Multiply by 1
	  } else {
	    outGrads[numBasis*nn +bb] = 0; 
	  }
	} else {              // multiply non-derivative factor
	  if (powers[numBasis*nn +bb] > 0) {
	    outGrads[numBasis*ii +bb] *= pow(xs[nn],powers[numBasis*nn +bb]);
	  } 
	}
      }
    }
  }
}

void Polynomial::gradAllDims(const vector<double>& xs, vector<double>& outGrads) {
  outGrads.resize(dims);
  fill(outGrads.begin(), outGrads.end(), 0);
  
  vector<double> basisGrads;
  basisGradients(xs, basisGrads);
  for(int jj = 0; jj < numBasis; jj++) {
    for(int dd = 0; dd < dims; ++dd) {
      outGrads[dd] += coefficients[jj]*basisGrads[dd*numBasis+jj];
    }
  }
}

/** Just a convenience function for getting the 1st derivative of the 1st dimension
 */
double Polynomial::grad1(const vector<double>& xs) {
  vector<double> outGrads;
  gradAllDims(xs, outGrads);

  return outGrads[0];

}

/** Just a convenience function for getting the 1st derivative of the 2nd dimension
 */
double Polynomial::grad2(const vector<double>& xs) {
  vector<double> outGrads;
  gradAllDims(xs, outGrads);

  return outGrads[1];

}

/** Just a convenience function for getting the 1st derivative of the 3rd dimension
 */
double Polynomial::grad3(const vector<double>& xs) {
  vector<double> outGrads;
  gradAllDims(xs, outGrads);

  return outGrads[2];
}

 /** Returns a string giving all the coefficients and their basis functions.  
   *  Helpful for debugging.  Limited to 3 dimensions (like everything.) 
   **/
string Polynomial::termsToString() {
  char varnames[3] = {'x', 'y', 'z'};
  std::stringstream termStr;
  for(int basis = 0; basis < numBasis; ++basis) {
    termStr << coefficients[basis] << " ";
    bool firstVar = true;
    for(int dd = 0; dd < dims; ++dd) {
      if(powers[dd*numBasis + basis] > 0) {  //If powers is 0 for this variable, it isn't used in this basis function
        if(!firstVar) {
          termStr << " * ";  //If we have multiple vars, they get multiplied together
        } else {
          firstVar = false;
        }
        termStr << varnames[dd];
        if(powers[dd*numBasis + basis] > 1) {  //First powers don't need an exponent
          termStr << "^" << powers[dd*numBasis + basis];
        }
      }
    }
    termStr << "\n";
  }
  return termStr.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// New Deriv Stuff
// ---------------
// Since derivatives are linear operators, we can apply them as many times as we like, in any order.
// They are written here in a very general way that allows us combine them however we need.
// For example, d2f/dx2 = dx(dx(...)), d2f/dy2 = dy(dy(...)), and so on.
////////////////////////////////////////////////////////////////////////////////////////////////////
// numBasis    = number of bais functions (aka "terms"), eg: [x^2 + x + 0] has 3 basis functions
// coefficient = the coefficients in front of every basis function
// powers      = (flattened) array of the powers of each basis function, 
//                eg: powers[n]=powers[n+numBasis]=1 means the n-th term is x^1*y^1 
// xs          = the coordinates: xs[0] = x, xs[1] = y
////////////////////////////////////////////////////////////////////////////////////////////////////
 
// dx modifies new_coefficients and new_powers such that it is a first derivative of x
void Polynomial::dx(vector<double>& new_coefficients, vector<double>& new_powers){
  int ndims = new_powers.size()/new_coefficients.size();
  for (int b=0;b<numBasis;b++) {
    if (new_coefficients[b]!=0 && new_powers[b] > 0){
      new_coefficients[b] = new_coefficients[b]*new_powers[b];
      new_powers[b] = new_powers[b]-1;                         // x changes
      if (ndims > 1) {
	new_powers[b+numBasis] = new_powers[b+numBasis];       // y
	if (ndims == 3) { 
	  new_powers[b+2*numBasis] = new_powers[b+2*numBasis]; // z
	}
      }
    } else { // we set everything to zero
      new_coefficients[b] = 0;
      new_powers[b] = 0;
      if (ndims > 1) {
	new_powers[b+numBasis] = 0;
	if (ndims == 3) {
	  new_powers[b+2*numBasis] = 0;
	}
      }
    }
  }
}

// dy modifies new_coefficients and new_powers such that it is a first derivative of y
void Polynomial::dy(vector<double>& new_coefficients, vector<double>& new_powers){
  int ndims = new_powers.size()/new_coefficients.size();
  for (int b=0;b<numBasis;b++) {
    if (new_coefficients[b]!=0 && new_powers[b+numBasis] > 0){
      new_coefficients[b] = new_coefficients[b]*new_powers[b+numBasis];
      new_powers[b] = new_powers[b];                         // x
      new_powers[b+numBasis] = new_powers[b+numBasis]-1;     // y changes
      if (ndims == 3) { 
	new_powers[b+2*numBasis] = new_powers[b+2*numBasis]; // z
      }
    } else { // we set everything to zero
      new_coefficients[b] = 0;
      new_powers[b] = 0;
      new_powers[b+numBasis] = 0;
      if (ndims == 3) {
	new_powers[b+2*numBasis] = 0;
      }
    }
  }
}

// dz modifies new_coefficients and new_powers such that it is a first derivative of z
void Polynomial::dz(vector<double>& new_coefficients, vector<double>& new_powers){
  for (int b=0;b<numBasis;b++) {
    if (new_coefficients[b]!=0 && new_powers[b+2*numBasis] > 0){
      new_coefficients[b] = new_coefficients[b]*new_powers[b+2*numBasis];
      new_powers[b] = new_powers[b];                         // x
      new_powers[b+numBasis] = new_powers[b+numBasis];       // y
      new_powers[b+2*numBasis] = new_powers[b+2*numBasis]-1; // z changes
    } else { // we set everything to zero
      new_coefficients[b] = 0;
      new_powers[b] = 0;
      new_powers[b+numBasis] = 0;
      new_powers[b+2*numBasis] = 0;
    }
  }
}

// Simple numerical evaluation of the polynomial after dx() or dy() is called. 
// There's another on above but I wrote my own so that it's consistent with my notation 
double Polynomial::evaluate(vector<double>& new_coefficients, vector<double>& new_powers, const vector<double>& xs){
  int N = xs.size();
  int total = 0;
  int mult;
  for (int b=0;b<numBasis;b++) {                      // for each basis function x,y,x^2,y^2....
    mult = new_coefficients[b];                       // where xs[0] = x, xs[1] = y....
    for (int i=0;i<N;i++) {                           // for each variable (aka dimension) x, y, z
      mult *=  pow(xs[i],new_powers[b+i*numBasis]);   // x^powers[b]*y^powers[b+numBasis] ...
    }
    total = total + mult;                             // then sum that multiple from each basis functions
  }
  return total;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Below are just combinations of dx and dy such that dxy = dx(dy(...)) etc
// There are some repeats (eg: derivXY = derivYX) but they make thing easier for the user
// and might be useful for tests
////////////////////////////////////////////////////////////////////////////////////////////////////
 
double Polynomial::derivX(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dx(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivY(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dy(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
} 

double Polynomial::derivZ(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dz(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivXX(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dx(new_coefficients,new_powers);
  dx(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivYY(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dy(new_coefficients,new_powers);
  dy(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}


double Polynomial::derivXY(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dx(new_coefficients,new_powers);
  dy(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivYX(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dy(new_coefficients,new_powers);
  dx(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivZZ(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dz(new_coefficients,new_powers);
  dz(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivXZ(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dx(new_coefficients,new_powers);
  dz(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivYZ(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dy(new_coefficients,new_powers);
  dz(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivZX(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dz(new_coefficients,new_powers);
  dx(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
}

double Polynomial::derivZY(const vector<double>& xs) {
  vector<double> new_coefficients;
  vector<double> new_powers;
  int value;
  
  new_coefficients = coefficients;
  new_powers = powers;
  dz(new_coefficients,new_powers);
  dy(new_coefficients,new_powers);
  value = evaluate(new_coefficients,new_powers,xs);
  return value;
  }
