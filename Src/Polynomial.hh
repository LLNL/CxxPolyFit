#ifndef POLYNOMIAL_HH
#define POLYNOMIAL_HH

#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

class Polynomial
{
  public:

  /** This constructor takes a set of points and does a least-squares fit to make a polynomial
   *  xs are the function inputs, so it's a [points x dims] vector
   *  ys are the answers to those inputs, so its a [points] vector
   *  in_order is the order of the polynomial to fit.
   **/
  Polynomial(const vector<double>& xs, const vector<double>& ys, int in_dims, int in_order);

  /* Constructor does NOT take pointer possession, it makes it's own copy. */
  //  Polynomial(const vector<double>& in_powers, const vector<double>& in_coefficients, double in_dims, double in_order);
  Polynomial(const vector<double>& in_coefficients, double in_dims, double in_order);
  Polynomial(const Polynomial& rhs);
  virtual ~Polynomial() { }
  
  static int getNumBasis(int in_dims, int in_order); 
  static void getBasisPowers(int in_dims, int in_order, vector<double>& out_powers);
  /** Evaluates each basis function at X.  Eval runs this to get the basis function values
   * and then multiples in the coefficients*/
  void basisEvals(const vector<double>& xs, vector<double>& outVals);
  
  /** Evaluates each basis derivitive at X.  grad runs this to get the basis function values
   * and then multiples in the coefficients
   * Out Grads has the same shape as powers*/
  void basisGradients(const vector<double>& xs, vector<double>& outGrads);

  /**Calcualte the 1st derivative of every dimension */
  void gradAllDims(const vector<double>& xs, vector<double>& outGrads);

  /**Get the 1st derivative of just the first dimension*/
  double grad1(const vector<double>& xs);
  /**Get the 1st derivative of just the second dimension*/
  double grad2(const vector<double>& xs);
  /**Get the 1st derivative of just the third dimension*/
  double grad3(const vector<double>& xs);

  //////////////////////////////////////////////////////////////////////////////////////////////////

  void dx(vector<double>& new_coefficients, vector<double>& new_powers);
  void dy(vector<double>& new_coefficients, vector<double>& new_powers);
  void dz(vector<double>& new_coefficients, vector<double>& new_powers);

  double evaluate(vector<double>& new_coefficients, vector<double>& new_powers, const vector<double>& xs);

  double derivX(const vector<double>& xs);
  double derivY(const vector<double>& xs);
  double derivZ(const vector<double>& xs);

  double derivXX(const vector<double>& xs);
  double derivYY(const vector<double>& xs);
  double derivZZ(const vector<double>& xs);

  double derivXY(const vector<double>& xs);
  double derivXZ(const vector<double>& xs);

  double derivYX(const vector<double>& xs);
  double derivYZ(const vector<double>& xs);

  double derivZX(const vector<double>& xs);
  double derivZY(const vector<double>& xs);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  /** Evaluate the function at the set of inputs */
  double eval(const vector<double>& xs);

  /** Returns a string giving all the coefficients and their basis functions.  
   *  Helpful for debugging.  Limited to 3 dimensions (like everything.) 
   **/
  string termsToString();

private:

  /** The number of dimensions in the polynomial.  Max is 3.  The is also the number of
   *  rows in the powers and coefficients arrays */
  int dims; 

  /** The order of the polynomial.  This should be less than the number of input points, 
   *  but otherwise can be quite large.  This indirectly determines the number of columns 
   *  in the powers and coefficients arrays.
   */
  int order;

  /** The number of columns in the powers and coefficients arrays.  
   *   Determined from dims and order. 
   */
  int numBasis;

  /** The powers of the the basis functions.  Each basis function is an entry in this array.
   *  pows has rows = dimension, column = powers of each basis so
   *  1 + x + y + x^2 + y^2 + xy =
   *  1 + x^1 * y^0 + y^1 * x^0 + x^2 * y^0 + x^1 * y^1 + y^2 * x^0 =
   *
   *  [ 0, 1, 0, 2, 1, 0
   *    0, 0, 1, 0, 1, 2]*pow(x[i],power[b])*
   */
  vector<double> powers;

  /** coefficients is a 1D array of the coefficients of each basis function.
   * Each basis function has 1 coefficent. The length is numcols.
   */
  vector<double> coefficients;
  

};

#endif
