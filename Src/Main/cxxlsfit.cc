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


#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fstream>
#include <string>
#include <cerrno>

#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>
#include <errno.h>
#include "Polynomial.hh"

using std::string;
using std::vector;
using std::endl;
using std::cout;

  // trim from start
  static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  }
  
  // trim from end
  static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  }
  
  // trim from both ends
  static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
  }

  static inline std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      if(item.size() > 0) {
        elems.push_back(item);
      }
    }
    return elems;
  }


  static inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
  }


//Reads a 1D function from a file.  The left column is assumed to be the axis, 2nd column the value.
//Reads the first two columns of a file into the back of 2 seperate vectors.
//If the vector passed in is not empty, the new values will just be pushed on behind the existing ones.
//It expects the first 2 columns to all be doubles.  (leading whitespace is OK.)
//A # as the first character on a line will be skipped.
void read1Dfunc(const string& filename, vector<double>& outvec1, vector<double>& outvec2) {

  std::ifstream in(filename.c_str());
  if(in.fail()) {
    std::cerr << "Opening of file " << filename << " failed." << std::endl;
    exit(1);
  }
  std::string line;
  int linenum = 0;

  while(std::getline(in, line)) {
    if(in.fail()) {
      std::cerr << "Failed to read file: " << filename << " line: " << linenum << std::endl;
      exit(1);
    }
    //Catches bad reads and comment lines
    if(line.size() > 0 && line[0] != '#') {
      std::istringstream iss(line);
      double val1 = 0;
      double val2 = 0;
      iss >> val1;
      iss >> val2;
      if(iss.fail()) {
        std::cout << "Error parsing line number: "<< linenum << std::endl;
        std::cout << "Bad line: "<< line << std::endl;
        exit(1);
      }
      outvec1.push_back(val1);
      outvec2.push_back(val2);
    }
    ++linenum;
  }
  
  in.close();
  return;
}
//Reads the first two columns of a file into the back of 2 seperate vectors.
//If the vector passed in is not empty, the new values will just be pushed on behind the existing ones.
//It expects the first 2 columns to all be doubles.  (leading whitespace is OK.)
//A # as the first character on a line will be skipped.
//
// the skip lines argument allows the reader to skip lines before attempting to read the axes.
//   so far this is mainly useful for ThomasFermi files, which don't use any comment indicators like '#', but
//   rather are headed by some lines to skip
void read2Daxis(const string& filename, vector<double>& outvec1, vector<double>& outvec2, int skiplines, int Xaxis, int Yaxis) {

  std::ifstream in(filename.c_str());
  if(in.fail()) {
    std::cerr << "Opening of file " << filename << " failed." << std::endl;
    exit(1);
  }
  std::string line;
  int linenum = 0;

  //The highest column we need to pull.
  int max_column = std::max(Xaxis, Yaxis);

  while(std::getline(in, line)) {
    if(in.fail()) {
      std::cerr << "Failed to read file: " << filename << " line: " << linenum << std::endl;
      exit(1);
    }
    ++linenum;
    line = trim(line);  //Trim whitespace

    //Skip the skiplines, Catches bad reads and comment lines, skip letter lines
    if(linenum > skiplines && line.size() > 0 && line[0] != '#' && line[0] >= '0' && line[0] <= '9') { 
      vector<double> vals;
      std::istringstream iss(line);
      for(int colnum = 0; colnum <= max_column; ++colnum) {
        double val = 0;
        iss >> val;
        vals.push_back(val);
        if(iss.fail()) {
          std::cout << "Error parsing line number: "<< linenum << std::endl;
          std::cout << "Bad line: "<< line << std::endl;
          exit(1);
        }
      }
      //Skip dups  (We may also want to sort and unique)
      if(std::find(outvec1.begin(), outvec1.end(), vals[Xaxis]) == outvec1.end()) {
        outvec1.push_back(vals[Xaxis]);
      }
      if(std::find(outvec2.begin(), outvec2.end(), vals[Yaxis]) == outvec2.end()) {
        outvec2.push_back(vals[Yaxis]);
      }
    }
  }
    
  in.close();
  return;
}


//Reads a file with 2 axes as the first two columns, then multiple column functions.
//Expects the first non-comment line to be column headers. (strings)
//All columns are expected to be doubles.
//The size of funcs determins the number of functions it expects to read.
//It expects row major order, where the X column (defaults to 0) is the one that varies fastest.
//If the vector passed in is not empty, the new values will just be pushed on behind the existing ones.
//
//Blank lines will be skipped.
//Lines starting with '#' are comments, and are skipped
//Params: [in] filename,
//        [out] vector<double> Xaxis, the Xaxis (fastest varying, row major order) 
//        [out] vector<double> Yaxis, the Yaxis (slower varying, row major order)
//        [out] vector< vector<double>* >& funcs, 2D vectors of values.  new'd in this function, must be deleted.
//        [in]  int skiplines   When reading skip the first n lines.  Useful for files that don't use '#' for comments
//        [in]  int Xaxis       Define the X axis column, defaults to 0
//        [in]  int Yaxis       Define the Y axis column, defaults to 1


void read2Dfuncs(const string& filename, vector<double>& Xaxis, vector<double>& Yaxis, vector<string>& headers, vector< vector<double>* >& funcs, int skiplines, int XaxisCol, int YaxisCol) {

  //Annoyingly, we can't guarantee much about the file format, 
  //by called read2Daxis we should get an accurate representation of what the
  //axes actually are, at the cost of opening and running through the file twice
  read2Daxis(filename, Xaxis, Yaxis, skiplines, XaxisCol, YaxisCol);

  std::ifstream in(filename.c_str());
  if(in.fail()) {
    std::cerr << "Opening of file " << filename << " failed." << std::endl;
    exit(1);
  }
  std::string line;
  
  //This is a little weird.  We need to figure out how many columns are in the file. 
  //(We assume all columns extend for the entire file.)
  //So we get lines until we find one that starts with numerals (skip comments, blanks,
  //and the column headers)
  //We then add a function for each column we find past the first two.  Then we reset the file
  //to the beginning.  We should now know how many functions we have, and can read them.
  size_t linenum = 0;
  int funcnum = 0;
  while(funcnum <= 0) {
    std::string line;
    std::getline(in, line);
    line = trim(line);
    ++linenum;
    if(in.fail()) {
      std::cerr << "Failed to read file: " << filename << " line: " << linenum << std::endl;
      exit(1);
    }
    
    //Catches bad reads and comment lines
    if(linenum > skiplines && line.size() > 0 && line[0] != '#' && line[0] >= '0' && line[0] <= '9') {
      std::istringstream iss(line);
      while(!iss.eof()) {

        double val1 = 0;
        iss >> val1;
        if(funcnum != XaxisCol && funcnum != YaxisCol) {  //If it's not an axis, it must be a func
          funcs.push_back(new vector<double>);
        }
        if(iss.eof())
          break;
        if(iss.fail()) {
          std::cout << "Error parsing line number: "<< linenum << std::endl;
          std::cout << "Bad line: "<< line << std::endl;
          exit(1);
        }
        funcnum++; 

      }  //end while
    }  
  }

  in.seekg(0, in.beg); //Reset to the beginning of the file so we can actually start reading functions
  linenum = 0;

  //Now pull out the headers
  bool foundHeaders = false;
  while(!foundHeaders) {
    std::string line;
    std::getline(in, line);
    line = trim(line);
    ++linenum;
    if(in.fail()) {
      std::cerr << "Failed to read file: " << filename << " line: " << linenum << std::endl;
      exit(1);
    }
    
    //Catches bad reads and comment lines
    if(linenum > skiplines && line.size() > 0 && line[0] != '#') {
      foundHeaders = true;
      split(line, ' ', headers);
      for(size_t ii = 0; ii < headers.size(); ++ii) {
        trim(headers[ii]);
      }
    }
  }

  //Now we should know exactly the length of the file, so just count up the lines
  for(size_t yy = 0; yy < Yaxis.size(); ++yy) {
    for(size_t xx = 0; xx < Xaxis.size(); ++xx) {
      std::string line;
      bool goodline = false;

      while(!goodline) { //Skip any blank lines or comments and don't count them

        std::getline(in, line);
        line = trim(line);
        ++linenum;
        if(in.fail()) {
          std::cerr << "Failed to read file: " << filename << " line: " << linenum << std::endl;
          exit(1);
        }

        //Catches bad reads and comment lines
        if(line.size() > 0 && line[0] != '#') {
          goodline = true;
        }
      }

      std::istringstream iss(line);
      double Xval = 0;
      double Yval = 0;
      funcnum = 0;
      int colnum = 0;
      for(size_t ii = 0; ii < funcs.size()+2; ++ii) { //#funcs + #axes (2)
        
        double val = 0;
        iss >> val;
        if(iss.fail()) {
          std::cerr << "Error parsing line number: "<< linenum << std::endl;
          std::cerr << "Bad line: "<< line << std::endl;
          exit(1);
        }

        if(colnum == XaxisCol) {
          Xval = val;
        }else if(colnum == YaxisCol) {
          Yval = val;
        } else {
          funcs[funcnum]->push_back(val);
          ++funcnum;
        }
        ++colnum;
      }

      if(Xval != Xaxis[xx]) {
        std::cerr << "ERROR on line: " << linenum << std::endl;
        std::cerr << "Found mismatched Xvalues.  For index " << xx << std::endl;
        std::cerr << "expected " << Xaxis[xx] << " got: " << Xval << std::endl;
        exit(1);
      }
      if(Yval != Yaxis[yy]) {
        std::cerr << "ERROR on line: " << linenum << std::endl;
        std::cerr << "Found mismatched Yvalues.  For index " << yy << std::endl;
        std::cerr << "expected " << Yaxis[yy] << " got: " << Yval << std::endl;
        exit(1);
      }    
    }
  }    
  in.close();
  return;
}


int main(int argc, char** argv) {

    if (argc < 4 || argc > 4) {
      cout << "cxxlsfit: Does a simple polynomial fit of a set of either 1D or 2D points" << endl;
      cout << "Error:  Incorrect number of arguments." << endl;
      cout << "Usage : cxxlsfit takes 3 arguments:" << endl;
      cout << "1. The number of dimensions in the input file" << endl;
      cout << "2. The order of the polynomial to fit." << endl;
      cout << "3. The input file, in column format." << endl;
      cout << endl;
      cout << "EXAMPLE: meos  1 2 input.dat" << endl;
      exit(-1);
    }

    int dims = atoi(argv[1]);
    int order = atoi(argv[2]);
    string infilename(argv[3]);

    if(dims > 2 || dims < 1) {
      cout << "This program can only fit 1 or 2D functions." << endl;
      exit(-1);
    }

    if(order < 1) {
      cout << "Order must be a positive non-zero integer" << endl;
      exit(-1);
    }


    Polynomial *poly = 0;

    if(dims == 1) {
      vector<double> xs;
      vector<double> ys;
      read1Dfunc(infilename, xs, ys);
      poly = new Polynomial(xs, ys, dims, order);
    } else {
      vector<double> xs1; //First dimension
      vector<double> xs2; //Second dimension
      vector<string> headers;
      vector< vector<double>* > data;

      read2Dfuncs(infilename, xs1, xs2, headers, data, 0, 0, 1);

      vector<double> xs; //Combined xs1 xs2 dimensions
      xs.resize(xs1.size() + xs2.size());
      for(int ii = 0; ii < xs.size(); ++ii) {
        xs[ii*2] = xs1[ii];
        xs[(ii*2)+1] = xs2[ii];
      }
      
      vector<double>& func = *data[0];
      poly = new Polynomial(xs, func, dims, order);
    }

    cout << poly->termsToString();
    delete poly;

}
