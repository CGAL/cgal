// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

#include <CGAL/Subdivision_surfaces_3.h>

#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Cartesian<double>            Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

using namespace std;
using namespace CGAL;

int main(int argc, char **argv) {
  if (argc != 3) { 
    cout << "Usage: Loop_subdivision filename d" << endl; 
    cout << "       filename: the input mash (.off)" << endl; 
    cout << "       d: the depth of the subdivision (0 < d < 10)" << endl; 
    exit(1);
  }

  ifstream in(argv[1]);
  int d = argv[2][0] - '0';

  Polyhedron P;
  in >> P; // read the .off

  Subdivision_surfaces_3<Polyhedron>::Loop_subdivision(P,d);

  cout << P; // write the .off
  
  return 0;
}
