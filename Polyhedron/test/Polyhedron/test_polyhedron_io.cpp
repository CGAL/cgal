// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test_polyhedron_io.C
// chapter       : $CGAL_Chapter: Support Library ... $
// package       : $CGAL_Package: Polyhedron_IO 2.11 (04 Feb 2000) $
// source        : polyhedron_io.fw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Polyhedral Surface File IO: Scanner and Writer
// ============================================================================

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cassert>

using namespace CGAL;

typedef Simple_cartesian<double> Kernel;
typedef Polyhedron_3<Kernel>     Polyhedron;

const char* triangle = "OFF\n"
                       "3 1 0\n"
                       "0 0 0\n"
                       "1 0 0\n"
                       "0 1 0\n"
                       "3  0 1 2\n";
const char* tetra =    "OFF\n"
                       "4 4 0\n"
                       "0 0 0.707107\n"
                       "1 1 0.707107\n"
                       "0 1 0\n"
                       "1 0 0\n"
                       "3  1 3 0\n"
                       "3  2 1 0\n"
                       "3  3 2 0\n"
                       "3  2 3 1\n";

const char* empty=    "OFF\n"
                      "0 0 0\n";

void test_file_IO_OFF()
{
  {
    Polyhedron P;
    std::istringstream in( triangle);
    in >> P;    /* 'in' is the stream where the object is read from. */
    assert(! in.fail());
    assert( P.is_triangle( P.halfedges_begin()));
    std::stringstream stream;
    stream << P << '\0';
    P = Polyhedron();
    scan_OFF( stream, P, true);
    assert( ! ! stream); // ! ! to fool VC7 CL1310
    assert( P.is_triangle( P.halfedges_begin()));

    std::stringstream stream_new;
    print_polyhedron_OFF( stream_new, P, true);
    stream_new << '\0';
    P = Polyhedron();
    stream_new >> P;
    assert( ! ! stream_new); // ! ! to fool VC7 CL1310
    assert( P.is_triangle( P.halfedges_begin()));
    {
      std::ofstream out2( "triangle_binary.off");
      print_polyhedron_OFF( out2, P, true);
    }
    std::ifstream filein( "triangle_binary.off");
    P = Polyhedron();
    filein >> P;
    assert(! filein.fail());
    assert( P.is_triangle( P.halfedges_begin()));
  }{
    Polyhedron P;
    std::istringstream in( tetra);
    in >> P;    /* 'in' is the stream where the object is read from. */
    assert( P.is_tetrahedron( P.halfedges_begin()));
    std::stringstream stream;
    stream << P << '\0';
    P = Polyhedron();
    stream >> P;
    assert( P.is_tetrahedron( P.halfedges_begin()));
    {
      std::ofstream out2( "tetra_binary.off");
      print_polyhedron_OFF( out2, P, true);
    }
    std::ifstream filein( "tetra_binary.off");
    P = Polyhedron();
    filein >> P;
    assert( P.is_tetrahedron( P.halfedges_begin()));
  }

  {
    Polyhedron P;
    std::istringstream in( ::empty);
    IO::read_OFF(in, P);
    assert(P.empty());
    assert(in);
  }
}

void test_file_IO_inventor()
{
  std::ofstream o("tmp");
  VRML_1_ostream out(o);
  Polyhedron P;
  out << P;
}

void test_file_IO_VRML_2()
{
  std::ofstream o("tmp");
  VRML_2_ostream out(o);
  Polyhedron P;
  out << P;
}

int main()
{
  test_file_IO_OFF();
  test_file_IO_inventor();
  test_file_IO_VRML_2();

  return 0;
}
