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
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Polyhedral Surface File IO: Scanner and Writer
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

// This is the test file for the new design. Skip new design test for 
// disabled compilers.
#ifndef CGAL_USE_POLYHEDRON_DESIGN_ONE

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/print_OFF.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/IO/print_inventor.h>
#include <CGAL/IO/print_VRML_1.h>
#include <CGAL/IO/print_VRML_2.h>

#include <cstddef>
#include <iostream>
#include <fstream>
#include <strstream>
#include <cstring>

using namespace CGAL;

char* triangle = "OFF\n"
                 "3 1 0\n"
                 "0 0 0\n"
                 "1 0 0\n"
                 "0 1 0\n"
                 "3  0 1 2\n";
char* tetra =    "OFF\n"
                 "4 4 0\n"
                 "0 0 0.707107\n"
                 "1 1 0.707107\n"
                 "0 1 0\n"
                 "1 0 0\n"
                 "3  1 3 0\n"
                 "3  2 1 0\n"
                 "3  3 2 0\n"
                 "3  2 3 1\n";


void test_file_IO_OFF() {
    typedef Cartesian<double>     Kernel;
    typedef Kernel::Point_3       Point;
    typedef Kernel::Plane_3       Plane;
    typedef Polyhedron_3<Kernel>  Polyhedron;
    {
        Polyhedron P;
        std::istrstream in( triangle, CGAL_CLIB_STD::strlen( triangle));
        in >> P;    /* 'in' is the stream where the object is read from. */
        CGAL_assertion( in);
        CGAL_assertion( P.is_triangle( P.halfedges_begin()));
        char* buffer = new char[100000];
        std::ostrstream out( buffer, 100000);
        out << P << '\0';
        std::istrstream bufin( buffer, CGAL_CLIB_STD::strlen(buffer));
        P = Polyhedron();
        scan_OFF( bufin, P, true);
        CGAL_assertion( bufin);
        CGAL_assertion( P.is_triangle( P.halfedges_begin()));

        std::ostrstream out_new( buffer, 100000);
        print_OFF( out_new, P, true);
        out_new << '\0';
        std::istrstream bufin_new( buffer, 100000);
        P = Polyhedron();
        bufin_new >> P;
        CGAL_assertion( bufin_new);
        CGAL_assertion( P.is_triangle( P.halfedges_begin()));
        {
            std::ofstream out2( "triangle_binary.off");
            print_OFF( out2, P, true);
        }
        std::ifstream filein( "triangle_binary.off");
        P = Polyhedron();
        filein >> P;
        CGAL_assertion( filein);
        CGAL_assertion( P.is_triangle( P.halfedges_begin()));
        delete[] buffer;
    }{
        Polyhedron P;
        std::istrstream in( tetra, CGAL_CLIB_STD::strlen( tetra));
        in >> P;    /* 'in' is the stream where the object is read from. */
        CGAL_assertion( P.is_tetrahedron( P.halfedges_begin()));
        char* buffer = new char[100000];
        std::ostrstream out( buffer, 100000);
        out << P << '\0';
        std::istrstream bufin( buffer, CGAL_CLIB_STD::strlen(buffer));
        P = Polyhedron();
        bufin >> P;
        CGAL_assertion( P.is_tetrahedron( P.halfedges_begin()));
        {
            std::ofstream out2( "tetra_binary.off");
            print_OFF( out2, P, true);
        }
        std::ifstream filein( "tetra_binary.off");
        P = Polyhedron();
        filein >> P;
        CGAL_assertion( P.is_tetrahedron( P.halfedges_begin()));
        delete[] buffer;
    }
}
void test_file_IO_wavefront() {}
void test_file_IO_inventor() {}
void test_file_IO_VRML_2() {}

#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //

int main(){
// This is the test file for the new design. Skip new design test for 
// disabled compilers.
#ifndef CGAL_USE_POLYHEDRON_DESIGN_ONE
    test_file_IO_OFF();
    test_file_IO_wavefront();
    test_file_IO_inventor();
    test_file_IO_VRML_2();
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
    return 0;
}
// EOF //
