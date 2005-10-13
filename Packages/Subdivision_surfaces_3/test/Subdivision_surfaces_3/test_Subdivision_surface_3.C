// ============================================================================
//
// Copyright (c) 2005 Le-Jeng Shiue
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
// file          : test/Subdivision_surface_3/test_Subdivision_surface_3.C
// package       : Subdivision_surface_3
// chapter       : Subdivision Surfaces
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// Test subdivision surfaces
// ============================================================================

#include <CGAL/Subdivision_surfaces_3.h>

#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

using namespace std;

using namespace CGAL;

#define TEST_DEPTH (3)

//#define TESTMESH_GENERAL   "data/??.off"

#define TESTMESH_QUAD      "data/corner.off"
#define TESTMESH_QUAD_OPEN "data/corner_with_hole.off"

#define TESTMESH_TRI       "data/quint_tris.off"
#define TESTMESH_TRI_OPEN  "data/quint_tris.off"

void test_Subdivision_surface_3() {
  typedef CGAL::Cartesian<double>            Kernel;
  typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

  // test Catmull-Clark subdivision on quad mesh
  {
    ifstream mesh(TESTMESH_QUAD);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::CatmullClark_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }

  // test Catmull-Clark subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::CatmullClark_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }


  // test Loop subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::Loop_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }

  // test Loop subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::Loop_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }

  // test Doo-Sabin subdivision on general mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::DooSabin_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }

  // test Sqrt-3 subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_surfaces_3<Polyhedron>::Sqrt3_subdivision(P,TEST_DEPTH);
    CGAL_assertion(P.is_valid());
  }
}


int main() {
    test_Subdivision_surface_3();
    return 0;
}
// EOF //
