// ============================================================================
//
// Copyright (c) 2005-2006, 2017 Le-Jeng Shiue
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
// file          : test/Subdivision_method_3/test_Subdivision_method_3.C
// package       : Subdivision_method_3
// chapter       : Subdivision Method
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Le-Jeng Shiue <Andy.Shiue@gmail.com>
//
// Test subdivision methods
// ============================================================================

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/subdivision_method_3.h>

#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;
using namespace CGAL;

#define TEST_DEPTH (3)

//#define TESTMESH_GENERAL   "data/??.off"

#define TESTMESH_QUAD      "data/corner.off"
#define TESTMESH_QUAD_OPEN "data/corner_with_hole.off"

#define TESTMESH_TRI       "data/quint_tris.off"
#define TESTMESH_TRI_OPEN  "data/nefertiti.off"

void test_Subdivision_surface_3() {
  typedef CGAL::Simple_cartesian<double>     Kernel;
  typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

  // test Catmull-Clark subdivision on quad mesh
  {
    ifstream mesh(TESTMESH_QUAD);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Catmull-Clark subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }


  // test Loop subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Loop subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Doo-Sabin subdivision on general mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::DooSabin_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Sqrt3_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }
}

void test_Subdivision_surface_3_SM() {
  typedef CGAL::Simple_cartesian<double>            Kernel;
  typedef CGAL::Surface_mesh<Kernel::Point_3>       Polyhedron;

  // test Catmull-Clark subdivision on quad mesh
  {
    ifstream mesh(TESTMESH_QUAD);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Catmull-Clark subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }


  // test Loop subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Loop subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Doo-Sabin subdivision on general mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::DooSabin_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P, true));
  }

  // test Doo-Sabin subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::DooSabin_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Sqrt3_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Sqrt3_subdivision(P);
    assert(CGAL::is_valid_polygon_mesh(P));
  }
}

void test_Subdivision_surface_3_SM_NP() {
  typedef CGAL::Simple_cartesian<double>            Kernel;
  typedef CGAL::Surface_mesh<Kernel::Point_3>       Polyhedron;

  // test Catmull-Clark subdivision on quad mesh
  {
    ifstream mesh(TESTMESH_QUAD);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                                   .number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Catmull-Clark subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::CatmullClark_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                                   .number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }


  // test Loop subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                           .number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Loop subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Loop_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                           .number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Doo-Sabin subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::DooSabin_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                               .number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Doo-Sabin subdivision on 'opened' quad mesh
  {
    ifstream mesh(TESTMESH_QUAD_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::DooSabin_subdivision(P,Subdivision_method_3::parameters::number_of_iterations(TEST_DEPTH));
    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on tri mesh
  {
    ifstream mesh(TESTMESH_TRI);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Sqrt3_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                            .number_of_iterations(TEST_DEPTH));

    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on 'opened' tri mesh
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    Subdivision_method_3::Sqrt3_subdivision(P,Subdivision_method_3::parameters::vertex_point_map(get(vertex_point, P))
                                            .number_of_iterations(TEST_DEPTH));

    std::ofstream out("out_0.off");
    out << P;

    assert(CGAL::is_valid_polygon_mesh(P));
  }

  // test Sqrt-3 subdivision on 'opened' tri mesh & with external property map
  {
    ifstream mesh(TESTMESH_TRI_OPEN);

    Polyhedron P;
    mesh >> P;

    typedef Kernel::Point_3                                               Point;
    typedef Kernel::Vector_3                                              Vector;

    typedef boost::graph_traits<Polyhedron>::vertex_descriptor            vertex_descriptor;
    typedef boost::unordered_map<vertex_descriptor, Kernel::Point_3>      Point_pmap;
    typedef boost::associative_property_map<Point_pmap>                   APM;
    typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type   VPM;

    Point_pmap um;
    APM apm(um);
    VPM vpm = get(vertex_point, P);

    // some arbitrary new coordinates (different from the internal vpm)
    for(vertex_descriptor vd : vertices(P)) {
      boost::property_traits<VPM>::reference pt = get(vpm, vd);
      Vector v = pt - Point(0., 0., -3.);
      put(apm, vd, pt + 0.5*v);
    }

    Subdivision_method_3::Sqrt3_subdivision(P,
                         Subdivision_method_3::parameters::vertex_point_map(apm)
                         .number_of_iterations(TEST_DEPTH));

    assert(CGAL::is_valid_polygon_mesh(P));
  }
}

int main() {
    test_Subdivision_surface_3();
    test_Subdivision_surface_3_SM();
    test_Subdivision_surface_3_SM_NP();
    std::cerr << "Done" << std::endl;
    return 0;
}
// EOF //
