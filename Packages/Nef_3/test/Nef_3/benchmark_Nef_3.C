// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
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
// file          : test/Nef_3/nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Demo program maintaining a stack of Nef polyhedra in the space and
// a manipulation language for stack ops, file loading and saving, etc.
// ============================================================================

#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>
//#include <CGAL/Extended_homogeneous_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Timer.h>

#include <CGAL/Nef_3/binop_intersection_tests.h>
#include <CGAL/Unique_hash_map.h>


template<typename Kernel>
class benchmark {

  typedef typename Kernel::RT                               RT;
  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::Plane_3                          Plane_3;
  typedef typename Kernel::Vector_3                         Vector_3;
  typedef typename Kernel::Segment_3                        Segment_3;
  typedef typename Kernel::Aff_transformation_3             Aff_transformation_3;
  typedef CGAL::Polyhedron_3<Kernel>                        Polyhedron;
  typedef CGAL::SNC_items<Kernel, bool>                     SNC_items;
  typedef CGAL::SNC_structure<SNC_items>                    SNC_structure;
  typedef typename SNC_structure::Vertex_const_iterator     Vertex_const_iterator;
  typedef typename SNC_structure::Vertex_handle             Vertex_handle;
  typedef typename SNC_structure::Vertex_const_handle       Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_iterator   Halfedge_const_iterator;
  typedef typename SNC_structure::Halfedge_handle           Halfedge_handle;
  typedef typename SNC_structure::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_iterator  Halffacet_const_iterator;
  typedef typename SNC_structure::Halffacet_handle          Halffacet_handle;
  typedef typename SNC_structure::Halffacet_const_handle    Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator
                                  Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Volume_handle             Volume_handle;
  typedef typename SNC_structure::Volume_const_iterator     Volume_const_iterator;
  typedef typename SNC_structure::Object_handle             Object_handle;
  typedef typename SNC_structure::SHalfedge_handle          SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle              SFace_handle;
  typedef typename SNC_structure::SFace_const_handle        SFace_const_handle;
  typedef typename SNC_structure::Shell_entry_const_iterator 
                                  Shell_entry_const_iterator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
                                  SHalfedge_around_facet_const_circulator;

  typedef CGAL::Nef_polyhedron_3<SNC_items>              Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_explorer          SNC_explorer;
  typedef typename Nef_polyhedron::SM_explorer           SM_explorer;
  typedef CGAL::SNC_intersection<SNC_structure>          SNC_intersection;
  typedef CGAL::SNC_decorator< SNC_structure >           SNC_decorator;
public:
  benchmark() {};

private:
    Nef_polyhedron& load_nef3( char* filename ) {
      typedef CGAL::Unique_hash_map< char*, Nef_polyhedron* > hash_map;
      static hash_map nef_map;
      Nef_polyhedron *p;

      if( !nef_map[ filename ] ) {
          std::cout << "loading " << filename << std::endl;
          p = new Nef_polyhedron( filename );
          nef_map[ filename ] =  p;
      } else
          p = nef_map[ filename ];
      return *p;
    }

    void test_intersection( char* filename_A, char* filename_B ) {
      CGAL::Timer t;

      CGAL::binop_intersection_tests_allpairs< SNC_decorator, typename Nef_polyhedron::AND, true > impl_allpairs;
      CGAL::binop_intersection_tests_segment_tree< SNC_decorator, typename Nef_polyhedron::AND, true > impl_segment_tree;

      Nef_polyhedron& A = load_nef3( filename_A );
      Nef_polyhedron& B = load_nef3( filename_B );

      std::cout << "intersecting " << filename_A << " against " << filename_B << " using all pairs ... " << std::flush;
      t.start();
      A.intersection( B, impl_allpairs );
      t.stop();
      std::cout << "overall time spent: " << t.time() << std::endl;
      t.reset();

      std::cout << "intersecting " << filename_A << " against " << filename_B << " using segment tree ... " << std::flush;
      t.start();
      A.intersection( B, impl_segment_tree );
      t.stop();
      std::cout << "overall time spent: " << t.time() << std::endl;

      //std::ofstream out("data/temp.nef3");
      //B.dump(out);
      //bool b = are_files_equal("data/temp.nef3",filename_B );
      //assert( b );
  }

public:
    void run_test(char* suffix) {
        /*char* sirp1 = "../../demo/Nef_3/sirpinski1.nef3";
        char* sirp2 = "../../demo/Nef_3/sirpinski2.nef3";
        char* sirp3 = "../../demo/Nef_3/sirpinski3.nef3";*/
        char* cube        = "data/cube.nef3.SH";
        char* seven_cubes = "data/seven-cubes.nef3.SH";
        char* many_cubes  = "data/49-cubes.nef3.SH";

        /*test_intersection( sirp1, sirp1 );
        test_intersection( sirp2, sirp2 );
        test_intersection( sirp3, sirp3 );
        test_intersection( sirp1, seven_cubes );
        test_intersection( sirp2, seven_cubes );
        test_intersection( sirp3, seven_cubes );*/
        test_intersection( cube, cube );
        test_intersection( cube, seven_cubes );
        test_intersection( seven_cubes, seven_cubes );
        test_intersection( many_cubes, seven_cubes );
        test_intersection( many_cubes, many_cubes );
    }

    void run_single_test( char* filename ) {
        test_intersection( filename, filename );
    }

};


int main( int argc, char** argv) {

  typedef CGAL::Gmpz                         NT;
  typedef CGAL::Simple_homogeneous<NT>       SH_Kernel;
  //typedef CGAL::Extended_homogeneous_3<NT>   EH_Kernel;
  
  CGAL::Timer t;
  t.start();

  benchmark<SH_Kernel> b;
  if( argc == 2 )
      b.run_single_test( argv[1] );
  else
      b.run_test(".SH");

  t.stop();
  std::cout << "Time " << t.time() << std::endl;
}

