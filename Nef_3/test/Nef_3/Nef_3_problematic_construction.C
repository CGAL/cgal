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
// revision      : $Revision: 1.39 $
// revision_date : $Date: 2006/01/20 15:18:36 $
//
// author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// ============================================================================
#include <CGAL/Nef_2/Nef_polynomial.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_items.h>
#include <CGAL/Timer.h>
#include <fstream>

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/OFF_to_nef_3.h>
typedef CGAL::Gmpz NT;
typedef CGAL::Gmpq FNT;
typedef CGAL::Quotient<NT> FNT2;

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
typedef leda_integer LNT;
typedef leda_rational LFNT;
typedef CGAL::Quotient<LNT> LFNT2;
#endif

template<typename Kernel>
class test {

  typedef CGAL::Nef_polyhedron_3<Kernel>                    Nef_polyhedron;
  typedef typename Nef_polyhedron::Mark                     Mark;
  typedef typename Kernel::RT                               RT;
  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::Plane_3                          Plane_3;
  typedef typename Kernel::Vector_3                         Vector_3;
  typedef typename Kernel::Segment_3                        Segment_3;
  typedef typename Kernel::Aff_transformation_3             Aff_transformation_3;
  typedef CGAL::Polyhedron_3<Kernel>                        Polyhedron;
  typedef CGAL::SNC_structure<Kernel,CGAL::SNC_items, Mark> SNC_structure;
  typedef typename SNC_structure::Vertex_const_iterator     Vertex_const_iterator;
  typedef typename SNC_structure::Vertex_const_handle       Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_iterator   Halfedge_const_iterator;
  typedef typename SNC_structure::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_iterator  Halffacet_const_iterator;
  typedef typename SNC_structure::Halffacet_const_handle    Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator
                                  Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Volume_const_iterator     Volume_const_iterator;
  typedef typename SNC_structure::Volume_const_handle       Volume_const_handle;
  typedef typename SNC_structure::Object_handle             Object_handle;
  typedef typename SNC_structure::SHalfedge_const_handle    SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfloop_const_handle    SHalfloop_const_handle;
  typedef typename SNC_structure::SFace_const_handle        SFace_const_handle;
  typedef typename SNC_structure::Shell_entry_const_iterator 
                                  Shell_entry_const_iterator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator
                                  SHalfedge_around_facet_const_circulator;

  typedef typename SNC_structure::Infi_box                  Infi_box;

  typedef typename Nef_polyhedron::SM_explorer           SM_explorer;
  typedef CGAL::SNC_intersection<SNC_structure>          SNC_intersection;
  typedef typename SNC_structure::Sphere_map             Sphere_map;
  typedef CGAL::SM_const_decorator<Sphere_map>           SM_decorator;
  typedef CGAL::SM_point_locator<SM_decorator>           SM_point_locator;

private:
  static const char* datadir;  

  bool are_files_equal(char* name1, char* name2) {
    std::ifstream in1(name1);
    std::ifstream in2(name2);
    std::string s1;
    std::string s2;
    while(in1) {
      in1 >> s1;
      in2 >> s2;
      if(s1 != s2) {
        std::cerr << s1 << std::endl;
	std::cerr << s2 << std::endl;
	return false;
      }
    }
    if(in2)
      return false;
    return true;
  }

  bool does_nef3_equals_file(Nef_polyhedron& N, char* name, char* suffix) {
    char* fullname = new char[strlen(datadir)+strlen(name)+strlen(suffix)+1];
    strcpy(fullname, datadir);
    strcat(fullname, name);
    strcat(fullname, suffix);
    std::ofstream out("data/temp.nef3");
    out << N;
    bool b = are_files_equal("data/temp.nef3",fullname);
    delete [] fullname;
    return b;
  }

  Nef_polyhedron built_nef_from_off(const char *name) {

     Nef_polyhedron N;

     char* fullname = new char[strlen(datadir)+strlen(name)+1];
     strcpy(fullname, datadir);
     strcat(fullname, name);

     std::ifstream off_file (fullname);
     CGAL_assertion(off_file != NULL);

     std::size_t discarded = CGAL::OFF_to_nef_3 (off_file, N, true);
     CGAL_assertion(discarded == 0);
     return N;
  }

public:
  void run_test(bool compare,char* suffix) {
    Nef_polyhedron N = built_nef_from_off( "nine_planes.off");
    if(compare)
      CGAL_assertion(does_nef3_equals_file(N,"nine_planes.nef3",suffix));
  }
};

template<typename Kernel>
const char* test<Kernel>::datadir="data/";

int main() {

  CGAL::Timer t;
  t.start();

  { typedef CGAL::Homogeneous<NT>              H_kernel;
    typedef CGAL::Cartesian<FNT>               C_kernel;
    //    typedef CGAL::Cartesian<FNT2>              Q_kernel;

    test<H_kernel>  test_H;
    test<C_kernel>  test_C;
    // test<Q_kernel>  test_Q;

    test_H.run_test(true,".H");
    test_C.run_test(false,".C");
    // test_Q.run_test();
  }

#ifdef CGAL_USE_LEDA
  { typedef CGAL::Homogeneous<LNT>              LH_kernel;
    typedef CGAL::Cartesian<LFNT>               LC_kernel;
    //    typedef CGAL::Cartesian<LFNT2>              LQ_kernel;
    
    test<LH_kernel>  test_LH;
    test<LC_kernel>  test_LC;
    //    test<LQ_kernel>  test_LQ;
    
    test_LH.run_test(true,".H");
    test_LC.run_test(false,".LC");
    //    test_LQ.run_test();
  }
#endif

  t.stop();
  std::cout << "Time " << t.time() << std::endl;
}

