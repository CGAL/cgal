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
// file          : test/Nef_3/include/Nef_3/test_with_homogeneous.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// ============================================================================
#ifndef CGAL_TEST_TEST_NEF_3
#define CGAL_TEST_TEST_NEF_3

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <fstream>

CGAL_BEGIN_NAMESPACE

template<typename Kernel>
class test_Nef_3 {

  typedef CGAL::Nef_polyhedron_3<Kernel>                    Nef_polyhedron;
  typedef typename Nef_polyhedron::Items                    Items;
  typedef typename Nef_polyhedron::Mark                     Mark;
  typedef typename Kernel::RT                               RT;
  typedef typename Kernel::FT                               FT;
  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::Plane_3                          Plane_3;
  typedef typename Kernel::Vector_3                         Vector_3;
  typedef typename Kernel::Segment_3                        Segment_3;
  typedef typename Kernel::Aff_transformation_3             Aff_transformation_3;
  typedef CGAL::Polyhedron_3<Kernel>                        Polyhedron;
  typedef CGAL::SNC_structure<Kernel, Items, Mark>          SNC_structure;

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
  bool cubes_tested;
  bool isolated_edge_tested;

public:
  test_Nef_3() {cubes_tested=false; isolated_edge_tested=false;};

private:

  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;

  struct Shell_explorer {
    bool first;
    const Nef_polyhedron& N;
    SFace_visited_hash& Done;
    Vertex_const_handle v_min;

    Shell_explorer(const Nef_polyhedron& NN, SFace_visited_hash& Vi) 
      : first(true), N(NN), Done(Vi) {}

    void visit(SFace_const_handle h) { 
      Done[h]=true;
    }

    void visit(Vertex_const_handle h) {
      if(first) {
	v_min = h;
	first=false; 
      } else if ( CGAL::lexicographically_xyz_smaller(
		h->point(),v_min->point()) ) 
	v_min = h; 
    }

    void visit(Halfedge_const_handle) {}
    void visit(Halffacet_const_handle) {}
    void visit(SHalfedge_const_handle) {}
    void visit(SHalfloop_const_handle) {}

    Vertex_const_handle& minimal_vertex() { return v_min; }
  };

  static const char* datadir;  

  bool are_files_equal(const char* name1, const char* name2) {
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

  bool does_nef3_equals_file(Nef_polyhedron& N, const char* name) {
    char* fullname = new char[std::strlen(datadir)+std::strlen(name)+1];
    std::strcpy(fullname, datadir);
    std::strcat(fullname, name);
    std::ofstream out("data/temp.nef3");
    out << N;
    bool b = are_files_equal("data/temp.nef3",fullname);
    delete [] fullname;
    return b;
  }

  Nef_polyhedron load_off(const char* name) {
    Polyhedron poly;
    std::ifstream off_file(name);
    CGAL_assertion(off_file != NULL);
    off_file >> poly;
    Nef_polyhedron N(poly);
    return N;
  }

  Nef_polyhedron load_nef3(const char* name) {
    char* fullname = new char[std::strlen(datadir)+std::strlen(name)+1];
    std::strcpy(fullname, datadir);
    std::strcat(fullname, name);
    std::ifstream input(fullname);
    CGAL_assertion_msg(input,"could not open file");
    Nef_polyhedron tmp;
    input >> tmp;
    delete[] fullname;
    return tmp;
  }

  void test_cubes() {
    if(!cubes_tested) {
      Nef_polyhedron N = load_off("data/cube.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));

      N = load_off("data/wrongly_oriented_cube.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));

      N = load_off("data/cube+v.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));
      
      N = load_off("data/cube+vee.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));
      
      N = load_off("data/cube+veeee.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH")); 

      N = load_off("data/cube+vONe.off");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));
      cubes_tested = true;
    }
  }

  void intersect_with_isolated_edge() {
    if(!isolated_edge_tested) {
      Nef_polyhedron N1 = load_nef3("single_edge.nef3.SH");
      N1.is_valid(0,0);
      Nef_polyhedron N2 = load_nef3("intersWithIsolatedEdge.nef3.SH");
      N2.is_valid(0,0);
      Nef_polyhedron N3 = N1.intersection(N2);
      N3.is_valid(0,0);
      CGAL_assertion(does_nef3_equals_file(N3,"intersWithIsolatedEdgeRef1.nef3.SH"));
      N3 = N1.symmetric_difference(N2);
      N3.is_valid(0,0);
      CGAL_assertion(does_nef3_equals_file(N3,"intersWithIsolatedEdgeRef2.nef3.SH"));
      isolated_edge_tested = true;
    }
  }
  
  void loadSave() {

    test_cubes();
    return;
    Polyhedron P;
    Nef_polyhedron N = load_off("data/cube.off");
    N.convert_to_Polyhedron(P);
    std::ofstream out("data/temp.off");
    out << P;
    N = load_off("data/temp.off");
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));

    N = load_nef3("topology.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"topology.nef3.SH"));

    if(Infi_box::extended_kernel()) {
      N = load_nef3("topology.nef3.EH");
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"topology.nef3.EH"));
    }
  }

  void newell() {
    Nef_polyhedron N = load_off("data/star.off");
    CGAL_assertion(N.is_valid(0,0));    
    CGAL_assertion(does_nef3_equals_file(N,"newell.nef3.SH"));
  }

  void transformation() {
    Nef_polyhedron C = load_off("data/cube.off");
    Nef_polyhedron N(C);
    
    N.transform(Aff_transformation_3(-1,0,0,
                                      0,1,0,
                                      0,0,1,1));
    N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(4,0,0,1)));
    
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3.SH"));
    CGAL_assertion(N == C);

    double alpha = CGAL_PI * 20 / 180.0;
    int tmpy = static_cast<int>(std::sin( alpha) * 128*256);
    int tmpx = static_cast<int>(std::cos( alpha) * 128*256);
    RT diry = tmpy;
    RT dirx = tmpx;
    RT sin_alpha;
    RT cos_alpha;
    RT w;
    CGAL::rational_rotation_approximation( dirx, diry, 
					   sin_alpha, cos_alpha, w,
					   RT(1), RT( 1000000));
    Aff_transformation_3 rotx20( w, RT(0), RT(0),
				 RT(0), cos_alpha,-sin_alpha,
				 RT(0), sin_alpha, cos_alpha,
				 w);
    N = C;
    N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,1,0,1)));
    N.transform(Aff_transformation_3( CGAL::SCALING, 6, 4));
    N.transform(Aff_transformation_3(0,-1,0,1,0,0,0,0,1,1));
    N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(0,-4,-1,1)));
    N.transform(Aff_transformation_3( CGAL::SCALING, 5, 12));
    N.transform(Aff_transformation_3(0,-1,0,1,0,0,0,0,1,1));
    N.transform(rotx20);
    N.transform(Aff_transformation_3( CGAL::SCALING, 8, 5));
    N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2,0,0,1)));
    N.transform(rotx20.inverse());
    N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(16,21,2,3)));
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(N == C);

    if(Infi_box::extended_kernel()) {    
      N = Nef_polyhedron(Plane_3(1,1,0,0));
      CGAL_assertion(N.is_valid(0,0));      
      Nef_polyhedron N1(Plane_3(1,-2,3,0));
      CGAL_assertion(N1.is_valid(0,0));
      N = N.join(N1);
      N1 = N;
      N.transform(rotx20);
      N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(1,0,0,1)));
      N.transform(Aff_transformation_3( CGAL::SCALING, 2, 1));
      N.transform(rotx20.inverse());
      N.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2,0,0,1)));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(N == N1);
    }
  }

  void construction() {

    Nef_polyhedron N = Nef_polyhedron(Nef_polyhedron::EMPTY);
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"empty.nef3.SH")); 
    N.clear();
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"empty.nef3.SH")); 

    N = Nef_polyhedron(Nef_polyhedron::COMPLETE);
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"complete.nef3.SH")); 
    
     if(Infi_box::extended_kernel()) {     
      N = Nef_polyhedron(Plane_3(3,4,5,0)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(3,4,5,31));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-31.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,29)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-29.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(1,2,3,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(1,2,3,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-23.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(1,2,4,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(1,2,4,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-23.nef3.EH")); 
         
      N = Nef_polyhedron(Plane_3(4,2,1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p4-2-1-17.nef3.EH"));
     
      N = Nef_polyhedron(Plane_3(2,4,1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-4-1-17.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(-4,-2,-1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p-4--2--1-17.nef3.EH")); 

      N = Nef_polyhedron(Plane_3(2,0,2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-0-2-0.nef3.EH")); 
     
      N = Nef_polyhedron(Plane_3(2,2,0,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-2-0-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(-2,0,-2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p-2-0--2-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(3,2,1,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-2-1-0.nef3.EH")); 
     
      N = Nef_polyhedron(Plane_3(2,3,1,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-3-1-0.nef3.EH")); 

      N = Nef_polyhedron(Plane_3(1,-1,5,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1--1-5-0.nef3.EH")); 
      
      N = Nef_polyhedron(Plane_3(0,0,1,10));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-0-1-10.nef3.EH")); 

      N = Nef_polyhedron(Plane_3(1,1,1,1));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-1-1-1.nef3.EH")); 
    }      
  }

  void point_location_SNC() {
    
    Volume_const_handle c;
    Object_handle o;
    Nef_polyhedron N;

    if(Infi_box::extended_kernel()) {
      point_location_SNC_in("point_location.nef3.EH");

      N = load_nef3("point_location.nef3.EH");
      o = N.locate(Point_3(-3,-3,0));
      CGAL_assertion(assign(c,o));
    }

    N = load_nef3("grid.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    o = N.locate(Point_3(2,2,-1));
    CGAL_assertion(assign(c,o));
    o = N.locate(Point_3(2,0,-1));
    CGAL_assertion(assign(c,o));
    
    N = load_nef3("single_vertex.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    o = N.locate(Point_3(4,1,3));
    CGAL_assertion(assign(c,o));
  }

  void point_location_SNC_in(const char* name) {

    Nef_polyhedron N = load_nef3(name);
    CGAL_assertion(N.is_valid(0,0));
    
    Vertex_const_iterator vin;
    Vertex_const_handle vout;
    CGAL_forall_vertices(vin,N) {
      Object_handle o = N.locate(vin->point());
      CGAL_assertion(assign(vout,o));
      CGAL_assertion(vin == vout);
    }
    
    Halfedge_const_iterator ein;
    Halfedge_const_handle eout;
    CGAL_forall_halfedges(ein,N) {
      Vector_3 d(ein->twin()->source()->point() - ein->source()->point());
      Point_3 s(ein->source()->point());
      Object_handle o = N.locate(s+(d/RT(2)));
      CGAL_assertion(assign(eout,o));
      CGAL_assertion(ein == eout || ein == eout->twin());
    }
    
    Halffacet_const_iterator fin;
    Halffacet_const_handle fout;
    CGAL_forall_halffacets(fin,N) {
      Halffacet_cycle_const_iterator fc(fin->facet_cycles_begin());
      SHalfedge_const_handle e;
      CGAL_assertion(fc.is_shalfedge());
      e = SHalfedge_const_handle(fc);
      SHalfedge_around_facet_const_circulator ec(e),ee(e);
      Vertex_const_handle v_min = (ec++)->source()->center_vertex();
      CGAL_For_all(ec,ee) {
	if (CGAL::lexicographically_xyz_smaller(ec->source()->center_vertex()->point(),
						v_min->point())) 
	  v_min = ec->source()->center_vertex(); 
      }
      Vector_3 orth = fin->plane().orthogonal_vector();
      Vector_3 vec(1-CGAL_NTS abs(orth.hx()),
		   1-CGAL_NTS abs(orth.hy()),
		   1-CGAL_NTS abs(orth.hz()));
      vec = vec / RT(2);
      Object_handle o = N.locate(v_min->point()+vec);
      CGAL_assertion(assign(fout,o));
      CGAL_assertion(fin->plane() == fout->plane() || 
			  fin->plane() == fout->twin()->plane());
    }
    
    Volume_const_iterator Cin;
    Volume_const_handle Cout;
    Vector_3 vec(1,1,1);
    vec = vec / RT(10);
    SFace_visited_hash Done(false);
    CGAL_forall_volumes(Cin,N) {
      Shell_explorer SE(N, Done);
      Shell_entry_const_iterator it;
      CGAL_forall_shells_of(it,Cin)
	N.visit_shell_objects(SFace_const_handle(it),SE);
      Point_3 p(SE.minimal_vertex()->point());
      Object_handle o;
      if(Cin == N.volumes_begin()) {
	o = N.locate(p-vec);
      }
      else {
	o = N.locate(p+vec);
      }
      CGAL_assertion(assign(Cout,o));
      CGAL_assertion(Cin == Cout);
    }
  }

  void intersection() {

    if(Infi_box::standard_kernel()) {
      Nef_polyhedron N = load_nef3("star.nef3.SH");
      SNC_intersection is(*N.sncp());
      
      Point_3 p;
    
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(0,0,0)));
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(2,0,0)));
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(3,0,0)));
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(-1,0,0)));
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(1,1,0)));
      CGAL_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(7,25,11)));
      CGAL_assertion(is.does_contain_internally(
			   Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			   Point_3(1,0,0)));

      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,0,0), Point_3(-1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(0,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,1)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,2,0), Point_3(0,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,3,0), Point_3(0,1,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,4,0), Point_3(0,2,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,5,0), Point_3(1,7,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,7,0), Point_3(1,5,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_assertion(is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_assertion(p == Point_3(0,0,0));

      Halffacet_const_iterator hf;
    
      CGAL_assertion(!is.does_contain_internally(
			    N.halffacets_begin(), Point_3(0,0,0)));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,0),Point_3(-31,15,1)),hf,p));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,-1),Point_3(-31,15,0)),hf,p));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,0),Point_3(-30,15,0)),hf,p));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-32,15,-1),Point_3(-32,15,0)),hf,p));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-16,8,-1),Point_3(-16,8,0)),hf,p));
      CGAL_forall_halffacets(hf, N)
	CGAL_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(0,0,-1),Point_3(0,0,1)),hf,p));    

      int i=0;
      CGAL_forall_halffacets(hf, N) {
	bool b = (i == 13 || i == 15);
	CGAL_assertion( b == is.does_intersect_internally(
				  Segment_3(Point_3(-31,15,-1),Point_3(-31,15,1)), 
                                  hf, p));
	if(b) CGAL_assertion(p == Point_3(-31,15,0));
	i++;
      }

      i=0;
      CGAL_forall_halffacets(hf, N) {
	bool b = (i == 14 || i == 16);
	CGAL_assertion( b == is.does_intersect_internally(
				   Segment_3(Point_3(-15,7,-1), Point_3(-15,7,1)), 
                                   hf, p));
	if(b) CGAL_assertion(p == Point_3(-15,7,0));
	i++;
      }
    }
  }

  void point_location_SM() {
    
    if(Infi_box::extended_kernel()) {  
     
      Nef_polyhedron N = load_nef3("marks_of_halfspheres.nef3.EH");
      CGAL_assertion(N.is_valid(0,0));
      Mark lower, upper;
      int i=0;
      Vertex_const_iterator vi = N.vertices_begin();
 
      do {vi++;} while(++i < 8); //  -1 1 0
      CGAL_assertion(vi->point() == Point_3(-1,1,0));
      SM_explorer SME(N.SMexplorer(vi));
      SM_point_locator PL(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 9); //  0 1 -1
      CGAL_assertion(vi->point() == Point_3(0,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 1);
      
      do {vi++;} while(++i < 10); //  1 1 0
      CGAL_assertion(vi->point() == Point_3(1,1,0));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 1);
      
      do {vi++;} while(++i < 11); //  0 1 1
      CGAL_assertion(vi->point() == Point_3(0,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 12); //  0 0 -1
      CGAL_assertion(vi->point() == Point_3(0,0,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 1);
      
      do {vi++;} while(++i < 13); //  0 0 1
      CGAL_assertion(vi->point() == Point_3(0,0,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 14); //  -1 0 0
      CGAL_assertion(vi->point() == Point_3(-1,0,0));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 15); //  1 0 0
      CGAL_assertion(vi->point() == Point_3(1,0,0));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 1);
      
      do {vi++;} while(++i < 17); //  -1 1 1
      CGAL_assertion(vi->point() == Point_3(-1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 18); //  1 1 1
      CGAL_assertion(vi->point() == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 21); //  -1 1 -1
      CGAL_assertion(vi->point() == Point_3(-1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);
      
      do {vi++;} while(++i < 22); //  1 1 -1
      CGAL_assertion(vi->point() == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 1);

      do {vi++;} while(++i < 29); //  -0.5 -0.5 -1
      CGAL_assertion(vi->point() == Point_3(-1,-1,-2,2));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);


      N = load_nef3("rotatedcube.nef3.SH");
      vi = N.vertices_begin();
      i=0;

      CGAL_assertion(vi->point() == Point_3(-57993689, 27367811, 6449, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(-45133849, -45554371, 6449, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(-6436271, 36459971, -52359431, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(-6423569, 36462211, 52359431, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(6423569, -36462211, -52359431, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(6436271, -36459971, 52359431, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(45133849, 45554371, -6449, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 1);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(57993689, -27367811, -6449, 37023709));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);
      

      N = load_nef3("two_edges.nef3.SH");
      vi = N.vertices_begin();
      i=0;

      CGAL_assertion(vi->point() == Point_3(1,-1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);
     

      N = load_nef3("plane-vertex.nef3.SH");
      vi = N.vertices_begin();
      i=0;

      CGAL_assertion(vi->point() == Point_3(-1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(-1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(0,1,0));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);


      N = load_nef3("is_vertex.nef3.SH");
      vi = N.vertices_begin();
      i=0;

      CGAL_assertion(vi->point() == Point_3(-1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(-1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(0,-2379,-8118, 5741));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(0,1,0));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);

      vi++; ++i;
      CGAL_assertion(vi->point() == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 0);
      CGAL_assertion(upper == 0);
    
      N = load_nef3("SMlocateOnSEdge.nef3.SH");
      CGAL_assertion(N.is_valid(0,0));
      vi = N.vertices_begin();
      i=0;

      do {vi++;} while(++i < 7); 
      CGAL_assertion(vi->point() == Point_3(0,0,-1,2));
      SME = N.SMexplorer(vi);
      PL = SM_point_locator(&*vi);
      PL.marks_of_halfspheres(lower,upper,2);
      CGAL_assertion(lower == 1);
      CGAL_assertion(upper == 1);

    }
  }

  void simplification_SNC() {

    test_cubes();
    
    Nef_polyhedron N = load_nef3("simplifySface.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    N = N.regularization();
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"simplifySfaceRef.nef3.SH"));    
    
    N = load_nef3("donotsimplify.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    Nef_polyhedron N1 = load_nef3("openSquare.nef3.SH");
    CGAL_assertion(N1.is_valid(0,0));
    N = N.join(N1);
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"donotsimplifyRef.nef3.SH"));

    N = load_nef3("unionfind.nef3.SH");
    CGAL_assertion(N.is_valid(0,0));
    N = N.regularization();
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"unionfindRef.nef3.SH"));

    N = load_nef3("nestedFCafterSimpl.nef3.SH");
    N = N.closure();
    CGAL_assertion(does_nef3_equals_file(N,"nestedFCafterSimplRef.nef3.SH"));
  }

  void simplification_SM() {

    Nef_polyhedron C = load_off("data/cube.off");
    Nef_polyhedron N  = C;
    Nef_polyhedron N1 = N;
    Nef_polyhedron N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,1,0,1)));
    CGAL_assertion(N2.is_valid(0,0));
    N1=N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2,-1,4,2)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+plane.nef3.SH"));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,0,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,2,1)));
    N1 = N.join(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+line.nef3.SH"));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,2,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,0,1)));
    N2 = N.symmetric_difference(N1);
    CGAL_assertion(N2.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N2,"cube+vertex1.nef3.SH"));

    N  = C;
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(0,0,-1,1)));
    N2 = N.symmetric_difference(N1);
    CGAL_assertion(N2.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N2,"cube+vertex2.nef3.SH"));
  }

  void synthesis() {

    if(Infi_box::extended_kernel()) {
      Nef_polyhedron N = Nef_polyhedron(Nef_polyhedron::COMPLETE);  
      CGAL_assertion(N.is_valid(0,0));               
      Nef_polyhedron N1 = Nef_polyhedron(Plane_3(1,0,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(-1,0,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,1,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,-1,0,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,0,1,-1));
      N = N.intersection(N1);
      CGAL_assertion(N.is_valid(0,0));
      N1 = Nef_polyhedron(Plane_3(0,0,-1,-1));
      N = N.intersection(N1);      
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube_created_from_halfspaces.nef3.SH"));
    }

    Nef_polyhedron N,N2,P,R,S,T;
    N = load_off("data/centered_cube.off");
    CGAL_assertion(N.is_valid(0,0));
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::SCALING, 4, 1));
    CGAL_assertion(N2.is_valid(0,0));
    P = N;
    P.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2,-2,0,1)));
    CGAL_assertion(P.is_valid(0,0));
    R = N;
    R.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-2, 2,0,1)));
    CGAL_assertion(R.is_valid(0,0));
    S = N;
    S.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3( 2,-2,0,1)));
    CGAL_assertion(S.is_valid(0,0));
    T = N;
    T.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3( 2, 2,0,1)));
    CGAL_assertion(T.is_valid(0,0));
    N2 = N2.difference(P);
    CGAL_assertion(N2.is_valid(0,0));
    N2 = N2.difference(T);
    CGAL_assertion(N2.is_valid(0,0));
    N2 = N2.difference(R);
    CGAL_assertion(N2.is_valid(0,0));
    N2 = N2.difference(S);
    CGAL_assertion(N2.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N2,"synthesis.nef3.SH"));

    intersect_with_isolated_edge();
  }

  void unary_operations() {
    
    Nef_polyhedron X = load_nef3("pl_update_test.nef3.SH");
    Nef_polyhedron Y = load_off("data/centered_cube.off");
    X=X.closure();
    X.difference(Y);
    CGAL_assertion(X.is_valid());

    if(Infi_box::extended_kernel()) {
      Nef_polyhedron T = load_nef3("topology.nef3.EH");
      Nef_polyhedron N = T;
      N = N.boundary();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"boundary.nef3.EH"));
    
      N = T;
      N = N.interior();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"interior.nef3.EH"));
      
      N = T;
      N = N.closure();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"closure.nef3.EH"));
      
      N = T;
      N = N.regularization();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"regularization.nef3.EH"));
      
      N = T; 
      N = N.complement();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"complement.nef3.EH"));
    }
  }
    
  void mark_evaluation() {
    
    if(Infi_box::standard_kernel()) {    
      Nef_polyhedron R;
      Nef_polyhedron P = load_nef3("mark_eval.nef3.SH");
      Nef_polyhedron N = load_nef3("mark_eval2.nef3.SH");
      Nef_polyhedron P2 = P;
      Nef_polyhedron N2 = N;
      P2.transform(Aff_transformation_3( CGAL::SCALING, 2, 1));
      N2.transform(Aff_transformation_3( CGAL::SCALING, 2, 1));
      
      R = P2.intersection(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePPinters.nef3.SH"));
      R = P2.intersection(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePNinters.nef3.SH"));
      R = N2.intersection(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNPinters.nef3.SH"));
      R = N2.intersection(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNNinters.nef3.SH"));
		    
      R = P2.join(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePPjoin.nef3.SH"));
      R = P2.join(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePNjoin.nef3.SH"));
      R = N2.join(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNPjoin.nef3.SH"));
      R = N2.join(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNNjoin.nef3.SH"));

      R = P2.difference(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePPdiff.nef3.SH"));
      R = P2.difference(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePNdiff.nef3.SH"));
      R = N2.difference(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNPdiff.nef3.SH"));
      R = N2.difference(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNNdiff.nef3.SH"));
		    
      R = P2.symmetric_difference(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePPsymdiff.nef3.SH"));
      R = P2.symmetric_difference(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "mePNsymdiff.nef3.SH"));
      R = N2.symmetric_difference(P);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNPsymdiff.nef3.SH"));
      R = N2.symmetric_difference(N);
      CGAL_assertion(R.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(R, "meNNsymdiff.nef3.SH"));
    }
  }		     
     
  void bug_scenarios() {
    
    Nef_polyhedron N = load_off("data/octa.off");    
    Nef_polyhedron N1 = N;
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3( 3, 4,5,9)));
    CGAL_assertion(N1.is_valid(0,0));
    N = N.symmetric_difference(N1);
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N, "octa_ref.nef3.SH"));

    N = load_off("data/2_cycles_on_halfsphere.off");
    N1 = load_off("data/2_cycles_on_halfsphere2.off");
    N = N1.difference(N);
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N, "2_cycles_on_halfsphere_ref.nef3.SH"));
  }
    
public:
  void run_test() {
    loadSave();
    newell();
    transformation();
    construction();
    point_location_SNC();
    intersection();
    point_location_SM();
    simplification_SNC();
    simplification_SM();
    synthesis();
    unary_operations();
    mark_evaluation();
    bug_scenarios();
  }

};

template<typename Kernel>
const char* test_Nef_3<Kernel>::datadir="data/";

CGAL_END_NAMESPACE
#endif // CGAL_TEST_NEF_3
