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

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Timer.h>

#include <fstream>

template<typename Kernel>
class test {

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

public:
  test() {};

private:

  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;

  struct Shell_explorer {
    bool first;
    const SNC_explorer& E;
    SFace_visited_hash& Done;
    Vertex_const_handle v_min;

    Shell_explorer(const SNC_explorer& EE, SFace_visited_hash& Vi) 
      : first(true), E(EE), Done(Vi) {}

    void visit(SFace_const_handle h) { 
      Done[h]=true;
    }

    void visit(Vertex_const_handle h) { 
      if(first) {
	v_min = h;
	first=false; 
      } else if ( CGAL::lexicographically_xyz_smaller(
		E.point(h),E.point(v_min)) ) 
	v_min = h; 
    }

    void visit(Halfedge_const_handle h) {}
    void visit(Halffacet_const_handle h) {}

    Vertex_const_handle& minimal_vertex() { return v_min; }
  };

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
    N.dump(out);
    bool b = are_files_equal("data/temp.nef3",fullname);
    delete [] fullname;
    return b;
  }

  Nef_polyhedron load_off( char* name) {
    Polyhedron poly;
    std::ifstream in(name);
    in >> poly;
    Nef_polyhedron N(poly);
    return N;
  }

  Nef_polyhedron load_nef3(char* name, char* suffix) {
    char* fullname = new char[strlen(datadir)+strlen(name)+strlen(suffix)+1];
    strcpy(fullname, datadir);
    strcat(fullname, name);
    strcat(fullname, suffix);
    Nef_polyhedron tmp(fullname);
    delete[] fullname;
    return tmp;
  }

  void loadSave(char* suffix) {
    Nef_polyhedron N;

    if(suffix[1] == 'E') {
      N = load_nef3("topologyE.nef3",suffix);
      CGAL_nef3_assertion(N.is_valid(0,0));
      CGAL_nef3_assertion(does_nef3_equals_file(N,"topologyE.nef3",suffix));
    }
    
    N = load_off("data/cube.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3",suffix));

    Polyhedron P;
    N.convert_to_Polyhedron(P);
    std::ofstream out("data/temp.off");
    out << P;
    CGAL_nef3_assertion(are_files_equal("data/temp.off","data/cube1.off"));
  }

  void construction(char* suffix) {

    if(suffix[1] == 'E') {
      
      Nef_polyhedron N = Nef_polyhedron(Nef_polyhedron::EMPTY);
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"empty.nef3",suffix)); 

      N = Nef_polyhedron(Nef_polyhedron::COMPLETE);
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"complete.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(3,4,5,0)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-0.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(3,4,5,31));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-4-5-31.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-0.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(0,2,2,29)); 
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-2-2-29.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(1,2,3,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(1,2,3,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-3-23.nef3",suffix)); 
     
      N = Nef_polyhedron(Plane_3(1,2,4,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-0.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(1,2,4,23));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1-2-4-23.nef3",suffix)); 
      
      N = Nef_polyhedron(Plane_3(4,2,1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p4-2-1-17.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(2,4,1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-4-1-17.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(-4,-2,-1,17));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p-4--2--1-17.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(2,0,2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-0-2-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(2,2,0,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-2-0-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(-2,0,-2,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p-2-0--2-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(3,2,1,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p3-2-1-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(2,3,1,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p2-3-1-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(1,-1,5,0));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p1--1-5-0.nef3",suffix)); 

      N = Nef_polyhedron(Plane_3(0,0,1,10));
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"cube+p0-0-1-10.nef3",suffix)); 
    }      
  }

  void point_location_SNC(char* suffix) {
    
    point_location_SNC_in("point_location_2.nef3",suffix);

    Volume_handle c;
    Nef_polyhedron N;
    if(suffix[1] == 'E') {
      N = load_nef3("grid.nef3",suffix);
      CGAL_nef3_assertion(N.is_valid(0,0));
      Object_handle o = N.locate(Point_3(2,2,-1));
      CGAL_nef3_assertion(assign(c,o));
      o = N.locate(Point_3(2,0,-1));
      CGAL_nef3_assertion(assign(c,o));

      N = load_nef3("point_location_2.nef3",suffix);
      o = N.locate(Point_3(-3,-3,0));
      CGAL_nef3_assertion(assign(c,o));
    
      N = load_nef3("single_vertex.nef3",suffix);
      CGAL_nef3_assertion(N.is_valid(0,0));
      o = N.locate(Point_3(4,1,3));
      CGAL_nef3_assertion(assign(c,o));
    }
  }

  void point_location_SNC_in(char* name, char* suffix) {

    if(suffix[1] == 'E') {
      Nef_polyhedron N = load_nef3(name,suffix);
      CGAL_nef3_assertion(N.is_valid(0,0));
      SNC_explorer E(N.SNCexplorer());

      Vertex_const_iterator vin;
      Vertex_handle vout;
      CGAL_nef3_forall_vertices(vin,E) {
	Object_handle o = N.locate(E.point(vin));
	CGAL_nef3_assertion(assign(vout,o));
	CGAL_nef3_assertion(vin == vout);
      }
      
      Halfedge_const_iterator ein;
      Halfedge_handle eout;
      CGAL_nef3_forall_halfedges(ein,E) {
	Vector_3 d(E.point(E.source(E.twin(ein))) - E.point(E.source(ein)));
	Point_3 s(E.point(E.source(ein)));
	Object_handle o = N.locate(s+(d/RT(2)));
	CGAL_nef3_assertion(assign(eout,o));
	CGAL_nef3_assertion(ein == eout || ein == E.twin(eout));
      }

      Halffacet_const_iterator fin;
      Halffacet_handle fout;
      CGAL_nef3_forall_halffacets(fin,E) {
	Halffacet_cycle_const_iterator fc(fin->facet_cycles_begin());
	SHalfedge_handle e;
	CGAL_nef3_assertion(assign(e, fc));
	SHalfedge_around_facet_const_circulator ec(e),ee(e);
	Vertex_const_handle v_min = E.vertex(ec++);
	CGAL_For_all(ec,ee) {
	  if (CGAL::lexicographically_xyz_smaller(E.point(E.vertex(ec)),
						  E.point(v_min))) 
	    v_min = E.vertex(ec); 
	}
	Vector_3 orth = E.plane(fin).orthogonal_vector();
	Vector_3 vec(1-CGAL_NTS abs(orth.hx()),
		     1-CGAL_NTS abs(orth.hy()),
		     1-CGAL_NTS abs(orth.hz()));
	vec = vec / RT(2);
	Object_handle o = N.locate(E.point(v_min)+vec);
	CGAL_nef3_assertion(assign(fout,o));
	CGAL_nef3_assertion(E.plane(fin) == E.plane(fout) || 
			    E.plane(fin) == E.plane(E.twin(fout)));
      }

      Volume_const_iterator Cin;
      Volume_handle Cout;
      Vector_3 vec(1,1,1);
      vec = vec / RT(2);
      SFace_visited_hash Done(false);
      CGAL_nef3_forall_volumes(Cin,E) {
	Shell_explorer SE(E,Done);
	Shell_entry_const_iterator it;
	CGAL_nef3_forall_shells_of(it,Cin)
	  E.visit_shell_objects(SFace_const_handle(it),SE);
	Point_3 p(E.point(SE.minimal_vertex()));
	Object_handle o;
	if(Cin == E.volumes_begin())
	  o = N.locate(p-vec);
	else
	  o = N.locate(p+vec);
	CGAL_nef3_assertion(assign(Cout,o));
	CGAL_nef3_assertion(Cin == Cout);
      }
    }
  }
  
  void intersection(char* suffix) {

    if(suffix[1] == 'E') {
      Nef_polyhedron N = load_nef3("star.nef3",suffix);
      SNC_explorer E(N.SNCexplorer());
      SNC_intersection is(*E.sncp());

      Point_3 p;
      
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(0,0,0)));
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(2,0,0)));
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(3,0,0)));
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(-1,0,0)));
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(1,1,0)));
      CGAL_nef3_assertion(!is.does_contain_internally(
			    Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			    Point_3(7,25,11)));
      CGAL_nef3_assertion(is.does_contain_internally(
			   Segment_3(Point_3(0,0,0), Point_3(2,0,0)),
			   Point_3(1,0,0)));

      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,0,0), Point_3(0,0,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,0,0), Point_3(0,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,0,0), Point_3(0,0,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,0,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,0,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(0,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,1)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,2,0), Point_3(0,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,3,0), Point_3(0,1,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(0,4,0), Point_3(0,2,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,5,0), Point_3(1,7,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,7,0), Point_3(1,5,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(3,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,3,0), Point_3(0,1,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(!is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,3,0)),
		       	   Segment_3(Point_3(3,0,0), Point_3(1,0,0)), p));
      CGAL_nef3_assertion(is.does_intersect_internally(
			   Segment_3(Point_3(0,1,0), Point_3(0,-1,0)),
		       	   Segment_3(Point_3(1,0,0), Point_3(-1,0,0)), p));
      CGAL_nef3_assertion(p == Point_3(0,0,0));

      Halffacet_const_iterator hf;
    
      CGAL_nef3_assertion(!is.does_contain_internally(
			    E.halffacets_begin(), Point_3(0,0,0)));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,0),Point_3(-31,15,0)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,0),Point_3(-31,15,1)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,-1),Point_3(-31,15,0)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-31,15,0),Point_3(-30,15,0)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-32,15,-1),Point_3(-32,15,0)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(-16,8,-1),Point_3(-16,8,0)),hf,p));
      CGAL_nef3_forall_halffacets(hf, E)
	CGAL_nef3_assertion(!is.does_intersect_internally(
			    Segment_3(Point_3(0,0,-1),Point_3(0,0,1)),hf,p));    

      int i=0;
      CGAL_nef3_forall_halffacets(hf, E) {
	bool b = (i == 2 || i == 3);
	CGAL_nef3_assertion( b == is.does_intersect_internally(
				  Segment_3(Point_3(-31,15,-1),Point_3(-31,15,1)), 
                                  hf, p));
	if(b) CGAL_nef3_assertion(p == Point_3(-31,15,0));
	i++;
      }

      i=0;
      CGAL_nef3_forall_halffacets(hf, E) {
	bool b = (i == 4 || i == 5);
	CGAL_nef3_assertion( b == is.does_intersect_internally(
				   Segment_3(Point_3(-15,7,-1), Point_3(-15,7,1)), 
                                   hf, p));
	if(b) CGAL_nef3_assertion(p == Point_3(-15,7,0));
	i++;
      }
    }
  }

  void point_location_SM(char* suffix) {

    if(suffix[1] == 'E') {  

      Nef_polyhedron N = load_nef3("marks_of_halfspheres2.nef3",suffix);
      CGAL_nef3_assertion(N.is_valid(0,0));
      int i=0;
      SNC_explorer E(N.SNCexplorer());
      Vertex_const_iterator vi = E.vertices_begin();
 
      do {vi++;} while(++i < 8); //  -1 1 0
      CGAL_assertion(E.point(vi) == Point_3(-1,1,0));
      SM_explorer SME(N.SMexplorer(vi));
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 9); //  0 1 -1
      CGAL_assertion(E.point(vi) == Point_3(0,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 10); //  1 1 0
      CGAL_assertion(E.point(vi) == Point_3(1,1,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 11); //  0 1 1
      CGAL_assertion(E.point(vi) == Point_3(0,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 12); //  0 0 -1
      CGAL_assertion(E.point(vi) == Point_3(0,0,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 13); //  0 0 1
      CGAL_assertion(E.point(vi) == Point_3(0,0,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 14); //  -1 0 0
      CGAL_assertion(E.point(vi) == Point_3(-1,0,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 15); //  1 0 0
      CGAL_assertion(E.point(vi) == Point_3(1,0,0));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);
      
      do {vi++;} while(++i < 17); //  -1 1 1
      CGAL_assertion(E.point(vi) == Point_3(-1,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 1);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 18); //  1 1 1
      CGAL_assertion(E.point(vi) == Point_3(1,1,1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 21); //  -1 1 -1
      CGAL_assertion(E.point(vi) == Point_3(-1,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);
      
      do {vi++;} while(++i < 22); //  1 1 -1
      CGAL_assertion(E.point(vi) == Point_3(1,1,-1));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 1);

      do {vi++;} while(++i < 29); //  -0.5 -0.5 -1
      CGAL_assertion(E.point(vi) == Point_3(-1,-1,-2,2));
      SME = N.SMexplorer(vi);
      CGAL_assertion(SME.mark_of_halfsphere(-1) == 0);
      CGAL_assertion(SME.mark_of_halfsphere(+1) == 0);

      N.transform(Aff_transformation_3( 0,-1, 0,
      					1,0, 0,
      					0,0, 1,
      					1));

      N = load_nef3("rotatedcube.nef3",suffix);
      N.transform(Aff_transformation_3( 0,-1, 0,
      					1,0, 0,
      					0,0, 1,
      					1));

      N = load_nef3("two_edges.nef3",suffix);
      N = load_nef3("plane-vertex.nef3",suffix);
      N = load_nef3("is_vertex.nef3",suffix);
    }
  }

  void simplification_SNC(char* suffix) {

    Nef_polyhedron C = load_off("data/cube.off");

    Nef_polyhedron N = load_off("data/cube+v.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube.nef3",suffix));

    N = load_off("data/cube+vee.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube1.nef3",suffix));

    N = load_off("data/cube+veeee.off");
    CGAL_assertion(N.is_valid(0,0));
    //    CGAL_assertion(does_nef3_equals_file(N,"cube2.nef3",suffix)); 

    N = load_off("data/cube+vONe.off");
    CGAL_assertion(N.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N,"cube3.nef3",suffix));

    //    Nef_polyhedron N = C;
    //    Nef_polyhedron N1 = N;
    //    Nef_polyhedron N2 = N;
    //    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,0,0,1)));
    //    N = N.symmetric_difference(N1);
    //    CGAL_assertion(N.is_valid(0,0));
    //    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(0,2,2,1)));
    //    N = N.join(N2);
    N = load_nef3("simplifySface.nef3",suffix);
    CGAL_assertion(N.is_valid(0,0));
    N = N.regularization();
    CGAL_assertion(N.is_valid(0,0));
    //    CGAL_assertion(does_nef3_equals_file(N,"simplifySfaceRef.nef3",suffix));    

    if(suffix[1] == 'E') {   
      N = C;
      CGAL_assertion(N.is_valid(0,0));
      Nef_polyhedron N1 = N;
      Nef_polyhedron N2 = N;
      N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,0,0,1)));
      CGAL_assertion(N2.is_valid(0,0));
      N1 = N1.intersection(N2);
      CGAL_assertion(N1.is_valid(0,0));
      N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,0,0,1)));
      CGAL_assertion(N1.is_valid(0,0));
      N2 = load_nef3("viereck.nef3", suffix);
      CGAL_assertion(N2.is_valid(0,0));
      N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,3,0,1)));
      CGAL_assertion(N2.is_valid(0,0));
      N2 = N1.difference(N2);
      CGAL_assertion(N2.is_valid(0,0));
      N = N.difference(N1);
      CGAL_assertion(N.is_valid(0,0));
      N = N.join(N2);
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"donotsimplify.nef3",suffix));
    } 
     
    if(suffix[1] == 'E') {   
      N = load_nef3("unionfind.nef3", suffix);
      CGAL_assertion(N.is_valid(0,0));
      Nef_polyhedron N1 =  N;
      N.transform(Aff_transformation_3(1, 0, 0, 
				       0, 0,-1,
				       0, 1, 0,
				       1));
      CGAL_assertion(N1.is_valid(0,0));
      N1.transform(Aff_transformation_3(1, 0, 0, 
					0, 0,-1,
					0, 1, 0,
					1));
      CGAL_assertion(N1.is_valid(0,0));

      Nef_polyhedron N2 = N.symmetric_difference(N1);
      CGAL_assertion(N2.is_valid(0,0));
      N2 = N2.regularization();
      CGAL_assertion(N2.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(N2,"unionfindref.nef3",suffix));

      N2 = N.join(N1);
      CGAL_assertion(N2.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N2,"unionfindref.nef3",suffix));
    }  

    if(suffix[1] == 'E') {
      N = load_nef3("topologyE.nef3",suffix);
      Nef_polyhedron N1 = load_nef3("single_vertex.nef3",suffix);
      N = N.join(N1);
      N = N.closure();
    }

    if(suffix[1] == 'E') {
      N = load_nef3("nestedFCafterSimpl.nef3",suffix);
      N = N.closure();
    }
  
  }

  void simplification_SM(char* suffix) {

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
    CGAL_assertion(does_nef3_equals_file(N1,"cube+plane.nef3",suffix));
 
    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,0,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,2,1)));
    N1 = N.join(N1);
    CGAL_assertion(N1.is_valid(0,0));
    //    CGAL_assertion(does_nef3_equals_file(N1,"cube+line.nef3",suffix));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,2,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,0,1)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    //    CGAL_assertion(does_nef3_equals_file(N1,"cube+vertex1.nef3",suffix));

    N  = C;
    N1 = N;
    N2 = N;
    N2.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(2,2,2,1)));
    N1 = N1.intersection(N2);
    CGAL_assertion(N1.is_valid(0,0));
    N1.transform(Aff_transformation_3( CGAL::TRANSLATION, Vector_3(-1,-1,-1,1)));
    N1 = N.symmetric_difference(N1);
    CGAL_assertion(N1.is_valid(0,0));
    CGAL_assertion(does_nef3_equals_file(N1,"cube+vertex2.nef3",suffix));
  }

  void pluecker_coordinates(char* suffix) {

    if(suffix[1] == 'E') {
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
      CGAL_assertion(does_nef3_equals_file(N,"cube_created_from_halfspaces.nef3", suffix));
    }
  }

  void unary_operations(char* suffix) {
    
    if(suffix[1] == 'E') {
      Nef_polyhedron T("data/topologyE.nef3");
      Nef_polyhedron N = T;
      N = N.boundary();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"boundary.nef3",suffix));
    
      N = T;
      N = N.interior();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"interior.nef3",suffix));
      
      N = T;
      N = N.closure();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"closure.nef3",suffix));
      
      N = T;
      N = N.regularization();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"regularization.nef3",suffix));
      
      N = T; 
      N = N.complement();
      CGAL_assertion(N.is_valid(0,0));
      CGAL_assertion(does_nef3_equals_file(N,"complement.nef3",suffix));
    }
  }
    
  void mark_evaluation(char* suffix) {
    
    if(suffix[1] == 'E') {    
      Nef_polyhedron R;
      Nef_polyhedron P = load_nef3("mark_eval.nef3", suffix);
      Nef_polyhedron N = load_nef3("mark_eval2.nef3", suffix);
      Nef_polyhedron P2 = P;
      Nef_polyhedron N2 = N;
      P2.transform(Aff_transformation_3( CGAL::SCALING, 2, 1));
      N2.transform(Aff_transformation_3( CGAL::SCALING, 2, 1));
      
      R = P2.intersection(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePPinters.nef3", suffix));
      R = P2.intersection(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePNinters.nef3", suffix));
      R = N2.intersection(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNPinters.nef3", suffix));
      R = N2.intersection(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNNinters.nef3", suffix));
      
      R = P2.join(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePPjoin.nef3", suffix));
      R = P2.join(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePNjoin.nef3", suffix));
      R = N2.join(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNPjoin.nef3", suffix));
      R = N2.join(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNNjoin.nef3", suffix));

      R = P2.difference(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePPdiff.nef3", suffix));
      R = P2.difference(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePNdiff.nef3", suffix));
      R = N2.difference(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNPdiff.nef3", suffix));
      R = N2.difference(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNNdiff.nef3", suffix));
      
      R = P2.symmetric_difference(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePPsymdiff.nef3", suffix));
      R = P2.symmetric_difference(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "mePNsymdiff.nef3", suffix));
      R = N2.symmetric_difference(P);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNPsymdiff.nef3", suffix));
      R = N2.symmetric_difference(N);
      CGAL_nef3_assertion(R.is_valid(0,0));
      //      CGAL_assertion(does_nef3_equals_file(R, "meNNsymdiff.nef3", suffix));
    }
  } 
    
public:
  void run_test(char* suffix) {
    
        loadSave(suffix);
        construction(suffix);
        point_location_SNC(suffix);
        intersection(suffix);
        point_location_SM(suffix);
        simplification_SNC(suffix);
        simplification_SM(suffix);
        pluecker_coordinates(suffix);
        unary_operations(suffix);
        mark_evaluation(suffix);
  }

};

template<typename Kernel>
const char* test<Kernel>::datadir="data/";

int main() {

  typedef CGAL::Gmpz                         NT;
  typedef CGAL::Simple_homogeneous<NT>       SH_Kernel;
  typedef CGAL::Extended_homogeneous_3<NT>   EH_Kernel;
  
  CGAL::Timer t;
  t.start();

  test<SH_Kernel> test_SH;
  test<EH_Kernel> test_EH;

  test_SH.run_test(".SH");
  test_EH.run_test(".EH");

  t.stop();
  std::cout << "Time " << t.time() << std::endl;
}

