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
// file          : include/CGAL/Nef_3/Nef_polyhedron_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Nef polyhedron in the space
// ============================================================================
#ifndef CGAL_NEF_POLYHEDRON_3_H
#define CGAL_NEF_POLYHEDRON_3_H
#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#define SNC_VISUALIZOR
// #define SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#ifdef SNC_VISUALIZOR
#include <CGAL/Nef_3/SNC_visualizor_OGL.h>
#endif // SNC_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_decorator.h>
#include <CGAL/Nef_3/SNC_SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_io_parser.h>
#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR

#include <CGAL/Nef_3/polyhedron_3_to_nef_3.h>
#include <CGAL/Polyhedron_3.h>

#include <list>

#undef _DEBUG
#define _DEBUG 11
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename T> class Nef_polyhedron_3;
template <typename T> class Nef_polyhedron_3_rep;

template <typename T>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_3<T>&); 
template <typename T>
std::istream& operator>>(std::istream&, Nef_polyhedron_3<T>&); 

template <typename T>
class Nef_polyhedron_3_rep : public Rep
{ 
  typedef Nef_polyhedron_3_rep<T>                     Self;
  friend class Nef_polyhedron_3<T>;
  typedef CGAL::SNC_structure<T>                       SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>           SNC_decorator;
  typedef CGAL::SNC_constructor<SNC_structure>         SNC_constructor;
  typedef CGAL::SNC_point_locator<SNC_structure>       SNC_point_locator;
  typedef CGAL::SNC_io_parser<SNC_structure>           SNC_io_parser;
#ifdef SNC_VISUALIZOR
  typedef CGAL::SNC_visualizor_OGL<SNC_structure>      SNC_visualizor;
#endif // SNC_VISUALIZOR
  typedef CGAL::SNC_SM_decorator<SNC_structure>        SM_decorator;
  typedef CGAL::SNC_SM_const_decorator<SNC_structure>  SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<SNC_structure>        SM_overlayer;
  typedef CGAL::SNC_SM_point_locator<SNC_structure>    SM_point_locator;
  typedef CGAL::SNC_SM_io_parser<SNC_structure>        SM_io_parser;
#ifdef SM_VISUALIZOR
  typedef CGAL::SNC_SM_visualizor<SNC_structure>       SM_visualizor;
#endif // SM_VISUALIZOR

  SNC_structure snc_;
  // SNC_point_locator* pl_;
  
 public:
  Nef_polyhedron_3_rep() : snc_() {}
  ~Nef_polyhedron_3_rep() { snc_.clear(); }
};

/*{\Manpage {Nef_polyhedron_3} {T} {Nef Polyhedra in Space}{N}}*/

/*{\Mdefinition
An instance of data type |\Mname| is a subset of 3-space which is the
result of forming complements and intersections starting from a set |H| of
halfspaces. |\Mtype| is closed under all binary set opertions |intersection|,
|union|, |difference|, |complement| and under the topological operations
|boundary|, |closure|, and |interior|.}*/

template <typename T>
class Nef_polyhedron_3 : public CGAL::Handle_for< Nef_polyhedron_3_rep<T> >
{ 
public:
typedef T Extended_kernel; 
static  T EK; // static extended kernel

  /*{\Mtypes 7}*/
  typedef Nef_polyhedron_3<T>    Self;
  typedef Handle_for< Nef_polyhedron_3_rep<T> > Base;
  typedef typename T::Kernel     Kernel;
  typedef typename T::Point_3    Point_3;
  typedef typename T::Plane_3    Plane_3;

  typedef bool Mark;
  /*{\Xtypemember marking set membership or exclusion.}*/

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  enum Content { EMPTY=0, COMPLETE=1 };
  /*{\Menum construction selection}*/

protected:
  struct AND { bool operator()(bool b1, bool b2) const { return b1&&b2; } };
  struct OR { bool operator()(bool b1, bool b2) const { return b1||b2; } };
  struct DIFF { bool operator()(bool b1, bool b2) const { return b1&&!b2; } };
  struct XOR { bool operator()(bool b1, bool b2) const 
    { return (b1&&!b2)||(!b1&&b2); } };

  typedef Nef_polyhedron_3_rep<T>               Nef_rep;
  typedef typename Nef_rep::SNC_structure       SNC_structure;
  typedef typename Nef_rep::SNC_decorator       SNC_decorator;
  typedef typename Nef_rep::SNC_constructor     SNC_constructor;
  typedef typename Nef_rep::SNC_point_locator   SNC_point_locator;
  typedef typename Nef_rep::SNC_io_parser       SNC_io_parser;
#ifdef SNC_VISUALIZOR
  typedef typename Nef_rep::SNC_visualizor      SNC_visualizor;
#endif // SNC_VISUALIZOR
  typedef typename Nef_rep::SM_decorator        SM_decorator;
  typedef typename Nef_rep::SM_const_decorator  SM_const_decorator;
  typedef typename Nef_rep::SM_overlayer        SM_overlayer;
  typedef typename Nef_rep::SM_point_locator    SM_point_locator;
  typedef typename Nef_rep::SM_io_parser        SM_io_parser;
#ifdef SM_VISUALIZOR
  typedef typename Nef_rep::SM_visualizor       SM_visualizor;
#endif // SM_VISUALIZOR

  SNC_structure& snc() { return ptr()->snc_; } 
  const SNC_structure& snc() const { return ptr()->snc_; } 

  friend std::ostream& operator<< <>
      (std::ostream& os, const Nef_polyhedron_3<T>& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Nef_polyhedron_3<T>& NP);

  typedef typename SNC_decorator::Vertex_handle    Vertex_handle;
  typedef typename SNC_decorator::Halfedge_handle  Halfedge_handle;
  typedef typename SNC_decorator::Halffacet_handle
                                                   Halffacet_handle;
  typedef typename SNC_decorator::Volume_handle    Volume_handle;
  typedef typename SNC_structure::Vertex_const_handle
                                                   Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle
                                                   Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle
                                                   Halffacet_const_handle;
  typedef typename SNC_structure::Volume_const_handle
                                                   Volume_const_handle;
  typedef typename SM_decorator::SVertex_handle    SVertex_handle;
  typedef typename SM_decorator::SHalfedge_handle  SHalfedge_handle;
  typedef typename SM_decorator::SFace_handle      SFace_handle;
  typedef typename SM_decorator::SVertex_const_handle
                                                   SVertex_const_handle;
  typedef typename SM_decorator::SHalfedge_const_handle
                                                   SHalfedge_const_handle;
  typedef typename SM_decorator::SFace_const_handle
                                                   SFace_const_handle;
  typedef typename SNC_decorator::Vertex_iterator  Vertex_iterator;
  typedef typename SNC_decorator::Halfedge_iterator
                                                   Halfedge_iterator;
  typedef typename SNC_decorator::Halffacet_iterator
                                                   Halffacet_iterator;
  typedef typename SNC_structure::Shell_entry_iterator
                                                   Shell_entry_iterator;
  typedef typename SNC_decorator::Volume_iterator  Volume_iterator;
  typedef typename SNC_structure::Vertex_const_iterator
                                                   Vertex_const_iterator;
  typedef typename SNC_structure::Halfedge_const_iterator 
                                                   Halfedge_const_iterator;
  typedef typename SNC_structure::Halffacet_const_iterator     
                                                   Halffacet_const_iterator;
  typedef typename SNC_structure::Volume_const_iterator     
                                                   Volume_const_iterator;
  typedef typename SNC_structure::Shell_entry_const_iterator
                                                   Shell_entry_const_iterator;
  typedef typename SM_decorator::SVertex_iterator  SVertex_iterator;
  typedef typename SM_decorator::SHalfedge_iterator
                                                   SHalfedge_iterator;
  typedef typename SM_decorator::SFace_iterator    SFace_iterator;
  typedef typename SM_decorator::SVertex_const_iterator
                                                   SVertex_const_iterator;
  typedef typename SM_decorator::SHalfedge_const_iterator 
                                                   SHalfedge_const_iterator;
  typedef typename SM_decorator::SFace_const_iterator     
                                                   SFace_const_iterator;

  typedef typename SNC_structure::Sphere_point     Sphere_point;
  typedef typename SNC_structure::Sphere_segment   Sphere_segment;
  typedef typename SNC_structure::Sphere_circle    Sphere_circle;
  
  Vertex_handle create_cube_corner(int x, int y, int z, bool space=false) {
    CGAL_assertion(CGAL_NTS abs(x) == 
		   CGAL_NTS abs(y) == 
		   CGAL_NTS abs(z) == 1);
    Vertex_handle v = snc().new_vertex();
    snc().point(v) = Point_3(x*IMMN, y*IMMN, z*IMMN); 
    /* TODO: to replace the IMMN constant by a real infimaximal number */
    SM_decorator D(v);
    Sphere_point sp[] = { Sphere_point(-x, 0, 0), 
			  Sphere_point(0, -y, 0), 
			  Sphere_point(0, 0, -z) };
    /* create box vertices */
    SVertex_handle sv[3];
    for(int vi=0; vi<3; ++vi)
      sv[vi] = D.new_vertex(sp[vi]);
    /* create facet's edge uses */
    Sphere_segment ss[3];
    SHalfedge_handle she[3];
    for(int si=0; si<3; ++si) {
      she[si] = D.new_edge_pair(sv[si], sv[(si+1)%3]);
      ss[si] = Sphere_segment(sp[si],sp[(si+1)%3]);
      D.circle(she[si]) = ss[si].sphere_circle();
      D.circle(D.twin(she[si])) = ss[si].opposite().sphere_circle();
    }
    /* create facets */
    SFace_handle fi = D.new_face();
    SFace_handle fe = D.new_face();
    D.link_as_face_cycle(she[0], fi);
    D.link_as_face_cycle(D.twin(she[0]), fe);
    /* set boundary marks */
    SHalfedge_iterator e = D.shalfedges_begin();
    SFace_handle f;
    Sphere_point p1 = D.point(D.source(e));
    Sphere_point p2 = D.point(D.target(e));
    Sphere_point p3 = D.point(D.target(D.next(e)));
    if ( spherical_orientation(p1,p2,p3) > 0 )
      f = D.face(e);
    else
      f = D.face(D.twin(e));
    D.mark(f) = space;
    CGAL_nef3_forall_sedges_of(e, v)
      D.mark(e) = D.mark(D.source(e)) = space;
    D.mark_of_halfsphere(-1) = (x<0 && y>0 && z>0);
    D.mark_of_halfsphere(+1) = (x>0 && y>0 && z<0);
    return v;
  }

  void initialize_simple_cube_vertices(Content space) {
    create_cube_corner( 1, 1, 1, space );
    create_cube_corner(-1, 1, 1, space );
    create_cube_corner( 1,-1, 1, space );
    create_cube_corner(-1,-1, 1, space );
    create_cube_corner( 1, 1,-1, space );
    create_cube_corner(-1, 1,-1, space );
    create_cube_corner( 1,-1,-1, space );
    create_cube_corner(-1,-1,-1, space );
  }

  void check_h_for_intersection_of_12_cube_edges_and_add_vertices
  (const Plane_3& p);
  void create_intersection_vertex_of_h_and_e();
  void init_cube_vertices_depending_on_h(const Plane_3& p);
  void add_h_to_local_view_of_v();
  
  void build_external_structure() {
    SNC_constructor C(snc());
    C.pair_up_halfedges();
    C.link_shalfedges_to_facet_cycles();
    C.categorize_facet_cycles_and_create_facets();
    C.create_volumes();
  }

  void clear_box_marks() 
    /* unset all frame marks */ {
    SNC_decorator D(snc());
    Volume_iterator c = D.volumes_begin(); 
    D.mark(c) = false;
    D.clear_outer_box_marks();
  }

public:
  public:
  /*{\Mcreation 3}*/

  Nef_polyhedron_3(Content space = EMPTY);
  /*{\Mcreate creates an instance |\Mvar| of type |\Mname|
  and initializes it to the empty set if |space == EMPTY|
  and to the whole space if |space == COMPLETE|.}*/

  Nef_polyhedron_3(const Plane_3& p, Boundary b = INCLUDED); 
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the
  halfspace on the positive side of |p| including |p| if |b==INCLUDED|,
  excluding |p| if |b==EXCLUDED|.}*/

  Nef_polyhedron_3(const Nef_polyhedron_3<T>& N1) : Base(N1) {}
  Nef_polyhedron_3& operator=(const Nef_polyhedron_3<T>& N1)
  { Base::operator=(N1); return (*this); }
  ~Nef_polyhedron_3() {}
  
  typedef Polyhedron_3< Kernel> Polyhedron;
  Nef_polyhedron_3( Polyhedron& P) {
    initialize_simple_cube_vertices(EMPTY);
    polyhedron_3_to_nef_3< Polyhedron, SNC_structure, SNC_constructor>
      ( P, snc() );
    build_external_structure();
    simplify();
  }
  
  void dump() { SNC_io_parser::dump( snc()); }

  void visualize() { 
#ifdef SNC_VISUALIZOR
    SNC_visualizor sncv( snc());
    sncv.draw();
    OGL::polyhedra_.back().debug();
    OGL::start_viewer();
#endif // SNC_VISUALIZOR
  }

 protected:
  void clone_rep() { *this = Nef_polyhedron_3<T>(snc()); }

  Nef_polyhedron_3(const SNC_structure& H, bool cloneit=true);
  /*{\Xcreate makes |\Mvar| a new object.  If |cloneit==true| then the
  underlying structure of |H| is copied into |\Mvar|.}*/

  /*{\Moperations 4 3 }*/

  void clear(Content space = EMPTY)
  { *this = Nef_polyhedron_3(space); }
  /*{\Mop makes |\Mvar| the empty set if |space == EMPTY| and the
  full space if |space == COMPLETE|.}*/

  bool is_empty() //const
  /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
  { //SNC_structure snc_non_const = snc(); // const_decorator not implemented
    SNC_decorator D(snc());
    Volume_iterator v = D.volumes_begin(); // it should be const_iterator
    return (D.number_of_vertices()==8 &&
	    D.number_of_edges()==12 &&
	    D.number_of_facets()==6 &&
	    D.number_of_volumes()==2 &&
	    D.mark(++v) == false);
  }

  bool is_space() //const
  /*{\Mop returns true if |\Mvar| is the whole space, false otherwise.}*/
  //{ SM_const_decorator D(pm());
  //  typename PM_const_decorator::Volume_const_iterator v = D.volumes_begin();
  { //SNC_structure snc_non_const = snc(); // const_decorator not implemented
    SNC_decorator D(snc());
    Volume_iterator v = D.volumes_begin(); // it should be const_iterator
    return (D.number_of_vertices()==8 &&
            D.number_of_edges()==12 &&
            D.number_of_facets()==6 && 
	    D.number_of_volumes()==2 && 
            D.mark(++v) == true);
  }

  /*{\Xtext \headerline{Destructive Operations}}*/

 private:
  void simplify() {
    snc().simplify();
  }
  
 public:
  void extract_complement();
  /*{\Xop converts |\Mvar| to its complement. }*/
  void extract_interior();
  /*{\Xop converts |\Mvar| to its interior. }*/
  void extract_boundary();
  /*{\Xop converts |\Mvar| to its boundary. }*/

  void extract_closure()
  /*{\Xop converts |\Mvar| to its closure. }*/
  { TRACEN("extract closure");
    if (refs()>1) *this = Nef_polyhedron_3<T>(snc()); // clone
    extract_complement();
    extract_interior();
    extract_complement();
  }

  void extract_regularization()
  /*{\Xop converts |\Mvar| to its regularization. }*/
  { TRACEN("extract regularization");
    if (refs()>1) *this = Nef_polyhedron_3<T>(snc()); // clone
    extract_interior();
    extract_closure();
  }

  /*{\Mtext \headerline{Constructive Operations}}*/

  Nef_polyhedron_3<T> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_complement();
    return res;
  }

  Nef_polyhedron_3<T> interior() const
  /*{\Mop    returns the interior of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_3<T> closure() const
  /*{\Mop returns the closure of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_3<T> boundary() const
  /*{\Mop returns the boundary of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_3<T> regularization() const
  /*{\Mop    returns the regularized polyhedron (closure of 
             the interior).}*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_regularization();
    return res;
  }


  /* HERE I AM!!!

  void binop( Self& P0, Self& P1) {
    // for each vertex x of P_i do
    //   qualify_x_with_respect_to P(1-i)
    //   biopn both local views Px(0,1): subdivide-slect-simplify...
    // for each edge-edge or edge-face interesection x do
    //   binop both local views Px(0,1): subdivide-select-simplify...
    // remove vertices whose local view is not that of a vertex
    // synthesis f spatial structure...
    } */

  Nef_polyhedron_3<T> intersection(const Nef_polyhedron_3<T>& N1) const
  /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/
  { Nef_polyhedron_3<T> res(snc(),false); // empty, no frame
  //EW_overlayer EWO(res.snc());
  //EWO.subdivide(snc(),N1.snc());
  //AND _and; EWO.select(_and);
  //EWO.simplify();
  //EWO.build_external_structure();
    res.clear_box_marks();
    return res;
  }

  Nef_polyhedron_3<T> join(const Nef_polyhedron_3<T>& N1) const
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/
  { Nef_polyhedron_3<T> res(snc(),false); // empty, no frame
  //EW_overlayer EWO(res.snc());
  //EWO.subdivide(snc(),N1.snc());
  //OR _or; EWO.select(_or);
  //EWO.simplify();
  //EWO.build_external_structure();
    res.clear_box_marks();
    return res;
  }

  Nef_polyhedron_3<T> difference(
    const Nef_polyhedron_3<T>& N1) const
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/
  { Nef_polyhedron_3<T> res(snc(),false); // empty, no frame
  //EW_overlayer EWO(res.snc());
  //EWO.subdivide(snc(),N1.snc());
  //DIFF _diff; EWO.select(_diff);
  //EWO.simplify();
  //EWO.build_external_structure();
    res.clear_box_marks();
    return res;
  }    

  Nef_polyhedron_3<T> symmetric_difference(
    const Nef_polyhedron_3<T>& N1) const
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/
  { Nef_polyhedron_3<T> res(snc(),false); // empty, no frame
  //EW_overlayer EWO(res.snc());
  //EWO.subdivide(snc(),N1.snc());
  //XOR _xor; EWO.select(_xor);
  //EWO.simplify();
  //EWO.build_external_structure();
    res.clear_box_marks();
    return res;
  }

  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement}. There are also the corresponding
  modification operations |*=,+=,-=,^=|.}*/

  Nef_polyhedron_3<T>  operator*(const Nef_polyhedron_3<T>& N1) const
  { return intersection(N1); }

  Nef_polyhedron_3<T>  operator+(const Nef_polyhedron_3<T>& N1) const
  { return join(N1); }

  Nef_polyhedron_3<T>  operator-(const Nef_polyhedron_3<T>& N1) const
  { return difference(N1); }

  Nef_polyhedron_3<T>  operator^(const Nef_polyhedron_3<T>& N1) const
  { return symmetric_difference(N1); }

  Nef_polyhedron_3<T>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_3<T>& operator*=(const Nef_polyhedron_3<T>& N1)
  { this = intersection(N1); return *this; }

  Nef_polyhedron_3<T>& operator+=(const Nef_polyhedron_3<T>& N1)
  { this = join(N1); return *this; }

  Nef_polyhedron_3<T>& operator-=(const Nef_polyhedron_3<T>& N1)
  { this = difference(N1); return *this; }

  Nef_polyhedron_3<T>& operator^=(const Nef_polyhedron_3<T>& N1)
  { this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality.}*/

  bool operator==(const Nef_polyhedron_3<T>& N1) const
  { return symmetric_difference(N1).is_empty(); }

  bool operator!=(const Nef_polyhedron_3<T>& N1) const
  { return !operator==(N1); }  

  bool operator<(const Nef_polyhedron_3<T>& N1) const
  { return !N1.difference(*this).is_empty() && difference(N1).is_empty(); } 

  bool operator>(const Nef_polyhedron_3<T>& N1) const   
  { return difference(*this).is_empty() && !difference(N1).is_empty(); } 

  bool operator<=(const Nef_polyhedron_3<T>& N1) const
  { return difference(N1).is_empty(); } 

  bool operator>=(const Nef_polyhedron_3<T>& N1) const
  { return N1.difference(*this).is_empty(); } 

  /*{\Mtext \headerline{Exploration}
  As Nef polyhedra are the result of forming complements 
  and intersections starting from a set |H| of halfspaces which are
  defined by oriented planes in three space. The corresponding 
  structure is represented by an extended wuerzburg structure 
  $W = (V,E,F,C)$. For topological queries within |W| the following 
  types and operations allow exploration access to this structure.}*/

  /*{\Mtypes 3}*/
  //typedef typename Nef_rep::SNC_decorator SNC_decorator; // defined above

  //typedef CGAL::SNC_explorer<SNC_decorator,T> SNC_explorer; //not implemented
  /*{\Mtypemember a decorator to examine the underlying plane map. 
  See the manual page of |EW_explorer|.}*/

  typedef typename SNC_point_locator::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying
  plane map. The kind of object |(vertex, halfedge, face)| can 
  be determined and the object can be assigned to a corresponding
  handle by the three functions:\\
  |bool assign(Vertex_const_handle& h, Object_handle)|\\
  |bool assign(Edge_const_handle& h, Object_handle)|\\
  |bool assign(Facet_const_handle& h, Object_handle)|\\
  |bool assign(Volume_const_handle& h, Object_handle)|\\
  where each function returns |true| iff the assignment to
  |h| was done.}*/

  /*{\Moperations 3 1 }*/

  bool contains(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the set
  represented by |\Mvar|.}*/
  { SNC_point_locator PL(snc()); return PL.mark(h); }

  bool contained_in_boundary(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the $2$-skeleton
  of |\Mvar|.}*/
  { Vertex_const_handle v;
    Halfedge_const_handle e;
    Halffacet_const_handle f;
    return  ( assign(v,h) || assign(e,h) || assign(f,h) );
  }

  Object_handle locate(const Point_3& p) const;
  /*{\Mop  returns a generic handle |h| to an object (vertex, edge, facet,
  volume) of the underlying SNC which contains the point |p| in its relative 
  interior. The point |p| is contained in the set represented by |\Mvar| if 
  |\Mvar.contains(h)| is true.}*/

  //EW_explorer explorer() const // not implemented, replaced by decorator?
  /*{\Mop returns a decorator object which allows read-only access of
  the underlying three-dimensional subdivision structure. 
  See the manual page |EW_explorer| for its usage.}*/
  //{ return EW_explorer(snc,EK); }

  /*{\Mimplementation Nef polyhedra are implemented on top of an
  extended Wuerzburg structure data structure (EWS) and use linear
  space in the number of vertices, edges and facets.  Operations like
  empty take constant time. The operations clear, complement, interior,
  closure, boundary, regularization, input and output take linear
  time. All binary set operations and comparison operations take time
  $O(N^2)$ where $N$ is the size of the output plus the size of the
  input.

  The point location operations run in linear query time without any 
  preprocessing.}*/

  /*{\Mexample Nef polyhedra are parameterized by a so called extended
  geometric kernel. There's currently only one such kernel based on a
  homogeneous representation of extended points called
  |Extended_homogeneous_3<NT>|.  The kernel is parameterized by a
  multiprecision integer type. The member types of |Nef_polyhedron_3<
  Extended_homogeneous_3<NT> >| map to corresponding types of the CGAL
  geometry kernel (e.g. |Nef_polyhedron::Plane_3| equals
  |CGAL::Homogeneous<leda_integer>::Plane_3| in the example below).
  \begin{Mverb}
  #include <CGAL/basic.h>
  #include <CGAL/leda_integer.h>
  #include <CGAL/Extended_homogeneous_3.h>
  #include <CGAL/Nef_polyhedron_3.h>

  using namespace CGAL;
  typedef  Extended_homogeneous_3<leda_integer> Extended_kernel;
  typedef  Nef_polyhedron_3<Extended_kernel>    Nef_polyhedron;
  typedef  Nef_polyhedron::Plane_3              Plane_3;

  int main()
  {
    Nef_polyhedron N1(Plane_3(1,0,0,0));
    Nef_polyhedron N2(Plane_3(0,1,0,0), Nef_polyhedron::EXCLUDED);
    Nef_polyhedron N3 = N1 * N2; // line (*)
    return 0;
  }
  \end{Mverb}
  After line (*) |N3| is the intersection of |N1| and |N2|.}*/


}; // end of Nef_polyhedron_3

template <typename T>
T Nef_polyhedron_3<T>::EK = T();

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3(Content space) : Base(Nef_rep()) {
  TRACEN("construction from empty or space.");
  initialize_simple_cube_vertices(space);
  build_external_structure();
}

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3(const Plane_3& h, Boundary b) : Base(Nef_rep()) {
  TRACEN("construction from plane "<<h);
  initialize_simple_cube_vertices(space);
  //add_box_corners(h, b);
  build_external_structure();
}

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3(const SNC_structure& W, bool clone) : Base(Nef_rep()) {
  TRACEN("construction from an existing nef3");
  if (clone) { 
    snc() = W;
  }
}

template <typename T>
void Nef_polyhedron_3<T>::extract_complement() {
  TRACEN("extract complement");
  if( is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_nef3_forall_vertices(v,D) {
    TRACE(D.mark(v));
    D.mark(v) = !D.mark(v); 
    TRACEN(" vcomplement "<<D.mark(v));
  }
  Halfedge_iterator e;
  CGAL_nef3_forall_edges(e,D) D.mark(e) = !D.mark(e);
  Halffacet_iterator f;
  CGAL_nef3_forall_facets(f,D) {
    TRACE(D.mark(f)<<f->mark_<<f->twin_->mark_);
    D.mark(f) = !D.mark(f); 
    TRACEN(" fcomplement "<<D.mark(f)<<", "<<&*f);
  }
  Volume_iterator c;
  CGAL_nef3_forall_volumes(c,D) D.mark(c) = !D.mark(c);
  clear_box_marks();
  dump();
}


template <typename T>
void Nef_polyhedron_3<T>::extract_interior() {
  TRACEN("extract interior");
  //if (refs()>1) *this = Nef_polyhedron_3<T>(snc()); // clone
  //EW_overlayer D(snc());
  Vertex_iterator v;
  //forall_vertices(v,D) D.mark(v) = false;
  Halfedge_iterator e;
  //forall_edges(e,D) D.mark(e) = false;
  Halffacet_iterator f;
  //forall_facets(f,D) D.mark(f) = false;
  Volume_iterator c;
    //forall_volumes(c,D) D.mark(c) = true;
  //D.simplify();
}

template <typename T>
void Nef_polyhedron_3<T>::extract_boundary() {
  TRACEN("extract boundary");
  //if (refs()>1) *this = Nef_polyhedron_3<T>(snc()); // clone
  //EW_overlayer D(snc());
  Vertex_iterator v;
  //forall_vertices(v,D) D.mark(v) = true;
  Halfedge_iterator e;
  //forall_edges(e,D) D.mark(e) = true;
  Halffacet_iterator f;
  //forall_facets(f,D) D.mark(f) = true;
  Volume_iterator c;
    //forall_volumes(c,D) D.mark(c) = false;
  clear_box_marks();
  //D.simplify();
}

template <typename T>
std::ostream& operator<<
 (std::ostream& os, const Nef_polyhedron_3<T>& NP)
{
  typedef typename Nef_polyhedron_3<T>::SNC_decorator SNC_decorator;
  CGAL::SNC_io_parser<SNC_decorator> O(os, NP.snc());
  O.print();
  return os;
}

template <typename T>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_3<T>& NP)
{
  typedef typename Nef_polyhedron_3<T>::SNC_decorator SNC_decorator;
  CGAL::SNC_io_parser<SNC_decorator> I(is, NP.snc());
  I.read();
  I.check_integrity();
  return is;
}

CGAL_END_NAMESPACE

#endif //CGAL_NEF_POLYHEDRON_3_H


