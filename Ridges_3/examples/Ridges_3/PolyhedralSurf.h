#ifndef _POLYHEDRALSURF_H_
#define _POLYHEDRALSURF_H_

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <algorithm>
#include <vector>
#include <list>

#include <stdlib.h>
#include <stdio.h>

#include "PolyhedralSurf_operations.h" 

//----------------------------------------------------------------
// A redefined items class for the Polyhedron_3 with 
// a refined vertex class that contains monge data and ring_tag
// a refined facet with a normal vector + tag is_visited
// a refined halfedge with length
template < class Refs, class Tag, class Pt, class FGeomTraits > 
class My_vertex:public CGAL::HalfedgeDS_vertex_base < Refs, Tag, Pt >
{
 protected:
  //to be def in Vertex_wrapper too from the Traits
  typedef typename FGeomTraits::FT FT;
  typedef typename FGeomTraits::Point_3 Point_3;
  typedef typename FGeomTraits::Vector_3 Vector_3; 
   
 protected:
  //monge data 
  Vector_3 m_d1; //max ppal dir
  Vector_3 m_d2; //min ppal dir; monge normal is then n_m=d1^d2, should be so that n_m.n_mesh>0
  FT m_k1, m_k2; //max/min ppal curv
  FT m_b0, m_b3; //blue/red extremalities
  FT m_P1, m_P2; //if fourth order quantities
    
  //this is for collecting i-th ring neighbours
  int ring_index;

 public:
  //monge data
  const Vector_3 d1() const { return m_d1; }
  Vector_3& d1() { return m_d1; }
  const Vector_3 d2() const { return d2; }
  Vector_3& d2() { return m_d2; }
  const FT k1() const { return m_k1; }
  FT& k1() { return m_k1; }
  const FT k2() const { return m_k2; }
  FT& k2() { return m_k2; }
  const FT b0() const { return m_b0; }
  FT& b0() { return m_b0; }
  const FT b3() const { return m_b3; }
  FT& b3() { return m_b3; }
  const FT P1() const { return m_P1; }
  FT& P1() { return m_P1; }//= 3*b1^2+(k1-k2)(c0-3k1^3)
  const FT P2() const { return m_P2; }
  FT& P2() { return m_P2; }//= 3*b2^2+(k2-k1)(c4-3k2^3)

  //this is for collecting i-th ring neighbours
  void setRingIndex(int i) { ring_index = i; }
  int getRingIndex() { return ring_index; }
  void resetRingIndex() { ring_index = -1; }
 
  My_vertex(const Point_3 & pt):
    CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt),
    ring_index(-1)    {}
  My_vertex()    {} 
};

//----------------------------------------------------------------
// Facet with normal and possibly more types. types are recovered
//from the FGeomTraits template arg
// tag for ridge computations
//----------------------------------------------------------------
template < class Refs, class Tag, class FGeomTraits >
class My_facet:public CGAL::HalfedgeDS_face_base < Refs, Tag >
{
public:
 typedef typename FGeomTraits::Vector_3 Vector_3;

protected:
  Vector_3 normal;
  //  int ring_index;
  bool m_is_visited;
public:
  My_facet() {}//: ring_index(-1) {}
  Vector_3 & getUnitNormal() { return normal; }
  void setNormal(Vector_3 n) { normal = n; }

/*   //this is for collecting i-th ring neighbours */
/*   void setRingIndex(int i) { ring_index = i; } */
/*   int getRingIndex() { return ring_index; } */
/*   void resetRingIndex() { ring_index = -1; } */

  //this is for following ridge lines
  void set_visited(bool b) { m_is_visited = b; }
  const bool is_visited() const { return m_is_visited;}
  void reset_is_visited() { m_is_visited = false; }
};

//----------------------------------------------------------------
// Halfedge with length
//----------------------------------------------------------------
template < class Refs, class Tprev, class Tvertex, class Tface > 
class My_halfedge:public CGAL::HalfedgeDS_halfedge_base < Refs, Tprev, Tvertex, Tface >
{
protected:
  //  int ring_index;
  double len;
public:
  My_halfedge() {}//: ring_index(-1) {}
 /*  void setRingIndex(int i) {	ring_index = i;    } */
/*   int getRingIndex() {return ring_index;    } */
/*   void resetRingIndex() {ring_index = -1;    } */

  void setLength(double l) { len = l; }
  double getLength() { return len; }
};

//------------------------------------------------
// Wrappers [Vertex, Face, Halfedge]
//------------------------------------------------
struct Wrappers_VFH:public CGAL::Polyhedron_items_3 {
  // wrap vertex
  template < class Refs, class Traits > struct Vertex_wrapper {
    typedef struct {
    public:
    //typedef typename Traits::Vector_3 Vector_3;
    //all types needed by the vertex...
      typedef typename Traits::FT FT;
      typedef typename Traits::Point_3 Point_3;
      typedef typename Traits::Vector_3 Vector_3;
    } FGeomTraits;
    typedef typename Traits::Point_3 Point_3;
    typedef My_vertex < Refs, CGAL::Tag_true, Point_3, FGeomTraits > Vertex;
  };

  // wrap face
  //NOTE: [HDS, Face] renamed [Polyhedron, Facet]
  template < class Refs, class Traits > struct Face_wrapper {
    //typedef typename Traits::Vector_3 Vector_3;
    //all types needed by the facet...
    typedef struct {
    public:
      //     typedef typename Traits::Point_3 Point_3;
      typedef typename Traits::Vector_3 Vector_3;
     } FGeomTraits;
    //custom type instantiated...
    typedef My_facet < Refs, CGAL::Tag_true, FGeomTraits > Face;
  };

  // wrap halfedge
  template < class Refs, class Traits > struct Halfedge_wrapper {
   typedef My_halfedge < Refs,
      CGAL::Tag_true,
      CGAL::Tag_true, CGAL::Tag_true > Halfedge;
  };
};

//------------------------------------------------
//PolyhedralSurf
//------------------------------------------------
typedef double                DFT;
typedef CGAL::Cartesian<DFT>  Data_Kernel;
typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH > Polyhedron;
typedef Data_Kernel::Vector_3 Vector_3;

class PolyhedralSurf:public Polyhedron {
public:
  PolyhedralSurf() {}
  
  static Vector_3 getHalfedge_vector(Halfedge * h);
  void compute_edges_length();
  double compute_mean_edges_length_around_vertex(Vertex * v);
  void compute_facets_normals();
  Vector_3 computeFacetsAverageUnitNormal(Vertex * v);
};

#endif
