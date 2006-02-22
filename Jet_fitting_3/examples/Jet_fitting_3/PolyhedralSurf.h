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

//???????????dont anderstand when and why typedef typename ??

//----------------------------------------------------------------
// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).
//----------------------------------------------------------------

template < class Refs, class Tag, class Pt, class FGeomTraits > 
class My_vertex:public CGAL::HalfedgeDS_vertex_base < Refs, Tag, Pt >
{
protected:
  typedef typename FGeomTraits::Point_3 Point_3;
  int ring_index;
public:
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
//----------------------------------------------------------------
template < class Refs, class Tag, class FGeomTraits >
class My_facet:public CGAL::HalfedgeDS_face_base < Refs, Tag >
{
public:
 typedef typename FGeomTraits::Vector_3 Vector_3;

protected:
  Vector_3 normal;
  int ring_index;

public:
  My_facet(): ring_index(-1) {}
  Vector_3 & getUnitNormal() { return normal; }
  void setNormal(Vector_3  n) { normal = n; }

  //this is for collecting i-th ring neighbours
  void setRingIndex(int i) { ring_index = i; }
  int getRingIndex() { return ring_index; }
  void resetRingIndex() { ring_index = -1; }
};

//----------------------------------------------------------------
// Halfedge
//----------------------------------------------------------------
template < class Refs, class Tprev, class Tvertex, class Tface>
class My_halfedge:public CGAL::HalfedgeDS_halfedge_base < Refs, Tprev, Tvertex, Tface >
{
protected:
  int ring_index;
  double len;
public:
  My_halfedge(): ring_index(-1) {}
  void setRingIndex(int i) {	ring_index = i;    }
  int getRingIndex() {return ring_index;    }
  void resetRingIndex() {ring_index = -1;    }
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
      typedef typename Traits::Point_3 Point_3;
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
      typedef typename Traits::Vector_3 Vector_3;
    } FGeomTraits;
    //custom type instantiated...
    typedef My_facet < Refs, CGAL::Tag_true, FGeomTraits > Face;
  };

  // wrap halfedge
  template < class Refs, class Traits > struct Halfedge_wrapper {
   typedef My_halfedge < Refs,
      CGAL::Tag_true,
      CGAL::Tag_true, CGAL::Tag_true>  Halfedge;
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
