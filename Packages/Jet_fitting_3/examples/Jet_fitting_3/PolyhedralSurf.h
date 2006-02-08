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
// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).
//----------------------------------------------------------------

template < class Refs, class Tag, class Pt, class FGeomTraits > 
class My_vertex:public CGAL::HalfedgeDS_vertex_base < Refs, Tag, Pt >
{
protected:
  typedef typename Refs::Vertex Vertex;
  typedef typename Refs::Halfedge Halfedge;
  typedef typename Refs::Face Facet;

  typedef typename FGeomTraits::Point_3 Point_3;
  typedef typename FGeomTraits::Vector_3 Vector_3;

public:
  char ring_tag;

  //this is for collecting i-th ring neighbours
protected:
  int ring_index;
public:
  void setRingIndex(int i) { ring_index = i; }
  int getRingIndex() { return ring_index; }
  void resetRingIndex() { ring_index = -1; }
  void setRingTag() {	ring_tag = 1; }
  char getRingTag() {	return ring_tag; }

  My_vertex(const Point_3 & pt):
    CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt),
    ring_tag(0), ring_index(-1)    {}
  My_vertex()    {}
};

//----------------------------------------------------------------
// Facet with normal and possibly more types. types are recovered
//from the FGeomTraits template arg
//----------------------------------------------------------------
template < class Refs, class Tag, class FGeomTraits >
class My_facet:public CGAL::HalfedgeDS_face_base < Refs, Tag >
{
  //protected:
public:
  typedef typename Refs::Vertex Vertex;
  typedef typename Refs::Vertex_handle Vertex_handle;
  typedef typename Refs::Halfedge Halfedge;
  typedef typename Refs::Halfedge_handle Halfedge_handle;

  typedef typename FGeomTraits::Vector_3 Vector_3;
  typedef typename FGeomTraits::Point_3 Point_3;

public:
  Vector_3 normal;

public:
  My_facet(): ring_index(-1) {}
  Vector_3 & getUnitNormal() { return normal; }
  void setNormal(Vector_3 & n) { normal = n; }

  //this is for collecting i-th ring neighbours
protected:
  int ring_index;
public:
  void setRingIndex(int i) { ring_index = i; }
  int getRingIndex() { return ring_index; }
  void resetRingIndex() { ring_index = -1; }
};

//----------------------------------------------------------------
// Halfedge
//----------------------------------------------------------------
template < class Refs, class Tprev, class Tvertex, class Tface, class FGeomTraits > 
class My_halfedge:public CGAL::HalfedgeDS_halfedge_base < Refs, Tprev, Tvertex, Tface >
{
protected:
  typedef typename FGeomTraits::Point_3 Point_3;
  typedef typename FGeomTraits::Vector_3 Vector_3;

protected:
  int ring_index;
  double len;
public:
  void setRingIndex(int i) {	ring_index = i;    }
  int getRingIndex() {return ring_index;    }
  void resetRingIndex() {ring_index = -1;    }
public:
  My_halfedge(): ring_index(-1) {}
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
      typedef typename Traits::Vector_3 Vector_3;
      typedef typename Traits::Plane_3 Plane_3;
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
      typedef typename Traits::Point_3 Point_3;
      typedef typename Traits::Vector_3 Vector_3;
      typedef typename Traits::Plane_3 Plane_3;
    } FGeomTraits;
    //custom type instantiated...
    typedef My_facet < Refs, CGAL::Tag_true, FGeomTraits > Face;
  };

  // wrap halfedge
  template < class Refs, class Traits > struct Halfedge_wrapper {
    typedef struct {
    public:
      typedef typename Traits::Point_3 Point_3;
      typedef typename Traits::Vector_3 Vector_3;
      typedef typename Traits::Plane_3 Plane_3;
    } FGeomTraits;
    typedef My_halfedge < Refs,
      CGAL::Tag_true,
      CGAL::Tag_true, CGAL::Tag_true, FGeomTraits > Halfedge;
  };
};

//------------------------------------------------
// Polyhedron
//------------------------------------------------
//using standard Cartesian kernel
typedef double                DFT;
typedef CGAL::Cartesian<DFT>  Data_Kernel;

typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH > Polyhedron;
typedef Polyhedron::Vertex Vertex;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Halfedge Halfedge;
typedef Polyhedron::Halfedge_handle Halfedge_handle;

typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::
Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
typedef Polyhedron::
Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Polyhedron::Facet_iterator Facet_iterator;

typedef Data_Kernel::Point_3 DPoint;
typedef Data_Kernel::Vector_3 Vector_3;
typedef Data_Kernel::Plane_3 Plane_3;
typedef Data_Kernel::Sphere_3 Sphere_3;

/////////////////////class PolyhedralSurf///////////////////////////
class PolyhedralSurf:public Polyhedron {
public:
  typedef Data_Kernel::Point_3 DPoint;
  typedef Data_Kernel::Vector_3 Vector_3;
  typedef Data_Kernel::Plane_3 Plane_3;

public:
  static Vector_3 getHalfedge_vector(Halfedge * h);

  PolyhedralSurf() {}
  void compute_edges_length();
  double compute_mean_edges_length_around_vertex(Vertex * v);
  void compute_facets_normals();
  Vector_3 computeFacetsAverageUnitNormal(Vertex * v);
};

#endif
