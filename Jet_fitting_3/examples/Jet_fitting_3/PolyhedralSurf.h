#ifndef _POLYHEDRALSURF_H_
#define _POLYHEDRALSURF_H_

#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <algorithm>
#include <vector>
#include <list>

#include <stdlib.h>
#include <stdio.h>



//----------------------------------------------------------------
// some functors 
//----------------------------------------------------------------

//the facet stores the normal
struct Facet_unit_normal {
  template < class Facet > 
  void operator() (Facet & f) 
    {
      typename Facet::Halfedge_handle h = f.halfedge();
      typename Facet::Vector_3 normal =
	CGAL::cross_product(h->vertex()->point() - 
			    h->opposite()->vertex()->point(),
			    h->next()->vertex()->point() -
			    h->opposite()->vertex()->point());
      f.setNormal( normal / CGAL::sqrt(normal * normal));
    }
};


//----------------------------------------------------------------
// A redefined items class for the Polyhedron_3 with a refined vertex
// class that contains a member for the normal vector and a refined
// facet with a normal vector instead of the plane equation (this is
// an alternative solution instead of using Polyhedron_traits_with_normals_3).
//----------------------------------------------------------------

template < class Refs, class Tag, class Pt, class FGeomTraits > 
class My_vertex:public CGAL::HalfedgeDS_vertex_base < Refs, Tag, Pt >
{
typedef typename FGeomTraits::Point_3 Point_3;

public:
 My_vertex(const Point_3 & pt):
   CGAL::HalfedgeDS_vertex_base < Refs, Tag, Point_3 > (pt){}
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
 //  int ring_index;
//   My_halfedge(): ring_index(-1) {}
//   void setRingIndex(int i) {	ring_index = i;    }
//   int getRingIndex() {return ring_index;    }
//   void resetRingIndex() {ring_index = -1;    }

template < class Refs, class Tprev, class Tvertex, class Tface>
class My_halfedge:public CGAL::HalfedgeDS_halfedge_base < Refs, Tprev, Tvertex, Tface >
{
public:
  double len;
public:
  My_halfedge(): len(-1) {}
  void set_length(double l) { len = l; }
  double get_length() const { return len; }
  double& get_length_ref()  { return len; }
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
template <class TPoly>
class THEdge_PM : 
  public boost::put_get_helper<double, THEdge_PM<TPoly> > //read_write
  //  public boost::put_get_helper<double&, THEdge_PM<TPoly> > //lvalue
{
public: 
  typedef boost::read_write_property_map_tag category;//read_write
  //typedef boost::lvalue_property_map_tag category;//lvalue

  typedef typename TPoly::Halfedge key_type;
  typedef typename TPoly::Traits::FT value_type;
  typedef typename TPoly::Traits::FT reference;//read_write
  //typedef typename TPoly::Traits::FT& reference;//lvalue
  
  THEdge_PM() {}
  reference operator[](key_type h) const {return h.len;}//get_length();}//read_write
  //reference operator[](key_type h)  {return h.len;}//lvalue
};


//use the std edge_weight_t tag...
template <class TPoly> 
THEdge_PM<TPoly> get(boost::edge_weight_t, TPoly& P) {return THEdge_PM<TPoly>();}


typedef double                DFT;
typedef CGAL::Cartesian<DFT>  Data_Kernel;
typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH > Polyhedron;
typedef Data_Kernel::Vector_3 Vector_3;

class PolyhedralSurf:public Polyhedron {
public:
  PolyhedralSurf() {}
  
  static Vector_3 getHalfedge_vector(Halfedge * h);
  
//   void compute_edges_length();
//   double compute_mean_edges_length_around_vertex(Vertex * v);
  
  void compute_facets_normals();
  Vector_3 computeFacetsAverageUnitNormal(Vertex * v);
};

#endif
