#ifndef CGAL_POLYHEDRALSURF_H_
#define CGAL_POLYHEDRALSURF_H_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/property_map.h>
#include <boost/graph/properties.hpp>


#include <algorithm>
#include <vector>
#include <list>

#include <cstdlib>
#include <cstdio>




//---------------------------------------------------------------- A
//redefined items class for the Polyhedron_3 with a refined vertex
//class that contains nothing more! (the _ring_tag is given by an
//externa std::map; a refined facet with a normal vector instead of
//the plane equation (this is an alternative solution instead of using
//Polyhedron_traits_with_normals_3). edges with the length
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
  //int ring_index;

public:
  const Vector_3& get_unit_normal() const { return normal; }
  Vector_3& get_unit_normal() { return normal; }

  //My_facet(): ring_index(-1) {}
  //void setNormal(Vector_3  n) { normal = n; }
//   //this is for collecting i-th ring neighbours
//   void setRingIndex(int i) { ring_index = i; }
//   int getRingIndex() { return ring_index; }
//   void resetRingIndex() { ring_index = -1; }
};


template <class TPoly>
class Facet_PM :
  public boost::put_get_helper<typename TPoly::Traits::Vector_3, Facet_PM<TPoly> >
{
public:

  //read_write
  typedef boost::read_write_property_map_tag category;
  typedef typename TPoly::Facet key_type;
  typedef typename TPoly::Traits::Vector_3 value_type;
  typedef typename TPoly::Traits::Vector_3& reference;

  Facet_PM() {}
  reference operator[](key_type f) const {return f.get_unit_normal();}
};

//XFC: we should have Facet instead of Vertex!
namespace boost{
  enum vertex_attribute_t        { vertex_attribute        = 1111 };
  //BOOST_INSTALL_PROPERTY(facet, attribute);
  BOOST_INSTALL_PROPERTY(vertex, attribute);

}

template <class TPoly>
Facet_PM<TPoly> get_fpm(boost::vertex_attribute_t, TPoly& ) {return Facet_PM<TPoly>();}



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
  double& get_length()  { return len; }
};

/*XFC: tentative ... failed so far...*/
//property map associated to the half edge
template <class TPoly>
class HEdge_PM :
  public boost::put_get_helper<typename TPoly::Traits::FT&, HEdge_PM<TPoly> >//double
{
public:
  //read_write or lvalue
  //typedef boost::read_write_property_map_tag category;
  typedef boost::lvalue_property_map_tag category;
  typedef typename TPoly::Halfedge key_type;
  typedef typename TPoly::Traits::FT value_type;
  typedef typename TPoly::Traits::FT& reference;

  HEdge_PM() {}
  reference operator[](key_type h) const {return h.len;}
};

//use the std edge_weight_t tag...
template <class TPoly>
HEdge_PM<TPoly> get_hepm(boost::edge_weight_t, TPoly& )
{return HEdge_PM<TPoly>();}



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
typedef CGAL::Simple_cartesian<DFT>  Data_Kernel;
typedef CGAL::Polyhedron_3 < Data_Kernel, Wrappers_VFH > Polyhedron;
typedef Data_Kernel::Vector_3 Vector_3;

class PolyhedralSurf:public Polyhedron {
public:
  PolyhedralSurf() {}
};

#endif
