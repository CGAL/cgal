#ifndef CGAL_POLYHEDRALSURF_H_
#define CGAL_POLYHEDRALSURF_H_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <list>

#include <boost/foreach.hpp>

//----------------------------------------------------------------
// A redefined items class for the Polyhedron_3 with
// a refined facet with a normal vector
//---------------------------------------------------------------

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
public:
  My_facet() {}
  const Vector_3 & getUnitNormal() const { return normal; }
  void setNormal(Vector_3 n) { normal = n; }
};

//------------------------------------------------
// Wrappers [Vertex, Face, Halfedge]
//------------------------------------------------
struct Wrappers_VFH:public CGAL::Polyhedron_items_3 {
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
};

//------------------------------------------------
//PolyhedralSurf with facet normal operations
//------------------------------------------------
typedef double                FT;
typedef CGAL::Simple_cartesian<FT>  Kernel;
typedef CGAL::Polyhedron_3 < Kernel, Wrappers_VFH > Polyhedron;
typedef Kernel::Vector_3 Vector_3;

class PolyhedralSurf;

namespace boost {
  template<>
  struct graph_traits<PolyhedralSurf> : public boost::graph_traits<Polyhedron>
  {};  template<>
  struct graph_traits<PolyhedralSurf const> : public boost::graph_traits<Polyhedron>
  {};
}

class PolyhedralSurf : public Polyhedron {
public:
  typedef boost::graph_traits<PolyhedralSurf>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<PolyhedralSurf>::face_descriptor face_descriptor;
  typedef boost::graph_traits<PolyhedralSurf>::halfedge_descriptor halfedge_descriptor;

  PolyhedralSurf() {}
  void compute_facets_normals();
  const Vector_3 computeFacetsAverageUnitNormal(const vertex_descriptor v);

};

void PolyhedralSurf::compute_facets_normals()
{
  BOOST_FOREACH(face_descriptor f, faces(*this)){
    halfedge_descriptor h = halfedge(f,*this);
      Vector_3 normal =
	CGAL::cross_product(target(h,*this)->point() -
			    target(opposite(h,*this),*this)->point(),
			    target(next(h,*this),*this)->point() -
			    target(opposite(h,*this),*this)->point());
      f->setNormal( normal / CGAL::sqrt(normal * normal));
  }
}

const Vector_3 PolyhedralSurf::computeFacetsAverageUnitNormal(vertex_descriptor v)
{
  halfedge_descriptor h;
  face_descriptor f;
  Vector_3 sum(0., 0., 0.), n;

  CGAL::Halfedge_around_target_circulator<PolyhedralSurf> hedgeb(halfedge(v,*this),*this), hedgee = hedgeb;

  do
    {
      h = *hedgeb;
      if (is_border_edge(h,*this))
	{
	  hedgeb++;
	  continue;
	}

      f =  face(h,*this);
      n = f->getUnitNormal();
      sum = (sum + n);
      hedgeb++;
    }
  while (hedgeb != hedgee);
  sum = sum / std::sqrt(sum * sum);
  return sum;
}


#endif
