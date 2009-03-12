///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: Enriched_polyhedron                                           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _POLYGON_MESH_
#define _POLYGON_MESH_

// CGAL
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <list>

// Forward declarations
struct Vertex_normal;
struct Facet_normal;


// a refined facet with a normal
template <class Refs, class Tag, class Point_3_, class Norm>
class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, Tag>
{
  // normal
  Norm m_normal;

public:
  // life cycle
  Enriched_facet()
  {
  }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }
};


// a refined halfedge with a general tag
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Enriched_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:

public:

  // life cycle
  Enriched_halfedge()
  {
  }
};


// a refined vertex with a normal, tag and camera
template <class Refs, class Tag, class Point_3_, class Norm>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, Tag, Point_3_>
{
  // normal
  Norm m_normal;

public:
  // life cycle
  Enriched_vertex()  {}
  // repeat mandatory constructors
  Enriched_vertex(const Point_3_& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, Tag, Point_3_>(pt)
  {
  }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }
};

// A redefined items class for the Polyhedron_3
// with a refined vertex class that contains a
// member for the normal vector and a refined
// facet with a normal vector instead of the
// plane equation (this is an alternative
// solution instead of using
// Polyhedron_traits_with_normals_3).

struct Enriched_items : public CGAL::Polyhedron_items_3
{
  // wrap vertex
  template <class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_vertex<Refs,
      CGAL::Tag_true,
      Point,
      Normal> Vertex;
  };

  // wrap face
  template <class Refs, class Traits>
  struct Face_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_facet<Refs,
      CGAL::Tag_true,
      Point,
      Normal> Face;
  };

  // wrap halfedge
  template <class Refs, class Traits>
  struct Halfedge_wrapper
  {
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_halfedge<Refs,
      CGAL::Tag_true,
      CGAL::Tag_true,
      CGAL::Tag_true,
      Normal> Halfedge;
  };
};


// Enriched polyhedron
template <class PolyhedronTraits_3,
          class PolyhedronItems_3 = Enriched_items>
class Enriched_polyhedron
  : public CGAL::Polyhedron_3<PolyhedronTraits_3, PolyhedronItems_3>
{
  // Private types
private:

  typedef CGAL::Polyhedron_3<PolyhedronTraits_3, PolyhedronItems_3> Base;

  // Public types
public:

  typedef typename PolyhedronTraits_3::FT FT;
  typedef typename PolyhedronTraits_3::Point_3 Point;

  // Repeat Polyhedron_3 public types
  typedef typename Base::HalfedgeDS               HalfedgeDS;
  typedef typename Base::Vertex                   Vertex;
  typedef typename Base::Halfedge                 Halfedge;
  typedef typename Base::Facet                    Facet;
  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Vertex_const_handle      Vertex_const_handle;
  typedef typename Base::Halfedge_handle          Halfedge_handle;
  typedef typename Base::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename Base::Facet_handle             Facet_handle;
  typedef typename Base::Facet_const_handle       Facet_const_handle;
  typedef typename Base::Vertex_iterator          Vertex_iterator;
  typedef typename Base::Vertex_const_iterator    Vertex_const_iterator;
  typedef typename Base::Halfedge_iterator        Halfedge_iterator;
  typedef typename Base::Halfedge_const_iterator  Halfedge_const_iterator;
  typedef typename Base::Facet_iterator           Facet_iterator;
  typedef typename Base::Facet_const_iterator     Facet_const_iterator;
  typedef typename Base::Point_iterator           Point_iterator;
  typedef typename Base::Point_const_iterator     Point_const_iterator;
  typedef typename Base::Halfedge_around_facet_circulator         Halfedge_around_facet_circulator;
  typedef typename Base::Halfedge_around_facet_const_circulator   Halfedge_around_facet_const_circulator;
  typedef typename Base::Halfedge_around_vertex_circulator        Halfedge_around_vertex_circulator;
  typedef typename Base::Halfedge_around_vertex_const_circulator  Halfedge_around_vertex_const_circulator;

public :

  // Default constructor, copy constructor and operator =() are fine.

  // Repeat Delaunay_triangulation_3 public methods
  Base::halfedges_begin;
  Base::halfedges_end;
  Base::facets_begin;
  Base::facets_end;
  Base::vertices_begin;
  Base::vertices_end;

  // Compute normals using mesh connectivity (per facet, then per vertex)
  void compute_normals_per_facet()
  {
    std::for_each(facets_begin(),facets_end(),Facet_normal());
  }
  void compute_normals_per_vertex()
  {
    std::for_each(vertices_begin(),vertices_end(),Vertex_normal());
  }
  void compute_normals()
  {
    compute_normals_per_facet();
    compute_normals_per_vertex();
  }
};

// compute facet normal (functor)
struct Facet_normal
{
  template <class Facet>
  void operator()(Facet& f)
  {
    typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
    typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
    do
    {
      typename Facet::Normal_3 normal = CGAL::cross_product(
        h->next()->vertex()->point() - h->vertex()->point(),
        h->next()->next()->vertex()->point() - h->next()->vertex()->point());
      double sqnorm = normal * normal;
      if(sqnorm != 0)
        normal = normal / (double)std::sqrt(sqnorm);
      sum = sum + normal;
    }
    while(++h != f.facet_begin());
    double sqnorm = sum * sum;
    if(sqnorm != 0.0)
      f.normal() = sum / std::sqrt(sqnorm);
    else
    {
      f.normal() = CGAL::NULL_VECTOR;
      std::cerr << "degenerate face\n";
    }
  }
};


// compute vertex normal (functor)
struct Vertex_normal
{
  template <class Vertex>
  void operator()(Vertex& v)
  {
    typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
    typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge = v.vertex_begin();
    typename Vertex::Halfedge_around_vertex_const_circulator begin = pHalfedge;
    CGAL_For_all(pHalfedge,begin)
      if(!pHalfedge->is_border())
        normal = normal + pHalfedge->facet()->normal();
    double sqnorm = normal * normal;
    if(sqnorm != 0.0f)
      v.normal() = normal / (float)std::sqrt(sqnorm);
    else
      v.normal() = CGAL::NULL_VECTOR;
  }
};


#endif // _POLYGON_MESH_
