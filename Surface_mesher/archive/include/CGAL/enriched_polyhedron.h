// -*- tab-width: 2 -*-

// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: Enriched_polyhedron                                           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _POLYGON_MESH_
#define _POLYGON_MESH_

// CGAL stuff
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_default.h>
#include <list>

// a refined facet with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
  // tag
  int m_tag;

  // normal
  Norm m_normal;

public:

  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory

  Enriched_facet()
  {
  }

  // tag
  const int& tag() const { return m_tag; }
  int& tag() { return m_tag; }
        void tag(const int& i) { m_tag = i; }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }
};

// a refined halfedge with a general tag and
// a binary tag to indicate whether it belongs
// to the control mesh or not
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Enriched_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:

  // general purpose tag
  int m_tag;

  // option for edge superimposing
  bool m_control_edge;
  bool m_sharp;

public:

  // life cycle
  Enriched_halfedge()
  {
    m_control_edge = true;
    m_sharp = false;
  }

  // tag
  const int& tag() const { return m_tag;  }
  int& tag() { return m_tag;  }
  void tag(const int& t)  { m_tag = t; }

  // control edge
  bool& control_edge()  { return m_control_edge; }
  const bool& control_edge()  const { return m_control_edge; }

  // sharp
  bool& sharp()  { return m_sharp; }
  const bool& sharp()  const { return m_sharp; }
};



// a refined vertex with a normal and a tag
template <class Refs, class T, class P, class Norm>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  // tag
  int m_tag;

  // normal
  Norm m_normal;

public:
  // life cycle
  Enriched_vertex()  {}
  // repeat mandatory constructors
  Enriched_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
  {
  }

  // normal
  typedef Norm Normal_3;
  Normal_3& normal() { return m_normal; }
  const Normal_3& normal() const { return m_normal; }

  // tag
  int& tag() {  return m_tag; }
  const int& tag() const {  return m_tag; }
  void tag(const int& t)  { m_tag = t; }
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
  template<class Refs, class Traits> struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_vertex<Refs,
    CGAL::Tag_true,
    Point,
    Normal> Vertex;
  };

  // wrap face
  template<class Refs, class Traits> struct Face_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_facet<Refs,
    CGAL::Tag_true,
    Point,
    Normal> Face;
  };

  // wrap halfedge
  template<class Refs, class Traits> struct Halfedge_wrapper
  {
    typedef typename Traits::Vector_3 Normal;
    typedef Enriched_halfedge<Refs,
    CGAL::Tag_true,
    CGAL::Tag_true,
    CGAL::Tag_true,
    Normal> Halfedge;
  };
};

static const double PI = 3.1415926535897932384626;

// compute facet normal
struct Facet_normal // (functor)
{
  template<class Facet> void operator()(Facet& f)
  {
    typename Facet::Normal_3 sum = CGAL::NULL_VECTOR;
    typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
    do
    {
      typename Facet::Normal_3 normal = CGAL::cross_product(h->next()->vertex()->point() - h->vertex()->point(), h->next()->next()->vertex()->point() - h->next()->vertex()->point());
      double sqnorm = normal * normal;
      if (sqnorm != 0)
        normal = normal / (float)std::sqrt(sqnorm);
      sum = sum + normal;
    } while (++h != f.facet_begin());
    float sqnorm = sum * sum;
    if (sqnorm != 0.0)
      f.normal() = sum / std::sqrt(sqnorm);
    else
    {
      f.normal() = CGAL::NULL_VECTOR;
      //      TRACE("degenerate face\n");
    }
  }
};

// compute vertex normal
struct Vertex_normal // (functor)
{
  template<class Vertex> void operator()(Vertex& v)
  {
    typename Vertex::Normal_3 normal = CGAL::NULL_VECTOR;
    typename Vertex::Halfedge_around_vertex_const_circulator pHalfedge =
        v.vertex_begin();
    typename Vertex::Halfedge_around_vertex_const_circulator begin =
        pHalfedge;
    CGAL_For_all(pHalfedge,begin)
    if(!pHalfedge->is_border())
    normal = normal + pHalfedge->facet()->normal();
    float sqnorm = normal * normal;
    if (sqnorm != 0.0f)
      v.normal() = normal / (float)std::sqrt(sqnorm);
    else
      v.normal() = CGAL::NULL_VECTOR;
  }
};

//*********************************************************
template <class kernel,
          class items,
          template < class T, class I, class A>
          class HDS = CGAL::HalfedgeDS_default >
class Enriched_polyhedron : public CGAL::Polyhedron_3<kernel,items,HDS>
{
public :
  typedef typename kernel::FT FT;
  typedef typename kernel::Point_3 Point;
  typedef typename kernel::Vector_3 Vector;
  typedef typename kernel::Iso_cuboid_3 Iso_cuboid;
  typedef CGAL::Polyhedron_3<kernel,items,HDS> Base;

  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Vertex_iterator Vertex_iterator;
  typedef typename Base::Halfedge_handle Halfedge_handle;
  typedef typename Base::Halfedge_iterator Halfedge_iterator;
  typedef typename Base::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
  typedef typename Base::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Base::Edge_iterator Edge_iterator;
  typedef typename Base::Facet_iterator Facet_iterator;
  typedef typename Base::Facet_handle Facet_handle;
  typedef typename Base::Facet Facet;

  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::edges_begin;
  using Base::edges_end;
  using Base::halfedges_begin;
  using Base::halfedges_end;
  using Base::facets_begin;
  using Base::facets_end;
  using Base::size_of_halfedges;
  using Base::size_of_facets;
  using Base::size_of_vertices;
private :

  // bounding box
  Iso_cuboid m_bbox;

  // type
  bool m_pure_quad;
  bool m_pure_triangle;

public :
  enum Vertex_type { SMOOTH,           // 0 sharp edge
                     DART,             // 1
                     CREASE_REGULAR,   // 2 - with two non-sharp edges on each side
                     CREASE_IRREGULAR, // 2
                     CORNER };         // 3 and more

public :

  // life cycle
  Enriched_polyhedron()
  {
    m_pure_quad = false;
    m_pure_triangle = false;
  }
  virtual ~Enriched_polyhedron()
  {
  }

  // type
  bool is_pure_triangle() { return m_pure_triangle; }
  bool is_pure_quad() { return m_pure_quad; }

  // normals (per facet, then per vertex)
  void compute_normals_per_facet()
  {
    std::for_each(facets_begin(),facets_end(),::Facet_normal());
  }
  void compute_normals_per_vertex()
  {
    std::for_each(vertices_begin(),vertices_end(),::Vertex_normal());
  }
  void compute_normals()
  {
    compute_normals_per_facet();
    compute_normals_per_vertex();
  }

  // bounding box
  Iso_cuboid& bbox() { return m_bbox; }
  const Iso_cuboid bbox() const { return m_bbox; }

  // compute bounding box
  void compute_bounding_box()
  {
    CGAL_assertion(size_of_vertices() != 0);

    FT xmin,xmax,ymin,ymax,zmin,zmax;
    Vertex_iterator pVertex = vertices_begin();
    xmin = xmax = pVertex->point().x();
    ymin = ymax = pVertex->point().y();
    zmin = zmax = pVertex->point().z();
    for(;
        pVertex !=  vertices_end();
        pVertex++)
    {
      const Point& p = pVertex->point();

      xmin =  (std::min)(xmin,p.x());
      ymin =  (std::min)(ymin,p.y());
      zmin =  (std::min)(zmin,p.z());

      xmax =  (std::max)(xmax,p.x());
      ymax =  (std::max)(ymax,p.y());
      zmax =  (std::max)(zmax,p.z());
    }
    m_bbox = Iso_cuboid(xmin,ymin,zmin,
                        xmax,ymax,zmax);
  }

  // bounding box
  FT xmin() { return m_bbox.xmin(); }
  FT xmax() { return m_bbox.xmax(); }
  FT ymin() { return m_bbox.ymin(); }
  FT ymax() { return m_bbox.ymax(); }
  FT zmin() { return m_bbox.zmin(); }
  FT zmax() { return m_bbox.zmax(); }

  Point center()
  {
    FT cx = (FT)0.5 * (xmin() + xmax());
    FT cy = (FT)0.5 * (ymin() + ymax());
    FT cz = (FT)0.5 * (zmin() + zmax());
    return Point(cx,cy,cz);
  }

  FT size()
  {
    FT dx = xmax() - xmin();
    FT dy = ymax() - ymin();
    FT dz = zmax() - zmin();
    return std::max(dx,std::max(dy,dz));
  }

  // degree of a face
  static unsigned int degree(Facet_handle pFace)
  {
    return CGAL::circulator_size(pFace->facet_begin());
  }

  // valence of a vertex
  static unsigned int valence(Vertex_handle pVertex)
  {
    return CGAL::circulator_size(pVertex->vertex_begin());
  }

  // check whether a vertex is on a boundary or not
  static bool is_border(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    if(pHalfEdge == NULL) // isolated vertex
      return true;
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->is_border())
        return true;
    return false;
  }

  // get any border halfedge attached to a vertex
  Halfedge_handle get_border_halfedge(Vertex_handle pVertex)
  {
    Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator d = pHalfEdge;
    CGAL_For_all(pHalfEdge,d)
      if(pHalfEdge->is_border())
        return pHalfEdge;
    return NULL;
  }

  // tag all vertices
  void tag_vertices(const int tag)
  {
    for(Vertex_iterator vit = vertices_begin(),
                        end = vertices_end();
        vit != end;
        ++vit)
      vit->tag(tag);
  }

  // tag all halfedges
  void tag_halfedges(const int tag)
  {
    for(Halfedge_iterator pHalfedge = halfedges_begin();
        pHalfedge != halfedges_end();
        pHalfedge++)
      pHalfedge->tag(tag);
  }

  // tag all facets
  void tag_facets(const int tag)
  {
    for(Facet_iterator pFacet = facets_begin();
        pFacet  != facets_end();
        pFacet++)
      pFacet->tag() = tag;
  }

  // set index for all vertices
  void set_index_vertices()
  {
    int index = 0;
    for(Vertex_iterator pVertex = vertices_begin();
        pVertex != vertices_end();
        pVertex++)
      pVertex->tag(index++);
  }

  // set index for all edges
  void set_index_edges()
  {
    int index = 0;
    for(Edge_iterator he = edges_begin();
        he != edges_end();
        he++)
    {
      he->tag(index);
      he->opposite()->tag(index);
      index++;
    }
  }

  // is pure degree ?
  bool is_pure_degree(unsigned int d)
  {
    for(Facet_iterator pFace  = facets_begin();
        pFace != facets_end();
        pFace++)
      if(degree(pFace) != d)
        return false;
    return true;
  }

  // compute type
  void compute_type()
  {
    m_pure_quad = is_pure_degree(4);
    m_pure_triangle = is_pure_degree(3);
  }

  // compute facet center
  void compute_facet_center(Facet_handle pFace,
                            Point& center)
  {
    Halfedge_around_facet_circulator pHalfEdge = pFace->facet_begin();
    Halfedge_around_facet_circulator end = pHalfEdge;
    Vector vec(0.0,0.0,0.0);
    int degree = 0;
    CGAL_For_all(pHalfEdge,end)
    {
      vec = vec + (pHalfEdge->vertex()->point()-CGAL::ORIGIN);
      degree++;
    }
    center = CGAL::ORIGIN + (vec/FT(degree));
  }

  bool is_sharp(Halfedge_handle he,
                const double angle_sharp)
  {
    Facet_handle f1 = he->facet();
    Facet_handle f2 = he->opposite()->facet();
    if(f1 == Facet_handle() || f2 == Facet_handle() )
      return false;
    const Vector& n1 = f1->normal();
    const Vector& n2 = f2->normal();
    if(angle_deg(n1,n2) > angle_sharp)
      return true;
    else
      return false;
  }

  Vertex_type type(Vertex_handle v)
  {
    unsigned int nb = nb_sharp_edges(v);
    switch(nb)
    {
      case 0:
        return SMOOTH;
      case 1:
        return DART;
      case 2: // crease vertex - may be regular or not
        return crease_type(v);
      default: // 3 and more
        return CORNER;
    }
  }

  // regular crease vertex must have valence 6, with
  // exactly two smooth edges on each side.
  Vertex_type crease_type(Vertex_handle v)
  {
    if(valence(v) != 6)
      return CREASE_IRREGULAR;

    // valence = 6 - let us check regularity

    // pick first sharp edge
    Halfedge_around_vertex_circulator he = v->vertex_begin();
    Halfedge_around_vertex_circulator end = he;
    CGAL_For_all(he,end)
      if(he->sharp())
        break;

    // next two must be smooth
    for(int i=0;i<2;i++)
      if(++he->sharp())
        return CREASE_IRREGULAR;

    // next one must be sharp
    if(!++he->sharp())
      return CREASE_IRREGULAR;

    // next two must be smooth
    for(int i=0;i<2;i++)
      if(++he->sharp())
        return CREASE_IRREGULAR;

    return CREASE_REGULAR;
  }

  // return true if succeds
  void incident_points_on_crease(Vertex_handle v,
                                 Point& a,
                                 Point& b)
  {
#ifdef _DEBUG
    Vertex_type vertex_type = type(v,angle_sharp);
    ASSERT(vertex_type == CREASE_IRREGULAR || vertex_type == CREASE_REGULAR);
#endif // _DEBUG

    // pick first sharp edge
    Halfedge_around_vertex_circulator he = v->vertex_begin();
    Halfedge_around_vertex_circulator end = he;
    CGAL_For_all(he,end)
      if(he->sharp())
        break;
    a = he->opposite()->vertex()->point();

    // pick next sharp edge
    he++;
    CGAL_For_all(he,end)
      if(he->sharp())
        break;
    b = he->opposite()->vertex()->point();
  }

  // return nb of sharp edges incident to v
  unsigned int nb_sharp_edges(Vertex_handle v) const
  {
    Halfedge_around_vertex_circulator he = v->vertex_begin();
    Halfedge_around_vertex_circulator end = he;
    unsigned int nb_sharp_edges = 0;
    CGAL_For_all(he,end)
      if(he->sharp())
        nb_sharp_edges++;
    return nb_sharp_edges;
  }

  // Angle between two vectors (in degrees)
  // we use this formula
  // uv = |u||v| cos(u,v)
  // u  ^ v  = w
  // |w| = |u||v| |sin(u,v)|
  //**************************************************
  static double angle_deg(const Vector &u,
                          const Vector &v)
  {
    static const double conv = 1.0/PI*180.0;
    return conv * angle_rad(u,v);
  }

  static FT len(const Vector &v)
  {
    return (FT)std::sqrt(CGAL_NTS to_double(v*v));
  }

  // Angle between two vectors (in rad)
  // uv = |u||v| cos(u,v)
  // u  ^ v  = w
  // |w| = |u||v| |sin(u,v)|
  //**************************************************
  static double angle_rad(const Vector &u,
                          const Vector &v)
  {
    // check
    double product = len(u)*len(v);
    if(product == 0)
      return 0.0;

    // cosine
    double dot = (u*v);
    double cosine = dot / product;

    // sine
    Vector w = CGAL::cross_product(u,v);
    double AbsSine = len(w) / product;

    if(cosine >= 0)
      return std::asin(fix_sine(AbsSine));
    else
      return PI-std::asin(fix_sine(AbsSine));
  }

  //**********************************************
  // fix sine
  //**********************************************
  static double fix_sine(double sine)
  {
    if(sine >= 1)
      return 1;
    else
      if(sine <= -1)
        return -1;
      else
        return sine;
  }

  unsigned int tag_sharp_edges(const double angle_sharp)
  {
    unsigned int nb = 0;
    for(Halfedge_iterator he = edges_begin();
        he != edges_end();
        he++)
    {
      const bool tag = is_sharp(he,angle_sharp);
      he->sharp() = tag;
      he->opposite()->sharp() = tag;
      nb += tag ? 1 : 0;
    }
    return nb;
  }


  // count #boundaries
  unsigned int nb_boundaries()
  {
    unsigned int nb = 0;
    tag_halfedges(0);
    for(Halfedge_iterator he = halfedges_begin();
        he != halfedges_end();
        he++)
    {
      if(he->is_border() && he->tag() == 0)
      {
        nb++;
        Halfedge_handle curr = he;
        do
        {
          curr  = curr->next();
          curr->tag(1);
        }
        while(curr != he);
      }
    }
    return nb;
  }

  // tag component
  void tag_component(Facet_handle pSeedFacet,
                     const int tag_free,
                     const int tag_done)
{
    pSeedFacet->tag() = tag_done;
    std::list<Facet_handle> facets;
    facets.push_front(pSeedFacet);
    while(!facets.empty())
    {
      Facet_handle pFacet = facets.front();
      facets.pop_front();
      pFacet->tag() = tag_done;
      Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
      Halfedge_around_facet_circulator end = pHalfedge;
      CGAL_For_all(pHalfedge,end)
      {
        Facet_handle pNFacet = pHalfedge->opposite()->facet();
        if(pNFacet != NULL && pNFacet->tag() == tag_free)
        {
          facets.push_front(pNFacet);
          pNFacet->tag() = tag_done;
        }
      }
    }
}

  // count #components
  unsigned int nb_components()
  {
    unsigned int nb = 0;
    tag_facets(0);
    for(Facet_iterator pFacet = facets_begin();
        pFacet != facets_end();
        pFacet++)
    {
      if(pFacet->tag() == 0)
      {
        nb++;
        tag_component(pFacet,0,1);
      }
    }
    return nb;
  }

  // compute the genus
  // V - E + F + B = 2 (C - G)
  // C -> #connected components
  // G : genus
  // B : #boundaries
  int genus()
  {
    int c = nb_components();
    int b = nb_boundaries();
    int v = size_of_vertices();
    int e = size_of_halfedges()/2;
    int f = size_of_facets();
    return genus(c,v,f,e,b);
  }
  int genus(int c,
            int v,
            int f,
            int e,
            int b)
  {
    return (2*c+e-b-f-v)/2;
  }

};


#endif // _POLYGON_MESH_
