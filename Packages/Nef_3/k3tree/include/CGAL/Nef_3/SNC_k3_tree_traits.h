#line 839 "k3_tree.nw"
#ifndef SNC_K3_TREE_TRAITS_H
#define SNC_K3_TREE_TRAITS_H

#include <CGAL/Nef_3/Bounding_box_3.h>
#include <list>

#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)

CGAL_BEGIN_NAMESPACE

template <class SNCstructure>
class Side_of_plane {
public:
  typedef typename SNCstructure::Vertex_handle Vertex_handle;
  typedef typename SNCstructure::Halfedge_handle Halfedge_handle;
  typedef typename SNCstructure::Halffacet_handle Halffacet_handle;
  typedef typename SNCstructure::Halffacet_triangle_handle Halffacet_triangle_handle;
  typedef typename SNCstructure::Object_handle Object_handle;

  typedef typename SNCstructure::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename SNCstructure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SNCstructure::SHalfedge_handle SHalfedge_handle;
  
  typedef typename SNCstructure::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  
  virtual Oriented_side operator()( const Plane_3& pl, Object_handle o) const;
  virtual Oriented_side operator()( const Plane_3& pl, Vertex_handle v) const;
  virtual Oriented_side operator()( const Plane_3& pl, Halfedge_handle e) const;
  virtual Oriented_side operator()( const Plane_3& pl, Halffacet_handle f) const;
  virtual Oriented_side operator()( const Plane_3& pl, Halffacet_triangle_handle f) const;

  typedef typename SNCstructure::SNC_decorator SNC_decorator;
  SNC_decorator D;
};

template <class SNCstructure>
class Objects_bbox {
public:
  typedef typename SNCstructure::Vertex_handle Vertex_handle;
  typedef typename SNCstructure::Halfedge_handle Halfedge_handle;
  typedef typename SNCstructure::Halffacet_handle Halffacet_handle;
  typedef typename SNCstructure::Halffacet_triangle_handle Halffacet_triangle_handle;
  typedef typename SNCstructure::Object_handle Object_handle;
  typedef typename SNCstructure::Object_list Object_list;

  typedef typename SNCstructure::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename SNCstructure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SNCstructure::SHalfedge_handle SHalfedge_handle;

  typedef typename SNCstructure::Kernel Kernel;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  typedef typename Kernel::RT RT;
  typedef typename Kernel::FT FT;
  typedef Bounding_box_3<FT> Bounding_box_3;
  
  virtual Bounding_box_3 operator()(const Object_list& o) const;
  virtual Bounding_box_3 operator()(Object_handle o) const;
  virtual Bounding_box_3 operator()(Vertex_handle v) const;
  virtual Bounding_box_3 operator()(Halfedge_handle e) const;
  virtual Bounding_box_3 operator()(Halffacet_handle f) const;
  virtual Bounding_box_3 operator()(Halffacet_triangle_handle f) const;

  typedef typename SNCstructure::SNC_decorator SNC_decorator;
  SNC_decorator D;
};

template <class SNCstructure>
class SNC_k3_tree_traits {
  typedef typename SNCstructure::Kernel Kernel;

public:
  typedef SNCstructure SNC_structure;
  typedef typename SNCstructure::Infi_box Infimaximal_box;
  typedef typename SNCstructure::Vertex_handle Vertex_handle;
  typedef typename SNCstructure::Halfedge_handle Halfedge_handle;
  typedef typename SNCstructure::Halffacet_handle Halffacet_handle;
  typedef typename SNCstructure::Halffacet_triangle_handle Halffacet_triangle_handle;
  typedef typename SNCstructure::SNC_decorator Explorer;

  typedef typename SNCstructure::Object_handle Object_handle;
  typedef typename SNCstructure::Object_list Object_list;

  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Direction_3 Direction_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  typedef typename Kernel::RT RT;
  typedef typename Kernel::FT FT;
  typedef Bounding_box_3<FT> Bounding_box_3;

  typedef typename Kernel::Intersect_3 Intersect;
  typedef Objects_bbox<SNCstructure> Objects_bbox;
  typedef Side_of_plane<SNCstructure> Side_of_plane;

  Intersect intersect_object() const {
    return Intersect();
  }

  Side_of_plane side_of_plane_object() const {
    return Side_of_plane();
  }

  Objects_bbox objects_bbox_object() const {
    return Objects_bbox();
  }
};

template <class SNCstructure>
Oriented_side 
Side_of_plane<SNCstructure>::operator()
  ( const Plane_3& pl, Object_handle o) const {
  Vertex_handle v;
  Halfedge_handle e;
  Halffacet_handle f;
  Halffacet_triangle_handle t;
  if( assign( v, o))
    return (*this)( pl, v);
  else if( assign( e, o))
    return (*this)( pl, e);
  else if( assign( f, o))
    return (*this)( pl, f);
  else if( assign( t, o))
    return (*this)( pl, t);
  else
    CGAL_assertion_msg( 0, "wrong handle");
  return Oriented_side(); // never reached
}

template <class SNCstructure>
Oriented_side 
Side_of_plane<SNCstructure>::operator()
( const Plane_3& pl, Vertex_handle v) const {
  return pl.oriented_side(D.point(v));
}

/* 
   An edge is considered intersecting a plane if its endpoints lie on the
   plane or if they lie on diferent sides.  Partial tangency is not considered
   as intersection, due the fact that a lower dimensional face (the vertex)
   should be already reported as an object intersecting the plane.
 */

template <class SNCstructure>
Oriented_side 
Side_of_plane<SNCstructure>::operator()
( const Plane_3& pl, Halfedge_handle e) const {
  Segment_3 s(D.segment(e));
  Oriented_side src_side = pl.oriented_side(s.source());
  Oriented_side tgt_side = pl.oriented_side(s.target());
  if( src_side == tgt_side)
    return src_side;
  if( src_side == ON_ORIENTED_BOUNDARY)
    return tgt_side;
  if( tgt_side == ON_ORIENTED_BOUNDARY)
    return src_side;
  return ON_ORIENTED_BOUNDARY;
}

template <typename SNCstructure>
Oriented_side
Side_of_plane<SNCstructure>::operator()
( const Plane_3& pl, Halffacet_triangle_handle t) const {
  bool on_positive_side = false, on_negative_side = false;
  Triangle_3 tr(t.get_triangle());
  for( int i = 0; i < 3; ++i) {
    Oriented_side side = pl.oriented_side(tr[i]);
    if( side == ON_POSITIVE_SIDE)
      on_positive_side = true;
    else if( side == ON_NEGATIVE_SIDE)
      on_negative_side = true;
  }
  if( on_positive_side && on_negative_side)
    return ON_ORIENTED_BOUNDARY;
  if( !on_positive_side && !on_negative_side)
    return ON_ORIENTED_BOUNDARY;
  if( on_positive_side) {
    CGAL_assertion( !on_negative_side);
    return ON_POSITIVE_SIDE;
  }
  CGAL_assertion( on_negative_side);
  return ON_NEGATIVE_SIDE;
}

/* 
   As for the edges, if a facet is tanget to the plane it is not considered as
   a interesection since lower dimensional faces, like the edges and vertices
   where the tangency occurrs, should be reported as the objects intersecting 
   the plane.
   So, an intersection is reported if all vertices of the facet lie on plane,
   for which it is only necessary to check three vertices, or if the facet 
   has vertices on both sides of the plane, so the intersection is known
   as far as two vertices located on different sides of the plane.
*/

template <class SNCstructure>
Oriented_side 
Side_of_plane<SNCstructure>::operator()
  ( const Plane_3& pl, Halffacet_handle f) const {
    CGAL_nef3_assertion( std::distance( f->facet_cycles_begin(), f->facet_cycles_end()) > 0);
  Halffacet_cycle_iterator fc(f->facet_cycles_begin());
  SHalfedge_handle e;
  CGAL_assertion( assign( e, fc));
  assign( e, fc);
  SHalfedge_around_facet_circulator sc(e), send(sc);
  //CGAL_assertion( iterator_distance( sc, send) >= 3); // TODO: facet with 2 vertices was found, is it possible?
  Oriented_side facet_side;
  do {
    facet_side = pl.oriented_side(D.point(D.vertex(sc)));
    ++sc;
  }
  while( facet_side == ON_ORIENTED_BOUNDARY && sc != send);
  if( facet_side == ON_ORIENTED_BOUNDARY)
    return ON_ORIENTED_BOUNDARY;
  CGAL_assertion( facet_side != ON_ORIENTED_BOUNDARY);
  while( sc != send) {
    Oriented_side point_side = pl.oriented_side(D.point(D.vertex(sc)));
    ++sc;
    if( point_side == ON_ORIENTED_BOUNDARY)
      continue;
    if( point_side != facet_side)
      return ON_ORIENTED_BOUNDARY;
  }
  return facet_side;
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT>
Objects_bbox<SNCstructure>::operator()
  ( const Object_list& O) const {
  CGAL_assertion( O.size() >= 0);
  Bounding_box_3 b(0,0,0,0,0,0);
  typename Object_list::const_iterator o;
  for( o = O.begin(); o != O.end(); ++o) {
    Vertex_handle v;
    if( assign( v, *o))
      b = b + (*this)(v);
  }
  return b;
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT> 
Objects_bbox<SNCstructure>::operator()
  (Object_handle o) const {
  Vertex_handle v;
  Halfedge_handle e;
  Halffacet_handle f;
  Halffacet_triangle_handle t;
  if( assign( v, o))
    return operator()(v);
  else if( assign( e, o))
    return operator()(e);
  else if( assign( f, o))
    return operator()(f);
  else if( assign( t, o))
    return operator()(t);
  else
    CGAL_assertion_msg( 0, "wrong handle");
  return Bounding_box_3(); // never reached
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT> 
Objects_bbox<SNCstructure>::operator()
  (Vertex_handle v) const {
  Point_3 p(D.point(v));
  return Bounding_box_3( p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT> 
Objects_bbox<SNCstructure>::operator()
  (Halfedge_handle e) const {
  return (operator()(D.vertex(e)) + operator()(D.vertex(D.twin(e))));
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT> 
Objects_bbox<SNCstructure>::operator()
  (Halffacet_triangle_handle t) const {
  Bounding_box_3 bbox(0,0,0,0,0,0);
  Triangle_3 tr(t.get_triangle());
  for( int i = 0; i < 3; ++i) {
    Point_3 p(tr[i]);
    bbox = bbox + Bounding_box_3( p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
  }
  return bbox;
}

template <class SNCstructure>
Bounding_box_3<typename SNCstructure::Kernel::FT> 
Objects_bbox<SNCstructure>::operator()
  (Halffacet_handle f) const { // TODO
  CGAL_assertion( f->facet_cycles_begin() != Halffacet_cycle_iterator());
  Halffacet_cycle_iterator fc(f->facet_cycles_begin());
  SHalfedge_handle e;
  CGAL_assertion( assign( e, fc));
  assign( e, fc);
  SHalfedge_around_facet_circulator sc(e), send(sc);
  CGAL_assertion( !is_empty_range( sc, send));
  Bounding_box_3 b(operator()(D.vertex(sc)));
  sc++;
  while( sc != send) {
    b = b + operator()(D.vertex(sc));
    sc++;
  }
  return b;
}

CGAL_END_NAMESPACE

#endif // SNC_K3_TREE_TRAITS_H

