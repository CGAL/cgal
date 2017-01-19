// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SNC_K3_TREE_TRAITS_H
#define CGAL_NEF_SNC_K3_TREE_TRAITS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_3/Bounding_box_3.h>
#include <CGAL/Lazy_kernel.h>
#include <list>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 503
#include <CGAL/Nef_2/debug.h>

namespace CGAL {


template <typename Kernel, typename Coordinate>
class ComparePoints {

  typedef typename Kernel::Point_3  Point_3;
 public:
  ComparePoints(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  CGAL::Comparison_result operator()(const Point_3 p1, const Point_3& p2) {
    switch(coord) {
    case 0: 
      CGAL_NEF_TRACEN("compare_x " << p1 << ", " << p2 << "=" << (int) CGAL::compare_x(p1, p2));
      return CGAL::compare_x(p1, p2);
    case 1: 
      CGAL_NEF_TRACEN("compare_y " << p1 << ", " << p2 << "=" << (int) CGAL::compare_y(p1, p2));
      return CGAL::compare_y(p1, p2);
    case 2: 
      CGAL_NEF_TRACEN("compare_z " << p1 << ", " << p2 << "=" << (int) CGAL::compare_z(p1, p2));
      return CGAL::compare_z(p1, p2);
    default: CGAL_error();
    }
    return CGAL::EQUAL;
  }
private:
  Coordinate coord;
};


template <typename Coordinate, typename EK>
class ComparePoints<CGAL::Lazy_kernel<EK>, Coordinate> {

  typedef CGAL::Lazy_kernel<EK>     Kernel;
  typedef typename Kernel::Point_3  Point_3;
 public:
  ComparePoints(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  CGAL::Comparison_result operator()( const Point_3 p1, const Point_3 p2) {
    switch(coord) {
    case 0: 
      if(CGAL::to_interval(p1.x()).second <
	 CGAL::to_interval(p2.x()).first)
	return CGAL::SMALLER;
      else if(CGAL::to_interval(p2.x()).second <
	      CGAL::to_interval(p1.x()).first)
	return CGAL::LARGER;
      return CGAL::EQUAL;
    case 1:  
      if(CGAL::to_interval(p1.y()).second <
	 CGAL::to_interval(p2.y()).first)
	return CGAL::SMALLER;
      else if(CGAL::to_interval(p2.y()).second <
	      CGAL::to_interval(p1.y()).first)
	return CGAL::LARGER;
      return CGAL::EQUAL;
    case 2: 
      if(CGAL::to_interval(p1.z()).second <
	 CGAL::to_interval(p2.z()).first)
	return CGAL::SMALLER;
      else if(CGAL::to_interval(p2.z()).second <
	      CGAL::to_interval(p1.z()).first)
	return CGAL::LARGER;
      return CGAL::EQUAL;
    default: CGAL_error();
    }
    return CGAL::EQUAL;
  }
private:
  Coordinate coord;
};

template <class SNC_decorator>
class Side_of_plane {
public:
  typedef typename SNC_decorator::SNC_structure SNC_structure;	
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;

#ifdef CGAL_NEF3_TRIANGULATE_FACETS
  typedef typename SNC_structure::Halffacet_triangle_handle 
                                  Halffacet_triangle_handle;
#endif
#ifdef CGAL_NEF3_FACET_WITH_BOX
  typedef typename SNC_structure::Partial_facet Partial_facet;
#endif
  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename Decorator_traits::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;
  
  typedef typename SNC_decorator::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Vector_3 Vector_3;  
  typedef typename Kernel::RT  RT;

  typedef ComparePoints<Kernel, int> ComparePoints_;
  
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  Side_of_plane(bool rc = false) : reference_counted(rc) {}
#else    
    Side_of_plane() {}
#endif


  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Object_handle o, Depth depth);
  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Vertex_handle v, Depth depth);
  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Halfedge_handle e, Depth depth);
  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Halffacet_handle f, Depth depth);
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Halffacet_triangle_handle f, Depth depth);
#endif
#ifdef CGAL_NEF3_FACET_WITH_BOX
  template<typename Depth> Oriented_side operator()
    ( const Point_3& pop, Partial_facet& f, Depth depth);
#endif
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  bool reference_counted;
#endif
  SNC_decorator D;
  Unique_hash_map<Vertex_handle, Oriented_side> OnSideMap;
  Unique_hash_map<const RT*, Oriented_side> OnSideMapRC;
};

template <class SNC_decorator>
class Objects_bbox {
public:
  typedef typename SNC_decorator::SNC_structure SNC_structure;	
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;

  typedef typename SNC_structure::Halffacet_triangle_handle Halffacet_triangle_handle;
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef std::vector<Object_handle> Object_list;

  typedef typename Decorator_traits::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;

  typedef typename SNC_decorator::Kernel Kernel;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  typedef typename Kernel::RT RT;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Kernel_tag Kernel_tag;
  typedef CGAL::Bounding_box_3<Tag_true, Kernel> 
    Bounding_box_3;
  //  typedef CGAL::Bounding_box_3
  //    <typename Is_extended_kernel<Kernel>::value_type, Kernel> 
  //    Bounding_box_3;
  
  Bounding_box_3 operator()( const Object_list& O) const {
    Bounding_box_3 b;
    typename Object_list::const_iterator o = O.begin();
    Vertex_handle v;
    while(o != O.end() && !CGAL::assign(v, *o)) ++o;
    if(o != O.end()) {
      FT q[3];
      q[0] = v->point().x();
      q[1] = v->point().y();
      q[2] = v->point().z();
      Bounding_box_3 b(q);
      for(++o; o != O.end(); ++o) {
	if( CGAL::assign( v, *o)) {
	  b.extend(v->point());
	}
      }
      return b;
    }
    FT q[3];
    q[0] = q[1] = q[2] = 0;
    return Bounding_box_3(q);
  }

  /*
  Bounding_box_3 operator()(Object_handle o) const {
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    if( CGAL::assign( v, o))
      return operator()(v);
    else if( CGAL::assign( e, o))
      return operator()(e);
    else if( CGAL::assign( f, o))
      return operator()(f);
    else {
      Halffacet_triangle_handle t;
      if( CGAL::assign( t, o))
	return operator()(t);
      else
	CGAL_error_msg( "wrong handle");
    }
    return Bounding_box_3(); // never reached
  }

  Bounding_box_3 operator()(Vertex_handle v) const {    
    Bounding_box_3 b;
    b.extend(v->point());
    return b;
  }

  Bounding_box_3 operator()(Halfedge_handle e) const {
    Bounding_box_3 b;
    b.extend(e->source()->point());
    b.extend(e->twin()->source()->point());
    return b;
  }

  Bounding_box_3 operator()(Halffacet_triangle_handle t) const {
    Bounding_box_3 bbox;
    Triangle_3 tr(t.get_triangle());
    for( int i = 0; i < 3; ++i)
      bbox.extend(tr[i]);
    return bbox;
  }

  Bounding_box_3 operator()(Halffacet_handle f) const {
    CGAL_assertion( f->facet_cycles_begin() != 
		    Halffacet_cycle_iterator());
    Halffacet_cycle_iterator fc(f->facet_cycles_begin());
    SHalfedge_handle e;
    CGAL_assertion(fc.is_shalfedge());
    e = SHalfedge_handle(fc);
    SHalfedge_around_facet_circulator sc(e), send(sc);
    CGAL_assertion( !is_empty_range( sc, send));
    Bounding_box_3 b;
    CGAL_For_all( sc, send)
      b.extend(sc->source()->source()->point());
    return b;
  }
  */
};

template <class Decorator>
class SNC_k3_tree_traits {

public:
  typedef Decorator SNC_decorator;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
  typedef typename SNC_structure::Kernel Kernel;

  typedef typename SNC_structure::Infi_box Infimaximal_box;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
  typedef typename SNC_structure::Halffacet_triangle_handle 
                                  Halffacet_triangle_handle;
#endif
#ifdef CGAL_NEF3_FACET_WITH_BOX
  typedef typename SNC_structure::Partial_facet Partial_facet;
#endif

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef std::vector<Object_handle> Object_list;

  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;

  typedef typename Kernel::RT RT;
  typedef typename Kernel::Kernel_tag Kernel_tag;
  //  typedef CGAL::Bounding_box_3
  //    <typename Is_extended_kernel<Kernel>::value_type, Kernel> 
  //    Bounding_box_3;
  typedef CGAL::Bounding_box_3<Tag_true, Kernel> 
    Bounding_box_3;

  typedef typename Kernel::Intersect_3 Intersect;
  typedef CGAL::Objects_bbox<SNC_decorator> Objects_bbox;
  typedef CGAL::Side_of_plane<SNC_decorator> Side_of_plane;

  Intersect intersect_object() const {
    return Intersect();
  }

  Objects_bbox objects_bbox_object() const {
    return Objects_bbox();
  }
};

template <class SNC_decorator>
template <typename Depth>
Oriented_side 
Side_of_plane<SNC_decorator>::operator()
  (const Point_3& pop, Object_handle o, Depth depth) {
  Vertex_handle v;
  Halfedge_handle e;
  Halffacet_handle f;
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
  Halffacet_triangle_handle t;
#endif
#ifdef CGAL_NEF3_FACET_WITH_BOX
  Partial_facet pf;
#endif
  if( CGAL::assign( v, o))
    return (*this)(pop, v, depth);
  else if( CGAL::assign( e, o))
    return (*this)(pop, e, depth);
  else if( CGAL::assign( f, o))
    return (*this)(pop, f, depth);
#ifdef CGAL_NEF3_FACET_WITH_BOX
  else if( CGAL::assign(pf, o))
    return (*this)(pop, pf, depth);
#endif
#ifdef CGAL_NEF3_TRIANGULATE_FACETS
  else if( CGAL::assign( t, o))
    return (*this)(pop, t, depth);
#endif
  else
    CGAL_error_msg( "wrong handle");
  
  return Oriented_side(); // never reached
}

template <class SNC_decorator>
template <typename Depth>
Oriented_side 
Side_of_plane<SNC_decorator>::operator()
( const Point_3& pop, Vertex_handle v, Depth depth) {
  Comparison_result cr;
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  if(reference_counted) {
    if(!OnSideMapRC.is_defined(&(v->point().hw())))
      switch(depth%3) {
      case 0: 
        cr = CGAL::compare_x(v->point(), pop);
        OnSideMapRC[&(v->point().hw())] = cr == LARGER ? ON_POSITIVE_SIDE :
                         cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      case 1:
        cr = CGAL::compare_y(v->point(), pop);
        OnSideMapRC[&(v->point().hw())] = cr == LARGER ? ON_POSITIVE_SIDE :
                         cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      case 2:
        cr = CGAL::compare_z(v->point(), pop);
        OnSideMapRC[&(v->point().hw())] = cr == LARGER ? ON_POSITIVE_SIDE :
                         cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      default: CGAL_error_msg( "wrong value");
      }
    return OnSideMapRC[&(v->point().hw())];
  } else {
#endif
    ComparePoints_ compare(depth%3);
    if(!OnSideMap.is_defined(v)) {
      cr = compare(v->point(), pop);
      OnSideMap[v] = cr == LARGER ? ON_POSITIVE_SIDE :
	cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
    }
    return OnSideMap[v];
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  }
#endif
  CGAL_error_msg( "should not be reached");
}

/* 
   An edge is considered intersecting a plane if its endpoints lie on the
   plane or if they lie on diferent sides.  Partial tangency is not considered
   as intersection, due the fact that a lower dimensional face (the vertex)
   should be already reported as an object intersecting the plane.
 */

template <class SNC_decorator>
template <typename Depth>
Oriented_side 
Side_of_plane<SNC_decorator>::operator()
( const Point_3& pop, Halfedge_handle e, Depth depth) {
  Vertex_handle v = e->source();
  Vertex_handle vt = e->twin()->source();
  /*
  Comparison_result cr;
  if(!OnSideMap.is_defined(v))
    switch(depth%3) {
    case 0: 
      cr = CGAL::compare_x(v->point(), pop);
      OnSideMap[v] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    case 1:
      cr = CGAL::compare_y(v->point(), pop);
      OnSideMap[v] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    case 2:
      cr = CGAL::compare_z(v->point(), pop);
      OnSideMap[v] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    default: CGAL_error_msg( "wrong value");
    }
  if(!OnSideMap.is_defined(vt))
    switch(depth%3) {
    case 0: 
      cr = CGAL::compare_x(vt->point(), pop);
      OnSideMap[vt] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    case 1:
      cr = CGAL::compare_y(vt->point(), pop);
      OnSideMap[vt] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    case 2:
      cr = CGAL::compare_z(vt->point(), pop);
      OnSideMap[vt] = cr == LARGER ? ON_POSITIVE_SIDE :
                       cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
      break;
    default: CGAL_error_msg( "wrong value");
    }
  Oriented_side src_side = OnSideMap[v];
  Oriented_side tgt_side = OnSideMap[vt];
*/
  Oriented_side src_side = (*this) (pop, v, depth);
  Oriented_side tgt_side = (*this) (pop, vt, depth);
  if( src_side == tgt_side)
    return src_side;
  if( src_side == ON_ORIENTED_BOUNDARY)
    return tgt_side;
  if( tgt_side == ON_ORIENTED_BOUNDARY)
    return src_side;
  return ON_ORIENTED_BOUNDARY;
}

#ifdef CGAL_NEF3_TRIANGULATE_FACETS
template <typename SNC_decorator>
template <typename Depth>
Oriented_side
Side_of_plane<SNC_decorator>::operator()
( const Point_3& pop, Halffacet_triangle_handle t, Depth depth) {
  bool on_positive_side = false, on_negative_side = false;
  Triangle_3 tr(t.get_triangle());
  for( int i = 0; i < 3; ++i) {
    Oriented_side side = ON_ORIENTED_BOUNDARY;
    Comparison_result cr;
    if(
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
       !reference_counted || 
#endif
       !OnSideMapRC.is_defined(&(tr[i].hw()))) {
      switch(depth%3) {
      case 0: 
        cr = CGAL::compare_x(tr[i], pop);
        side = cr == LARGER ? ON_POSITIVE_SIDE :
                 cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      case 1:
        cr = CGAL::compare_y(tr[i], pop);
        side = cr == LARGER ? ON_POSITIVE_SIDE :
                 cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      case 2:
        cr = CGAL::compare_z(tr[i], pop);
        side = cr == LARGER ? ON_POSITIVE_SIDE :
                 cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
        break;
      default: CGAL_error_msg( "wrong value");
      }
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
      if(reference_counted) 
        OnSideMapRC[&(tr[i].hw())] = side;
    } else if(reference_counted)
      side = OnSideMapRC[&(tr[i].hw())];
#endif
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
#endif

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

#ifdef CGAL_NEF3_FACET_WITH_BOX
template <class SNC_decorator>
template <typename Depth>
Oriented_side
Side_of_plane<SNC_decorator>::operator()
  (const Point_3& pop, Partial_facet& pf, Depth depth) {
  CGAL_error_msg( "not implemented yet");
  
  return ON_ORIENTED_BOUNDARY;
}
#endif

template <class SNC_decorator>
template <typename Depth>
Oriented_side 
Side_of_plane<SNC_decorator>::operator()
  (const Point_3& pop, Halffacet_handle f, Depth depth) {
    CGAL_assertion( std::distance( f->facet_cycles_begin(), f->facet_cycles_end()) > 0);
    /*
#ifdef CGAL_NEF3_FACET_WITH_BOX
    switch(depth%3) {
    case 0:
      if(f->b.min_coord(0) > pop.x())
	return ON_POSITIVE_SIDE;
      if(f->b.max_coord(0) < pop.x())
	return ON_NEGATIVE_SIDE;
      break;
    case 1:
      if(f->b.min_coord(1) > pop.y())
	return ON_POSITIVE_SIDE;
      if(f->b.max_coord(1) < pop.y())
	return ON_NEGATIVE_SIDE;
      break;
    case 2:
      if(f->b.min_coord(2) > pop.z())
	return ON_POSITIVE_SIDE;
      if(f->b.max_coord(2) < pop.z())
	return ON_NEGATIVE_SIDE;
      break;
    default: CGAL_error_msg( "wrong value");
    }    
    return ON_ORIENTED_BOUNDARY;
#else  
    */
  Halffacet_cycle_iterator fc(f->facet_cycles_begin());
  SHalfedge_handle e;
  CGAL_assertion(fc.is_shalfedge());
  e = SHalfedge_handle(fc);
  SHalfedge_around_facet_circulator sc(e), send(sc);
  //CGAL_assertion( iterator_distance( sc, send) >= 3); // TODO: facet with 2 vertices was found, is it possible?

  Oriented_side facet_side;
  Vertex_handle v;
  do {
    v = sc->source()->center_vertex();
    facet_side = (*this) (pop, v, depth);
    ++sc;
  }
  while( facet_side == ON_ORIENTED_BOUNDARY && sc != send);
  if( facet_side == ON_ORIENTED_BOUNDARY)
    return ON_ORIENTED_BOUNDARY;
  CGAL_assertion( facet_side != ON_ORIENTED_BOUNDARY);
  Oriented_side point_side;
  while( sc != send) {
    v = sc->source()->center_vertex();
    point_side = (*this) (pop, v, depth);
    ++sc;
    if( point_side == ON_ORIENTED_BOUNDARY)
      continue;
    if( point_side != facet_side)
      return ON_ORIENTED_BOUNDARY;
  }
  return facet_side;
  //#endif
}

} //namespace CGAL

#endif // CGAL_NEF_SNC_K3_TREE_TRAITS_H
