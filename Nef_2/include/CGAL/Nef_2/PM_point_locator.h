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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_PM_POINT_LOCATOR_H
#define CGAL_PM_POINT_LOCATOR_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Constrained_triang_traits.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Circulator_project.h>
#include <CGAL/Polygon_2_algorithms.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 17
#include <CGAL/Nef_2/debug.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif

#ifdef CGAL_USE_LEDA_LIBRARY
#include <CGAL/LEDA_basic.h>
# if __LEDA__ > 410 && __LEDA__ < 441
#  define CGAL_USING_PPL
#  include <CGAL/Nef_2/PM_persistent_PL.h>
# endif
#endif

namespace CGAL {

enum object_kind { VERTEX, EDGE_CROSSING, EDGE_COLLINEAR };

template < class Node, class Object>
struct Project_halfedge_point {
  typedef Node         argument_type;
  typedef Object       result_type;
  Object& operator()( Node& x) const   {
    return x.vertex()->point();
  }
  const Object& operator()( const Node& x) const   {
    return x.vertex()->point();
  }
};

/*{\Moptions print_title=yes }*/
/*{\Msubst
PM_decorator_#PMD
Geometry_#GEO
}*/
/*{\Manpage {PM_naive_point_locator}{PMD,GEO}
{Naive point location in plane maps}{PL}}*/
/*{\Mdefinition An instance |\Mvar| of data type |\Mname|
encapsulates naive point location queries within a plane map |P|.  The
two template parameters are specified via concepts. |PM_decorator_|
must be a model of the concept |PMDecorator| as described in the
appendix.  |Geometry_| must be a model of the concept
|AffineGeometryTraits_2| as described in the appendix. For a
specification of plane maps see also the concept of
|PMConstDecorator|.}*/

/*{\Mgeneralization PMD}*/

template <typename PM_decorator_, typename Geometry_>
class PM_naive_point_locator : public PM_decorator_ {
protected:
  typedef PM_decorator_ Base;
  typedef PM_naive_point_locator<PM_decorator_,Geometry_> Self;

  const Geometry_& K;
public:
  /*{\Mtypes 5}*/
  typedef PM_decorator_                 Decorator;
  /*{\Mtypemember equals |PM_decorator_|.}*/
  typedef typename Decorator::Plane_map Plane_map;
  /*{\Mtypemember the plane map type decorated by |Decorator|.}*/
  typedef typename Decorator::Mark      Mark;
  /*{\Mtypemember the attribute of all objects (vertices, edges,
  faces).}*/

  typedef Geometry_                       Geometry;
  /*{\Mtypemember equals |Geometry_|.}*/
  typedef typename Geometry_::Point_2     Point;
  /*{\Mtypemember the point type of the geometry kernel.\\
  \require |Geometry::Point_2| equals |Plane_map::Point|.}*/
  typedef typename Geometry_::Segment_2   Segment;
  /*{\Mtypemember the segment type of the geometry kernel.}*/
  typedef typename Geometry_::Direction_2 Direction;

  /*{\Mtext Local types are handles, iterators and circulators of the
  following kind: |Vertex_const_handle|, |Vertex_const_iterator|,
  |Halfedge_const_handle|, |Halfedge_const_iterator|, |Face_const_handle|,
  |Face_const_iterator|.}*/

  typedef CGAL::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying plane
  map. The kind of the object |(vertex, halfedge,face)| can be determined and
  the object assigned by the three functions:\\
  |bool assign(Vertex_const_handle& h, Object_handle o)|\\
  |bool assign(Halfedge_const_handle& h, Object_handle o)|\\
  |bool assign(Face_const_handle& h, Object_handle o)|\\ where each
  function returns |true| iff the assignment of |o| to |h| was valid.}*/

   typedef typename PM_decorator_::Vertex_handle Vertex_handle;
   typedef typename PM_decorator_::Halfedge_handle Halfedge_handle;
   typedef typename PM_decorator_::Face_handle Face_handle;
   typedef typename PM_decorator_::Vertex_const_handle Vertex_const_handle;
   typedef typename PM_decorator_::Halfedge_const_handle Halfedge_const_handle;
   typedef typename PM_decorator_::Face_const_handle Face_const_handle;
   typedef typename PM_decorator_::Vertex_iterator Vertex_iterator;
   typedef typename PM_decorator_::Halfedge_iterator Halfedge_iterator;
   typedef typename PM_decorator_::Face_iterator Face_iterator;
   typedef typename PM_decorator_::Vertex_const_iterator Vertex_const_iterator;
   typedef typename PM_decorator_::Halfedge_const_iterator Halfedge_const_iterator;
   typedef typename PM_decorator_::Face_const_iterator Face_const_iterator;
   typedef typename PM_decorator_::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
   typedef typename PM_decorator_::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
   typedef typename PM_decorator_::Halfedge_around_face_circulator Halfedge_around_face_circulator;
   typedef typename PM_decorator_::Halfedge_around_face_const_circulator Halfedge_around_face_const_circulator;

  using Base::clear;
  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::halfedges_begin;
  using Base::halfedges_end;
  using Base::faces_begin;
  using Base::faces_end;
  using Base::number_of_vertices;
  using Base::number_of_halfedges;
  using Base::number_of_faces;
  using Base::info;
  using Base::is_closed_at_source;
  using Base::source;
  using Base::target;
  using Base::cyclic_adj_succ;
  using Base::mark;
  using Base::twin;
  using Base::flip_diagonal;
  using Base::is_isolated;
  using Base::first_out_edge;
  using Base::next;
  using Base::previous;
  using Base::face;
  using Base::point;

  Halfedge_const_handle out_wedge(Vertex_const_handle v,
    const Direction& d, bool& collinear) const
  /*{\Xop returns a halfedge |e| bounding a wedge in between two
  neighbored edges in the adjacency list of |v| which contains |d|.
  If |d| extends along a edge then |e| is this edge. If |d| extends
  into the interior of such a wedge then |e| is the first edge hit
  when |d| is rotated clockwise. \precond |v| is not isolated.}*/
  { CGAL_NEF_TRACEN("out_wedge "<<PV(v));
    CGAL_assertion(!is_isolated(v));
    collinear=false;
    Point p = point(v);
    Halfedge_const_handle e_res = first_out_edge(v);
    Direction d_res = direction(e_res);
    Halfedge_around_vertex_const_circulator el(e_res),ee(el);
    CGAL_For_all(el,ee) {
      if ( K.strictly_ordered_ccw(d_res, direction(el), d) )
        e_res = el; d_res = direction(e_res);
    }
    CGAL_NEF_TRACEN("  determined "<<PE(e_res)<<" "<<d_res);
    if ( direction(cyclic_adj_succ(e_res)) == d ) {
      e_res = cyclic_adj_succ(e_res);
      collinear=true;
    }
    CGAL_NEF_TRACEN("  wedge = "<<PE(e_res)<<" "<<collinear);
    return e_res;
  }

  Segment segment(Halfedge_const_handle e) const
  { return K.construct_segment(point(source(e)), point(target(e))); }

  Direction direction(Halfedge_const_handle e) const
  { return K.construct_direction(point(source(e)),point(target(e))); }

  /*{\Mcreation 3}*/

  /*{\Moptions constref=yes}*/
  PM_naive_point_locator(const Plane_map& P, const Geometry& k = Geometry()) :
    Base(const_cast<Plane_map&>(P)), K(k) {}
  /*{\Mcreate constructs a point locator working on |P|.}*/
  /*{\Moptions constref=no}*/
  /*{\Moperations 2.5 0.5}*/

  const Mark& mark(Object_handle h) const
  /*{\Mop returns the mark associated to the object |h|.}*/
  { Vertex_const_handle v;
    Halfedge_const_handle e;
    Face_const_handle f;
    if ( assign(v,h) ) return mark(v);
    if ( assign(e,h) ) return mark(e);
    if ( assign(f,h) ) return mark(f);
    CGAL_error_msg("PM_point_locator::mark: Object_handle holds no object.");
    return mark(v); // never reached
  }


  Object_handle locate(const Segment& s) const
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of the underlying plane map |P| which contains the point |p =
  s.source()| in its relative interior. |s.target()| must be a point
  such that |s| intersects the $1$-skeleton of |P|.}*/
  { CGAL_NEF_TRACEN("locate naivly "<<s);
    if (this->number_of_vertices() == 0)
      CGAL_error_msg("PM_naive_point_locator: plane map is empty.");
    Point p = K.source(s);
    Vertex_const_iterator vit;
    for(vit = this->vertices_begin(); vit != this->vertices_end(); ++vit) {
      if ( p == point(vit) ) return make_object(vit);
    }
    Halfedge_const_iterator eit;
    for(eit = this->halfedges_begin(); eit != this->halfedges_end(); ++(++eit)) {
      // we only have to check each second halfedge
      if ( K.contains(segment(eit),p) )
        return make_object(eit);
    }
    Vertex_const_handle v_res;
    Halfedge_const_handle e_res;
    Segment ss = s; // we shorten the segment iteratively
    Direction dso = K.construct_direction(K.target(s),p), d_res;
    CGAL::Unique_hash_map<Halfedge_const_handle,bool> visited(false);
    for(vit = this->vertices_begin(); vit != this->vertices_end(); ++vit) {
      Point p_res, vp = point(vit);
      if ( K.contains(ss,vp) ) {
        CGAL_NEF_TRACEN(" location via vertex at "<<vp);
        ss = K.construct_segment(p,vp); // we shrink the segment
        if ( is_isolated(vit) ) {
          v_res = vit; e_res = Halfedge_const_handle();
        } else { // not isolated
          bool dummy;
          e_res = out_wedge(vit,dso,dummy);
          Halfedge_around_vertex_const_circulator el(e_res),ee(el);
          CGAL_For_all(el,ee)
            visited[el] = visited[twin(el)] = true;
          /* e_res is now the counterclockwise maximal halfedge out
             of v just before s */
          if ( K.orientation(p,vp,point(target(e_res))) < 0 ) // right turn
            e_res = previous(e_res);
          // correction to make e_res visible from p
          CGAL_NEF_TRACEN("  determined "<<PE(e_res));
        }
      }
    }

    for (eit = this->halfedges_begin(); eit != this->halfedges_end(); ++eit) {
      if ( visited[eit] ) continue;
      Point se = point(source(eit)),
            te = point(target(eit));
      int o1 = K.orientation(ss,se);
      int o2 = K.orientation(ss,te);
      if ( o1 == -o2 && // internal intersection
           K.orientation(se,te,K.source(ss)) !=
           K.orientation(se,te,K.target(ss)) ) {
          CGAL_NEF_TRACEN(" location via halfedge "<<segment(eit));
        Point p_res = K.intersection(s,segment(eit));
        ss = K.construct_segment(p,p_res);
        e_res = (o2 > 0 ? eit : twin(eit));
        // o2>0 => te left of s and se right of s => p left of e
        visited[eit] = visited[twin(eit)] = true;
        CGAL_NEF_TRACEN("  determined "<<PE(e_res)<<" "<<mark(e_res));
        CGAL_NEF_TRACEN("             "<<mark(face(e_res)));
      }
    }

    if ( e_res != Halfedge_const_handle() )
      return make_object((Face_const_handle)(face(e_res)));
    else
      return make_object((Face_const_handle)(face(v_res)));
  }


  template <typename Object_predicate>
  Object_handle ray_shoot(const Segment& s, const Object_predicate& M) const
  /*{\Mop returns an |Object_handle o| which can be converted to a
  |Vertex_const_handle|, |Halfedge_const_handle|, |Face_const_handle|
  |h| as described above.  The object predicate |M| has to have function
  operators \\
  |bool operator() (const Vertex_/Halfedge_/Face_const_handle&)|.\\
  The object returned is intersected by the segment |s| and has
  minimal distance to |s.source()| and |M(h)| holds on the converted
  object. The operation returns the null handle |NULL| if the ray shoot
  along |s| does not hit any object |h| of |P| with |M(h)|.}*/
  { CGAL_NEF_TRACEN("naive ray_shoot "<<s);
    CGAL_assertion( !K.is_degenerate(s) );
    Point p = K.source(s);
    Segment ss(s);
    Direction d = K.construct_direction(K.source(s),K.target(s));
    Object_handle h = locate(s);
    Vertex_const_handle v;
    Halfedge_const_handle e;
    Face_const_handle f;
    if ( ( assign(v,h) && M(v) ) ||
         ( assign(e,h) && M(e) ) ||
         ( assign(f,h) && M(f) ) ) return h;
    h = Object_handle();
    CGAL_NEF_TRACEN("not contained");
    for (v = this->vertices_begin(); v != this->vertices_end(); ++v) {
      Point pv = point(v);
      if ( !K.contains(ss,pv) ) continue;
      CGAL_NEF_TRACEN("candidate "<<pv);
      if ( M(v) ) {
        h = make_object(v);     // store vertex
        ss = K.construct_segment(p,pv); // shorten
        continue;
      }
      // now we know that v is not marked but on s
      bool collinear;
      Halfedge_const_handle e = out_wedge(v,d,collinear);
      if ( collinear ) {
        if ( M(e) ) {
          h = make_object(e);
          ss = K.construct_segment(p,pv);
        }
        continue;
      }
      if ( M(face(e)) ) {
        h = make_object(face(e));
        ss = K.construct_segment(p,pv);
      }
    } // all vertices

    Halfedge_const_iterator e_res;
    for(e = this->halfedges_begin(); e != this->halfedges_end(); ++(++e)) {
      Segment es = segment(e);
      int o1 = K.orientation(ss,K.source(es));
      int o2 = K.orientation(ss,K.target(es));
      if ( o1 == -o2 && o1 != 0 &&
           K.orientation(es, K.source(ss)) ==
           - K.orientation(es, K.target(ss)) ) {
        // internal intersection
        CGAL_NEF_TRACEN("candidate "<<es);
        Point p_res = K.intersection(s,es);
        e_res = (o2 > 0 ? e : twin(e));
        // o2 > 0 => te left of s and se right of s => p left of e
        if ( M(e_res) ) {
          h = make_object(e_res);
          ss = K.construct_segment(p,p_res);
        } else if ( M(face(twin(e_res))) ) {
          h = make_object(face(twin(e_res)));
          ss = K.construct_segment(p,p_res);
        }
      }
    }

    return h;
  }


  // C++ is really friendly:
  #define USECMARK(t) const Mark& mark(t h) const { return Base::mark(h); }
  #define USEMARK(t)  Mark& mark(t h) const { return Base::mark(h); }
  USEMARK(Vertex_handle)
  USEMARK(Halfedge_handle)
  USEMARK(Face_handle)
  USECMARK(Vertex_const_handle)
  USECMARK(Halfedge_const_handle)
  USECMARK(Face_const_handle)
  #undef USEMARK
  #undef USECMARK
  /*{\Mimplementation Naive query operations are realized by checking
  the intersection points of the $1$-skeleton of the plane map |P| with
  the query segments $s$. This method takes time linear in the size $n$
  of the underlying plane map without any preprocessing.}*/
}; // PM_naive_point_locator<PM_decorator_,Geometry_>


/*{\Moptions print_title=yes }*/
/*{\Msubst
PM_decorator_#PMD
Geometry_#GEO
}*/
/*{\Manpage {PM_point_locator}{PMD,GEO}
{Point location in plane maps via LMWT}{PL}}*/
/*{\Mdefinition An instance |\Mvar| of data type |\Mname|
encapsulates point location queries within a plane map |P|. The two
template parameters are specified via concepts. |PMD| must be a model
of the concept |PMDecorator| as described in the appendix.  |GEO| must
be a model of the concept |AffineGeometryTraits_2| as described in the
appendix. For a specification of plane maps see also the concept of
|PMConstDecorator|.}*/

/*{\Mgeneralization PMD^#PM_naive_point_locator<PMD,GEO>}*/

template <typename PM_decorator_, typename Geometry_>
class PM_point_locator : public
  PM_naive_point_locator<PM_decorator_,Geometry_> {
protected:
  typedef PM_naive_point_locator<PM_decorator_,Geometry_> Base;
  typedef PM_point_locator<PM_decorator_,Geometry_> Self;
  Base CT;
  #ifdef CGAL_USING_PPL
  typedef PM_persistent_PL_traits<Base>  PMPPLT;
  typedef PointLocator<PMPPLT>           PMPP_locator;
  PMPP_locator* pPPL;
  #define LOCATE_IN_TRIANGULATION pPPL->locate_down
  #else
  #define LOCATE_IN_TRIANGULATION walk_in_triangulation
  #endif

public:

  typedef typename Base::Decorator Decorator;
  typedef typename Base::Plane_map Plane_map;
  typedef typename Base::Mark Mark;
  typedef typename Base::Geometry Geometry;
  typedef typename Base::Point Point;
  typedef typename Base::Segment Segment;
  typedef typename Base::Direction Direction;
  typedef typename Base::Object_handle Object_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Halfedge_handle Halfedge_handle;
  typedef typename Base::Face_handle Face_handle;
  typedef typename Base::Vertex_const_handle Vertex_const_handle;
  typedef typename Base::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Base::Face_const_handle Face_const_handle;
  typedef typename Base::Vertex_iterator Vertex_iterator;
  typedef typename Base::Halfedge_iterator Halfedge_iterator;
  typedef typename Base::Face_iterator Face_iterator;
  typedef typename Base::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Base::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Base::Face_const_iterator Face_const_iterator;
  typedef typename Base::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Base::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;
  typedef typename Base::Halfedge_around_face_circulator Halfedge_around_face_circulator;
  typedef typename Base::Halfedge_around_face_const_circulator Halfedge_around_face_const_circulator;

  using Base::K;
  using Base::number_of_vertices;
  using Base::faces_begin;
  using Base::info;
  using Base::flip_diagonal;
  using Base::twin; 
  using Base::next; 
  using Base::previous; 
  using Base::source; 
  using Base::target; 
  using Base::point; 
  using Base::segment; 
  using Base::face;

  /*{\Mtypes 2}*/
  /*{\Mtext All local types of |PM_naive_point_locator| are inherited.}*/
  typedef std::pair<Vertex_const_handle,Face_const_handle>   VF_pair;
  typedef std::pair<Halfedge_const_handle,Face_const_handle> EF_pair;

  struct CT_link_to_original : Decorator { // CT decorator
    const Decorator& Po;
    using Decorator::info;

    CT_link_to_original(const Decorator& P, const Decorator& Poi)
      : Decorator(P), Po(Poi) {}

    void operator()(Vertex_handle vn, Vertex_const_handle vo) const
    { Face_const_handle f;
      if ( Po.is_isolated(vo) ) f = Po.face(vo);
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<VF_pair>::create(info(vn));
      geninfo<VF_pair>::access(info(vn)) = VF_pair(vo,f);
      #else
      info(vn) = VF_pair(vo,f);      
      #endif
      CGAL_NEF_TRACEN("linking to org "<<PV(vn));
    }

    void operator()(Halfedge_handle hn, Halfedge_const_handle ho) const
    { 
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<EF_pair>::create(info(hn));
      geninfo<EF_pair>::access(info(hn)) = EF_pair(ho,Po.face(ho));
      #else
      info(hn) = EF_pair(ho,Po.face(ho));      
      #endif
      CGAL_NEF_TRACEN("linking to org "<<PE(hn));
    }
  };

protected:
  Vertex_const_handle input_vertex(Vertex_const_handle v) const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<VF_pair>::const_access(CT.info(v)).first; 
    #else
    return 
      boost::any_cast<VF_pair>(CT.info(v)).first; 
    #endif
  }

  Halfedge_const_handle input_halfedge(Halfedge_const_handle e) const
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<EF_pair>::const_access(CT.info(e)).first; 
    #else
    return 
      boost::any_cast<EF_pair>(CT.info(e)).first; 
    #endif
  }

  Face_const_handle input_face(Halfedge_const_handle e) const
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    return geninfo<EF_pair>::const_access(CT.info(e)).second;
    #else
    return 
      boost::any_cast<EF_pair>(CT.info(e)).second;
    #endif
  }


  Object_handle input_object(Vertex_const_handle v) const
  { return make_object(input_vertex(v)); }

  Object_handle input_object(Halfedge_const_handle e) const
  { Halfedge_const_handle e_org = input_halfedge(e);
    if ( e_org != Halfedge_const_handle() )
      return make_object( e_org );
    // now e_org is not existing
    return make_object( input_face(e) );
  }

  /*{\Mimplementation
  The efficiency of this point location module is mostly based on
  heuristics. Therefore worst case bounds are not very expressive. The
  query operations take up to linear time for subsequent query
  operations though they are better in practise. They trigger a one-time
  initialization which needs worst case $O(n^2)$ time though runtime
  tests often show subquadratic results. The necessary space for the
  query structure is subsumed in the storage space $O(n)$ of the input
  plane map. The query times are configuration dependent. If LEDA is
  present then point location is done via the slap method based on
  persistent dictionaries.  Then $T_{pl}(n) = O( \log(n) )$. If CGAL is
  not configured to use LEDA then point location is done via a segment
  walk in the underlying convex subdivision of $P$. In this case
  $T_{pl}(n)$ is the number of triangles crossed by a walk from the
  boundary of the structure to the query point. The time for the ray
  shooting operation $T_{rs}(n)$ is the time for the point location
  $T_{pl}(n)$ plus the time for the walk in the triangulation that is
  superimposed to the plane map. Let's consider the plane map edges as
  obstacles and the additional triangulation edges as non-obstacle
  edges. Let's call the sum of the lengths of all edges of the
  triangulation its weight. If the calculated triangulation
  approximates\cgalFootnote{The calculation of general
  minimum-weight-triangulations is conjectured to be NP-complete and
  locally-minimum-weight-triangulations that we use are considered good
  approximations.} the minimum weight triangulation of the obstacle set
  then the stepping quotient\cgalFootnote {The number of non-obstacle edges
  crossed until an obstacle edge is hit.} for a random direction of the
  ray shot is expected to be $O( \sqrt{n} )$.}*/


  struct CT_new_edge : Decorator {
    const Decorator& _DP;
    using Decorator::mark;
    using Decorator::previous;
    using Decorator::is_closed_at_source;
    using Decorator::info;
    using Decorator::source;
    using Decorator::twin;

    CT_new_edge(const Decorator& CT, const Decorator& DP) :
      Decorator(CT), _DP(DP) {}
    void operator()(Halfedge_handle& e) const
    { Halfedge_handle e_from = previous(e);
      Face_const_handle f;
      if ( is_closed_at_source(e) ) // source(e) was isolated before
        #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
        f = geninfo<VF_pair>::access(info(source(e))).second;
      else
        f = geninfo<EF_pair>::access(info(e_from)).second;      
        #else
        f = 
          boost::any_cast<VF_pair>(info(source(e))).second;
      else
        f = 
          boost::any_cast<EF_pair>(info(e_from)).second;              
        #endif
      mark(e) = _DP.mark(f);
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<EF_pair>::create(info(e));
      geninfo<EF_pair>::create(info(twin(e)));

      geninfo<EF_pair>::access(info(e)).first =
      geninfo<EF_pair>::access(info(twin(e))).first =
        Halfedge_const_handle();

      geninfo<EF_pair>::access(info(e)).second =
      geninfo<EF_pair>::access(info(twin(e))).second = f;
      #else
      info(e)=EF_pair(Halfedge_const_handle(),f);
      info(twin(e))=EF_pair(Halfedge_const_handle(),f);
      #endif
      CGAL_NEF_TRACEN("CT_new_edge "<<PE(e));
    }
  };

  void triangulate_CT() const
  {
    CGAL_NEF_TRACEN("triangulate_CT");
    typedef CGAL::Constrained_triang_traits<
      Decorator,Geometry,CT_new_edge> NCTT;
    typedef CGAL::generic_sweep<NCTT> Constrained_triang_sweep;
    CT_new_edge NE(CT,*this);
    Constrained_triang_sweep T(NE,CT.plane_map(),this->K); T.sweep();
  }

  void minimize_weight_CT() const
  { CGAL_NEF_TRACEN("minimize_weight_CT");
    if ( this->number_of_vertices() < 2 ) return;
    std::list<Halfedge_handle> S;
    /* We maintain a stack |S| of edges containing diagonals
       which might have to be flipped. */
    int flip_count = 0;
    Halfedge_iterator e;
    for (e = CT.halfedges_begin(); e != CT.halfedges_end(); ++(++e)) {
      Halfedge_const_handle e_org = input_halfedge(e);
      if ( e_org != Halfedge_const_handle() )
        continue;
      S.push_back(e);
    }

    while ( !S.empty() ) {
      Halfedge_handle e = S.front(); S.pop_front();
      Halfedge_handle r = twin(e);
      Halfedge_const_handle e_org = input_halfedge(e);
      if ( e_org != Halfedge_const_handle() )
        continue;
      Halfedge_handle e1 = next(r);
      Halfedge_handle e3 = next(e);
      // e1,e3: edges of quadrilateral with diagonal e

      Point a = point(source(e1));
      Point b = point(target(e1));
      Point c = point(source(e3));
      Point d = point(target(e3));

      if (! (this->K.orientation(b,d,a) > 0 && // left_turn
             this->K.orientation(b,d,c) < 0) ) // right_turn
        continue;

      if ( this->K.first_pair_closer_than_second(b,d,a,c) ) { // flip
        CGAL_NEF_TRACEN("flipping diagonal of quadilateral"<<a<<b<<c<<d);
        Halfedge_handle e2 = next(e1);
        Halfedge_handle e4 = next(e3);
        S.push_back(e1);
        S.push_back(e2);
        S.push_back(e3);
        S.push_back(e4);
        flip_diagonal(e);
        flip_count++;
      }


    }
    CGAL_NEF_TRACEN("  flipped "<<flip_count);
  }

public:
  /*{\Mcreation 3}*/

  PM_point_locator() {
    #ifdef CGAL_USING_PPL
    pPPL = 0;
    #endif

  }

  /*{\Moptions constref=yes}*/
  PM_point_locator(const Plane_map& P, const Geometry& k = Geometry());
  /*{\Mcreate constructs a point locator working on |P|.}*/

  ~PM_point_locator();

  /*{\Moperations 2.5 0.5}*/

  const Decorator& triangulation() const { return CT; }
  /*{\Mop access to the constrained triangulation structure that
  is superimposed to |P|.}*/
  /*{\Moptions constref=no}*/


  Object_handle locate(const Point& p) const
  /*{\Mop returns a generic handle |h| to an object (vertex, halfedge,
  face) of |P| which contains the point |p| in its relative
  interior.}*/
  {
    Object_handle h = LOCATE_IN_TRIANGULATION(p);
    Vertex_const_handle v_triang;
    if ( assign(v_triang,h) ) {
      return input_object(v_triang);
    }
    Halfedge_const_handle e_triang;
    if ( assign(e_triang,h) ) {
      Halfedge_const_handle e = input_halfedge(e_triang);
      if ( e == Halfedge_const_handle() ) // inserted during triangulation
        return make_object(input_face(e_triang));
      int orientation_ = this->K.orientation(segment(e),p);
      if ( orientation_ == 0 ) return make_object(e);
      if ( orientation_ < 0 )  return make_object(face(twin(e)));
      if ( orientation_ > 0 )  return make_object(face(e));
    }
    CGAL_assertion(!check_tag(typename Is_extended_kernel<Geometry>::value_type()));
    return make_object(Face_const_handle(faces_begin()));
    //    CGAL_error(); return h; // compiler warning
  }

  bool ray_shoot_from_outer_facet(Segment& , object_kind& ,
				  Vertex_const_handle &,
				  Halfedge_const_handle& ,
				  const Tag_true& ) const {
    return false;
  }

  bool ray_shoot_from_outer_facet(Segment& s, object_kind& current,
				  Vertex_const_handle &v,
				  Halfedge_const_handle& e,
				  const Tag_false& ) const {
    CGAL_NEF_TRACEN("target on outer facet");
    Point p = this->K.source(s);
    Vertex_const_handle v1 = CT.vertices_begin();
    Halfedge_const_handle e1 = CT.twin(CT.first_out_edge(v1));
    Halfedge_around_face_const_circulator circ(e1), end(circ);
    Point i;
    Segment seg;
    bool found = false;
    CGAL_For_all(circ, end) {
      //	std::cerr << s << std::endl;
      //	std::cerr << point(source(circ)) << "->" << point(target(circ)) << std::endl;
      Object o = intersection(s, Segment(point(source(circ)),
					 point(target(circ))));

      if(assign(i,o)) {
	CGAL_NEF_TRACEN("intersection in point " << i);
	found = true;
	s = Segment(p,i);
	if(i == point(source(circ))) {
	  current = VERTEX;
	  v = source(circ);
	} else if(i == point(target(circ))) {
	  current = VERTEX;
	  v = target(circ);
	} else {
	  current = EDGE_CROSSING;
	  e = circ;
	}
      } else if(assign(seg,o)) {
	found = true;
	CGAL_NEF_TRACEN("overlap of segments");
	current = EDGE_COLLINEAR;
	e = circ;
      }
    }
    return found;
  }

  template <typename Object_predicate>
  Object_handle ray_shoot(const Segment& ss, const Object_predicate& M) const
  /*{\Mop returns an |Object_handle o| which can be converted to a
  |Vertex_const_handle|, |Halfedge_const_handle|, |Face_const_handle|
  |h| as described above.  The object predicate |M| has to have
  function operators\\
  |bool operator() (const Vertex_/ Halfedge_/Face_const_handle&) const|.\\
  The object returned is intersected by the segment |s| and has minimal
  distance to |s.source()| and |M(h)| holds on the converted object. The
  operation returns the null handle |NULL| if the ray shoot along |s|
  does not hit any object |h| of |P| with |M(h)|.}*/
  { Segment s(ss);
    CGAL_NEF_TRACEN("ray_shoot "<<s);
    CGAL_assertion( !this->K.is_degenerate(s) );
    Point p = this->K.source(s);
    Direction d = this->K.construct_direction(p,s.target());
    Vertex_const_handle v;
    Halfedge_const_handle e;
    object_kind current;
    Object_handle h = LOCATE_IN_TRIANGULATION(p);
    if ( assign(v,h) ) {
      CGAL_NEF_TRACEN("located vertex "<<PV(v));
      current = VERTEX;
    } else if ( assign(e,h) ) {
      CGAL_NEF_TRACEN("located edge "<<PE(e));
      int orientation_ = this->K.orientation( segment(e), p);
      if ( orientation_ == 0 ) { // p on segment
        CGAL_NEF_TRACEN("on edge "<<PE(e));
        if ( d == CT.direction(e) )
        { current = EDGE_COLLINEAR; }
        else if ( d == CT.direction(CT.twin(e)) )
        { e = CT.twin(e); current = EDGE_COLLINEAR; }
        else { // crossing
          current = EDGE_CROSSING;
          if ( !(this->K.orientation(CT.segment(e),s.target())>0) ) // not left_turn
            e = CT.twin(e);
        }

      } else { // p not on segment, thus in triangle
        if ( orientation_ < 0  ) e = CT.twin(e);
        // now p left of e
        CGAL_NEF_TRACEN("in face at "<<PE(e));
        if ( M(input_face(e)) ) // face mark
          return make_object(input_face(e));

        Point p1 = CT.point(CT.source(e)),
              p2 = CT.point(CT.target(e)),
              p3 = CT.point(CT.target(next(e)));
        int or1 = this->K.orientation(p,s.target(),p1);
        int or2 = this->K.orientation(p,s.target(),p2);
        int or3 = this->K.orientation(p,s.target(),p3);
        if ( or1 == 0 && !this->K.left_turn(p1,p2,s.target()) )
        { v = CT.source(e); current = VERTEX; }
        else if ( or2 == 0 && !this->K.left_turn(p2,p3,s.target()) )
        { v = CT.target(e); current = VERTEX; }
        else if ( or3 == 0 && !this->K.left_turn(p3,p1,s.target()) )
        { v = CT.target(CT.next(e)); current = VERTEX; }
        else if ( or2 > 0 && or1 < 0 && !this->K.left_turn(p1,p2,s.target()) )
        { e = CT.twin(e); current = EDGE_CROSSING; }
        else if ( or3 > 0 && or2 < 0 && !this->K.left_turn(p2,p3,s.target()) )
        { e = CT.twin(CT.next(e)); current = EDGE_CROSSING; }
        else if ( or1 > 0 && or3 < 0 && !this->K.left_turn(p3,p1,s.target()) )
        { e = CT.twin(CT.previous(e)); current = EDGE_CROSSING; }
        else return Object_handle();

      }
    } else {

      if(check_tag(typename Is_extended_kernel<Geometry>::value_type())) {
	CGAL_error_msg( "code is only for Bounded_kernel");
      }
      if(!ray_shoot_from_outer_facet(s,current,v,e,typename Is_extended_kernel<Geometry>::value_type()))
	return Object_handle();
    }

    CGAL_NEF_TRACEN("current = " << current);
    if(current == VERTEX){
      CGAL_NEF_TRACEN(point(v));
    }
    while (true) switch ( current ) {
      case VERTEX:
        { CGAL_NEF_TRACEN("vertex "<<CT.point(v));
          Vertex_const_handle v_org = input_vertex(v);
          if ( M(v_org) ) return make_object(v_org);
          if ( CT.point(v) == s.target() ) return Object_handle();
          // stop walking at s.target(), or determine next object on s:
          bool collinear;
          Halfedge_const_handle e_out = CT.out_wedge(v,d,collinear);
          if (collinear) // ray shoot via e_out
          { e = e_out; current = EDGE_COLLINEAR; }
          else { // ray shoot in wedge left of e_out
            if ( M(input_face(e_out)) )
              return make_object(input_face(e_out));
            e = CT.twin(CT.next(e_out)); current = EDGE_CROSSING;
          }
        }

        break;
      case EDGE_CROSSING:
        { CGAL_NEF_TRACEN("crossing edge "<<segment(e));
          if ( this->K.orientation(CT.segment(e),s.target()) == 0 )
            return Object_handle();
          Halfedge_const_handle e_org = input_halfedge(e);
          if ( e_org != Halfedge_const_handle() ) { // not a CT edge
            if ( M(e_org) ) return make_object(e_org);
            if ( M(face(e_org)) ) return make_object(face(e_org));
          }
          Vertex_const_handle v_cand = CT.target(CT.next(e));
          CGAL_NEF_TRACEN("v_cand "<<PV(v_cand));
          int orientation_ = this->K.orientation(p,s.target(),CT.point(v_cand));
          switch( orientation_ ) {
            case 0:
              v = v_cand; current = VERTEX; break;
            case +1:
              e = CT.twin(CT.next(e)); current = EDGE_CROSSING; break;
            case -1:
              e = CT.twin(CT.previous(e)); current = EDGE_CROSSING; break;
          }
        }

        break;
      case EDGE_COLLINEAR:
        { CGAL_NEF_TRACEN("collinear edge "<<CT.segment(e));
          Halfedge_const_handle e_org = input_halfedge(e);
          if ( e_org == Halfedge_const_handle() ) { // a CT edge
            if ( M(input_face(e)) )
              return make_object(input_face(e));
          } else { // e_org is not a CT edge
            if ( M(e_org) )
              return make_object(e_org);
          }
          if ( this->K.strictly_ordered_along_line(
                 CT.point(CT.source(e)),s.target(),CT.point(CT.target(e))) )
            return Object_handle();
          v = CT.target(e); current = VERTEX;
        }

        break;
    }
    // CGAL_error(); return h; // compiler warning
  }

  bool within_outer_cycle(Vertex_const_handle ,
			  const Point& , const Tag_true& ) const {
    return true;
  }

  bool within_outer_cycle(Vertex_const_handle v,
			  const Point& q, const Tag_false& ) const {
    typedef Project_halfedge_point<typename Decorator::Halfedge, Point> Project;
    typedef Circulator_project<Halfedge_around_face_const_circulator,
      Project, const Point&, const Point*> Circulator;
    typedef Container_from_circulator<Circulator> Container;

    Halfedge_const_handle e_min = CT.twin(CT.first_out_edge(v));
    Halfedge_around_face_const_circulator circ(e_min);
    Circulator c(circ);
    Container ct(c);
    if(is_empty_range(ct.begin(), ct.end()) ||
       bounded_side_2(ct.begin(), ct.end(),q) == CGAL::ON_UNBOUNDED_SIDE)
      return false;

    return true;
  }

  Object_handle walk_in_triangulation(const Point& p) const;

}; // PM_point_locator<PM_decorator_,Geometry_>


#ifdef CGAL_USING_PPL
static const char* const pointlocationversion ="point location via pers dicts";
#else
static const char* const pointlocationversion ="point location via seg walks";
#endif

template <typename PMD, typename GEO>
PM_point_locator<PMD,GEO>::
PM_point_locator(const Plane_map& P, const Geometry& k) :
  Base(P,k), CT(*(new Plane_map),k)

{ CGAL_NEF_TRACEN("PM_point_locator construction");
  CT.clone_skeleton(P,CT_link_to_original(CT,*this));
  triangulate_CT();
  minimize_weight_CT();
  #ifdef CGAL_USING_PPL
  pPPL = new PMPP_locator(CT,PMPPLT(K));
  #endif

}

template <typename PMD, typename GEO>
PM_point_locator<PMD,GEO>::
~PM_point_locator()
{ CGAL_NEF_TRACEN("clear_static_point_locator");
  Vertex_iterator vit, vend = CT.vertices_end();
  for (vit = CT.vertices_begin(); vit != vend; ++vit) {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<VF_pair>::clear(CT.info(vit));
    #else
    CT.info(vit)=boost::any();
    #endif
  }
  Halfedge_iterator eit, eend = CT.halfedges_end();
  for (eit = CT.halfedges_begin(); eit != eend; ++eit) {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    geninfo<EF_pair>::clear(CT.info(eit));
    #else
    CT.info(eit)=boost::any();
    #endif
  }
  CT.clear();
  delete &(CT.plane_map());
  #ifdef CGAL_USING_PPL
  delete pPPL; pPPL=0;
  #endif
}

template <typename PMD, typename GEO>
typename PM_point_locator<PMD,GEO>::Object_handle
PM_point_locator<PMD,GEO>::walk_in_triangulation(const Point& q) const
{
  CGAL_NEF_TRACEN("walk in triangulation "<<q);

  Vertex_const_handle v = CT.vertices_begin();

  if(!check_tag(typename Is_extended_kernel<GEO>::value_type()))
    if(!within_outer_cycle(v,q,typename Is_extended_kernel<Geometry>::value_type()))
      return Object_handle();

  Halfedge_const_handle e;
  Point p = CT.point(v);
  if ( p == q ) return make_object(v);
  //  Segment s = this->K.construct_segment(p,q);
  Direction dir = this->K.construct_direction(p,q);
  object_kind current = VERTEX;
  while (true) switch ( current ) {
    case VERTEX:
      {
        CGAL_NEF_TRACEN("vertex "<<CT.point(v));
        if ( CT.point(v) == q )
          return make_object(v); // stop walking at q
        bool collinear;
        Halfedge_const_handle e_out = CT.out_wedge(v,dir,collinear);
        if (collinear) // ray shoot via e_out
        { e = e_out; current = EDGE_COLLINEAR; }
        else  // ray shoot in wedge left of e_out
        { e = CT.twin(CT.next(e_out)); current = EDGE_CROSSING; }
      }

      break;
    case EDGE_CROSSING:
      { CGAL_NEF_TRACEN("crossing edge "<<CT.segment(e));
        if ( !(this->K.orientation(CT.segment(e),q) > 0) ) // q not left of e
          return make_object(e);
        Vertex_const_handle v_cand = CT.target(CT.next(e));
        int orientation_ = this->K.orientation(p,q,CT.point(v_cand));
        switch( orientation_ ) {
          case 0:  // collinear
            if ( this->K.strictly_ordered_along_line(p,q,CT.point(v_cand)) )
              return make_object(e);
            v = v_cand; current = VERTEX; break;
          case +1: // left_turn
            e = twin(next(e)); current = EDGE_CROSSING; break;
          case -1:
            e = twin(previous(e)); current = EDGE_CROSSING; break;
        }
      }

      break;
    case EDGE_COLLINEAR:
      { CGAL_NEF_TRACEN("collinear edge "<<CT.segment(e));
        if ( this->K.strictly_ordered_along_line(
               CT.point(CT.source(e)),q,CT.point(CT.target(e))) )
          return make_object(e);
        v = CT.target(e); current = VERTEX;
      }

      break;
  }
  return Object_handle(); // never reached warning acceptable
}

} //namespace CGAL

#endif // CGAL_PM_POINT_LOCATOR_H
