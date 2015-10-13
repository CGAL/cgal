// Copyright (c) 2014 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s) : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_HELPERS_H
#define CGAL_BOOST_GRAPH_HELPERS_H


#include <boost/foreach.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>


namespace CGAL {

  namespace Euler {

    template< typename Graph>
    void fill_hole(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                   Graph& g);
  }

/*!
   \ingroup PkgBGLHelperFct
    returns `true` if the halfedge `hd` is on a border. 
  */
template <typename FaceGraph>
bool is_border(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return face(hd,g) == boost::graph_traits<FaceGraph>::null_face();
}

 /*!
   \ingroup PkgBGLHelperFct
    returns `true` if the halfedge `hd` or the opposite halfedge is on a border. 
  */
template <typename FaceGraph>
bool is_border_edge(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return is_border(hd, g) || is_border(opposite(hd,g), g);
}

 /*!
   \ingroup PkgBGLHelperFct
    returns `true` if the edge `e` is on a border. 
  */
template <typename FaceGraph>
bool is_border(typename boost::graph_traits<FaceGraph>::edge_descriptor ed, const FaceGraph& g)
{
  return is_border_edge(halfedge(ed,g), g);
}

 /*!
   \ingroup PkgBGLHelperFct
    returns a halfedge which is on a border and whose target vertex is `vd`, if such a halfedge exists. 
  */
template <typename FaceGraph>
boost::optional<typename boost::graph_traits<FaceGraph>::halfedge_descriptor>
is_border(typename boost::graph_traits<FaceGraph>::vertex_descriptor vd,
          const FaceGraph& g)
{
  CGAL::Halfedge_around_target_iterator<FaceGraph> havib, havie;
  for(boost::tie(havib, havie) = halfedges_around_target(halfedge(vd, g), g); havib != havie; ++havib) {
    if(is_border(*havib,g)) {
      typename boost::graph_traits<FaceGraph>::halfedge_descriptor h = *havib;
      return h;
    }
  }
  // empty
  return boost::optional<typename boost::graph_traits<FaceGraph>::halfedge_descriptor>();
}


 /*!
   \ingroup PkgBGLHelperFct
    returns `true` if there are no border edges. 
  */
template <typename FaceGraph>
bool is_closed(const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(g)){
    if(is_border(hd,g)){
      return false;
    }
  }
  return true;
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if the target of `hd` has exactly two incident edges. 
  */ 
template <typename FaceGraph>
bool is_bivalent(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == opposite(next(opposite(next(hd,g),g),g),g);
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if all vertices have exactly two incident edges. 
  */ 
template <typename FaceGraph>
  bool is_bivalent_mesh(const FaceGraph& g)  
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    halfedge_descriptor hd = halfedge(vd,g);
    if((hd == boost::graph_traits<FaceGraph>::null_halfedge()) ||
       (! is_bivalent(hd,g))){
      return false;
    }
  }
  return true;
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if the target of `hd` has exactly three incident edges. 
  */ 
template <typename FaceGraph>
bool is_trivalent(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == opposite(next(opposite(next(opposite(next(hd,g),g),g),g),g),g);
}
	
  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if all 
    vertices have exactly three incident edges. 
  */ 
template <typename FaceGraph>
  bool is_trivalent_mesh(const FaceGraph& g)  
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    halfedge_descriptor hd = halfedge(vd,g);
    if((hd == boost::graph_traits<FaceGraph>::null_halfedge()) ||
       (! is_trivalent(halfedge(hd,g),g))){
      return false;
    }
  }
  return true;
}

 /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the connected component denoted by `hd` is a triangle. 
    \pre `g` must be valid.
  */ 
template <typename FaceGraph>
  bool is_isolated_triangle(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)  
{ 
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor beg = hd;
  if(is_border(hd,g)) return false;
  for(int i=0; i<3;i++){
    if(! is_border(opposite(hd,g),g)) return false;
    hd = next(hd,g);
  }
  return hd == beg;
}

 /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the face denoted by `hd` is a triangle, that is it has three incident halfedges. 
 */
template <typename FaceGraph>
bool is_triangle(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == next(next(next(hd,g),g),g);
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if all faces are triangles. 
  */ 
template <typename FaceGraph>
  bool is_triangle_mesh(const FaceGraph& g)  
{
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(g)){
    if(! is_triangle(halfedge(fd,g),g)){
      return false;
    }
  }
  return true;
}

/*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the connected component denoted by `hd` is a quadrilateral. 
  */
template <typename FaceGraph>
bool is_isolated_quad(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
 typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor beg = hd;
  if(is_border(hd,g)) return false;
  for(int i=0; i<4;i++){
    if(! is_border(opposite(hd,g),g)) return false;
    hd = next(hd,g);
  }
  return hd == beg;
}


 /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the face denoted by `hd` is a quad, that is it has four incident halfedges. 
 */
template <typename FaceGraph>
bool is_quad(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == next(next(next(next(hd,g),g),g),g);
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` if all faces are quadrilaterals. 
  */ 
template <typename FaceGraph>
  bool is_quad_mesh(const FaceGraph& g)  
{
    typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(g)){
    if(! is_quad(halfedge(fd,g),g)){
      return false;
    }
  }
  return true;
}
 
  /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the connected component denoted by `hd` is a tetrahedron. 
  */ 
template <typename FaceGraph>
bool is_tetrahedron( typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)   
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h1 = hd;
  if(is_border(h1,g)) return false;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h2 = next(h1,g);
  halfedge_descriptor h3 = next(h2,g);
  halfedge_descriptor h4 = next(opposite(h1,g),g );
  halfedge_descriptor h5 = next(opposite(h2,g),g );
  halfedge_descriptor h6 = next(opposite(h3,g),g );
  // check halfedge combinatorics.
  // at least three edges at vertices 1, 2, 3.
  if ( h4 == opposite(h3,g) ) return false;
  if ( h5 == opposite(h1,g) ) return false;
  if ( h6 == opposite(h2,g) ) return false;
  // exact three edges at vertices 1, 2, 3.
  if ( next(opposite(h4,g),g) != opposite(h3,g) ) return false;
  if ( next(opposite(h5,g),g) != opposite(h1,g) ) return false;
  if ( next(opposite(h6,g),g) != opposite(h2,g) ) return false;
  // three edges at v4.
  if ( opposite(next(h4,g),g) != h5 ) return false;
  if ( opposite(next(h5,g),g) != h6 ) return false;
  if ( opposite(next(h6,g),g) != h4 ) return false;
  // All facets are triangles.
  if ( next(next(next(h1,g),g),g) != h1 ) return false;
  if ( next(next(next(h4,g),g),g) != h4 ) return false;
  if ( next(next(next(h5,g),g),g) != h5 ) return false;
  if ( next(next(next(h6,g),g),g) != h6 ) return false;
  // all edges are non-border edges.
  if ( is_border(h1,g) ) return false;  // implies h2 and h3
  if ( is_border(h4,g) ) return false;
  if ( is_border(h5,g) ) return false;
  if ( is_border(h6,g) ) return false;
  return true;
  }

template <typename FaceGraph>
bool is_valid_halfedge_descriptor( typename boost::graph_traits<FaceGraph>::halfedge_descriptor h, const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  face_descriptor f = face(h,g);
  halfedge_descriptor done(h);
  do{
    if(face(h,g) != f){
      std::cerr << "halfedge " << h << " is invalid\n";
      return false;
    }
    halfedge_descriptor hn = h;
    hn = next(h,g);
    if(prev(hn,g) != h){
      std::cerr << "halfedge " << h << " is invalid\n";
      return false;
    }
    h = hn;
  } while(h != done);
  return true;
}

template <typename FaceGraph>
bool is_valid_vertex_descriptor( typename boost::graph_traits<FaceGraph>::vertex_descriptor v, const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h = halfedge(v,g), done(h);
  if(h == boost::graph_traits<FaceGraph>::null_halfedge()){
    return true;
  }
  do{
    if(target(h,g) != v){
      std::cerr << "vertex " << v << " is invalid\n";
      return false;
    }
    h = opposite(next(h,g),g);
  }while(h != done);
  return true;
}

template <typename FaceGraph>
bool is_valid_face_descriptor( typename boost::graph_traits<FaceGraph>::face_descriptor f, const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(f,g);
  if(face(h,g) != f){
    std::cerr << "face " << f << " is invalid\n";
    return false;
  }
  return true;
}


template <typename FaceGraph>
bool is_valid_polygon_mesh(const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor     face_descriptor;
  BOOST_FOREACH(vertex_descriptor v, vertices(g)){
    if(! is_valid_vertex_descriptor(v,g)){
      return false;
    }
  }
  BOOST_FOREACH(halfedge_descriptor h, halfedges(g)){
    if(! is_valid_halfedge_descriptor(h,g)){
      return false;
    }
  }
  BOOST_FOREACH(face_descriptor f, faces(g)){
    if(! is_valid_face_descriptor(f,g)){
      return false;
    }
  }
  return true;
}

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the connected component denoted by `hd` is a hexahedron. 
  */ 
template <typename FaceGraph>
bool is_hexahedron( typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)   
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h1 = hd;
  if(is_border(h1,g)) return false;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor h2 = next(h1,g);
  halfedge_descriptor h3 = next(h2,g);
  halfedge_descriptor h4 = next(h3,g);
  halfedge_descriptor h1o = opposite(h1,g);
  halfedge_descriptor h2o = opposite(h2,g);
  halfedge_descriptor h3o = opposite(h3,g);
  halfedge_descriptor h4o = opposite(h4,g);
  if(opposite(next(h2o,g),g) != prev(h1o,g)) return false;
  if(opposite(next(h3o,g),g) != prev(h2o,g)) return false;
  if(opposite(next(h4o,g),g) != prev(h3o,g)) return false;
  if(opposite(next(h1o,g),g) != prev(h4o,g)) return false;
  if(! is_quad(h1,g)) return false;
  if(! is_quad(h1o,g)) return false;
  if(! is_quad(h2o,g)) return false;
  if(! is_quad(h3o,g)) return false;
  if(! is_quad(h4o,g)) return false;
  h1o =next(next(h1o,g),g);
  h2o =next(next(h2o,g),g);
  h3o =next(next(h3o,g),g);
  h4o =next(next(h4o,g),g);
  if(next(opposite(h2o,g),g) != opposite(h1o,g)) return false;
  if(next(opposite(h3o,g),g) != opposite(h2o,g)) return false;
  if(next(opposite(h4o,g),g) != opposite(h3o,g)) return false;
  if(next(opposite(h1o,g),g) != opposite(h4o,g)) return false;

  if(! is_quad(opposite(h4o,g),g)) return false;
  return true;
}



/** 
 * \ingroup PkgBGLHelperFct
 * Creates an isolated triangle with border edges in `g` having `p0`, `p1`, and `p2` as points and adds it to the graph `g`.
 * \returns the non-border halfedge which has the target vertex associated with `p0`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_triangle(const P& p0, const P& p1, const P& p2, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
  typedef typename Traits::face_descriptor                 face_descriptor;
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  Point_property_map ppmap = get(CGAL::vertex_point, g);
  vertex_descriptor v0, v1, v2;
  v0 = add_vertex(g);
  v1 = add_vertex(g);
  v2 = add_vertex(g);

  ppmap[v0] = p0;
  ppmap[v1] = p1;
  ppmap[v2] = p2;
  halfedge_descriptor h0 = halfedge(add_edge(g),g);
  halfedge_descriptor h1 = halfedge(add_edge(g),g);
  halfedge_descriptor h2 = halfedge(add_edge(g),g);
  set_next(h0, h1, g);
  set_next(h1, h2, g);
  set_next(h2, h0, g);
  set_target(h0, v1, g);
  set_target(h1, v2, g);
  set_target(h2, v0, g);
  set_halfedge(v1, h0, g);
  set_halfedge(v2, h1, g);
  set_halfedge(v0, h2, g);
  face_descriptor f = add_face(g);
  set_face(h0,f,g);
  set_face(h1,f,g);
  set_face(h2,f,g);
  set_halfedge(f,h0,g);
  h0 = opposite(h0,g);
  h1 = opposite(h1,g);
  h2 = opposite(h2,g);
  set_next(h0, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  set_face(h0, boost::graph_traits<Graph>::null_face(),g);
  set_face(h1, boost::graph_traits<Graph>::null_face(),g);
  set_face(h2, boost::graph_traits<Graph>::null_face(),g);
  return opposite(h2,g);
}

/** 
 * \ingroup PkgBGLHelperFct
 * Creates an isolated quad with border edges in `g` having `p0`, `p1`, `p2`, and `p3` as points and adds it to the graph `g`.
 * \returns the non-border halfedge which has the target vertex associated with `p0`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_quad(const P& p0, const P& p1, const P& p2, const P& p3, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
  typedef typename Traits::face_descriptor                 face_descriptor;
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  Point_property_map ppmap = get(CGAL::vertex_point, g);
  vertex_descriptor v0, v1, v2, v3;
  v0 = add_vertex(g);
  v1 = add_vertex(g);
  v2 = add_vertex(g);
  v3 = add_vertex(g);

  ppmap[v0] = p0;
  ppmap[v1] = p1;
  ppmap[v2] = p2;
  ppmap[v3] = p3;
  halfedge_descriptor h0 = halfedge(add_edge(g),g);
  halfedge_descriptor h1 = halfedge(add_edge(g),g);
  halfedge_descriptor h2 = halfedge(add_edge(g),g);
  halfedge_descriptor h3 = halfedge(add_edge(g),g);
  set_next(h0, h1, g);
  set_next(h1, h2, g);
  set_next(h2, h3, g);
  set_next(h3, h0, g);
  set_target(h0, v1, g);
  set_target(h1, v2, g);
  set_target(h2, v3, g);
  set_target(h3, v0, g);
  set_halfedge(v1, h0, g);
  set_halfedge(v2, h1, g);
  set_halfedge(v3, h2, g);
  set_halfedge(v0, h3, g);
  face_descriptor f = add_face(g);
  set_face(h0,f,g);
  set_face(h1,f,g);
  set_face(h2,f,g);
  set_face(h3,f,g);
  set_halfedge(f,h0,g);
  h0 = opposite(h0,g);
  h1 = opposite(h1,g);
  h2 = opposite(h2,g);
  h3 = opposite(h3,g);
  set_next(h0, h3, g);
  set_next(h3, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  set_target(h3, v3, g);
  set_face(h0, boost::graph_traits<Graph>::null_face(),g);
  set_face(h1, boost::graph_traits<Graph>::null_face(),g);
  set_face(h2, boost::graph_traits<Graph>::null_face(),g);
  set_face(h3, boost::graph_traits<Graph>::null_face(),g);
  return opposite(h3,g);
}

/** 
 * \ingroup PkgBGLHelperFct
 * Creates an isolated hexahedron in `g` having `p0`, `p1`, ...\ , and `p7` as points and adds it to the graph `g`.
 * \returns the halfedge which has the target vertex associated with `p0`, in the face with the vertices with the points `p0`, `p1`, `p2`, and `p3`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_hexahedron(const P& p0, const P& p1, const P& p2, const P& p3,
                const P& p4, const P& p5, const P& p6, const P& p7, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;

  halfedge_descriptor hb = make_quad(p0, p1, p2, p3, g);
  halfedge_descriptor ht = prev(make_quad(p4, p7, p6, p5, g),g);
  for(int i=0; i <4; i++){
    halfedge_descriptor h = halfedge(add_edge(g),g);
    set_target(h,target(hb,g),g);
    set_next(h,opposite(hb,g),g);
    set_next(opposite(next(ht,g),g),h,g);
    h = opposite(h,g);
    set_target(h,target(ht,g),g);
    set_next(h,opposite(ht,g),g);
    set_next(opposite(next(hb,g),g),h,g);
    hb = next(hb,g);
    ht = prev(ht,g);
  }
  for(int i=0; i <4; i++){
    Euler::fill_hole(opposite(hb,g),g);
    hb = next(hb,g);
  }
  return hb;
}
/** 
 * \ingroup PkgBGLHelperFct
 * Creates an isolated tetrahedron in `g` having `p0`, `p1`, `p2`, and `p3` as points and adds it to the graph `g`.
 * \returns the halfedge which has the target vertex associated with `p0`, in the face with the vertices with the points `p0`, `p1`, and `p2`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_tetrahedron(const P& p0, const P& p1, const P& p2, const P& p3, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
  typedef typename Traits::face_descriptor                 face_descriptor;
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;

  Point_property_map ppmap = get(CGAL::vertex_point, g);
  vertex_descriptor v0, v1, v2, v3;
  v0 = add_vertex(g);
  v1 = add_vertex(g);
  v2 = add_vertex(g);
  v3 = add_vertex(g);

  ppmap[v0] = p0;
  ppmap[v1] = p1;
  ppmap[v2] = p2;
  ppmap[v3] = p3;
  halfedge_descriptor h0 = halfedge(add_edge(g),g);
  halfedge_descriptor h1 = halfedge(add_edge(g),g);
  halfedge_descriptor h2 = halfedge(add_edge(g),g);
  set_next(h0, h1, g);
  set_next(h1, h2, g);
  set_next(h2, h0, g);
  set_target(h0, v1, g);
  set_target(h1, v2, g);
  set_target(h2, v0, g);
  set_halfedge(v1, h0, g);
  set_halfedge(v2, h1, g);
  set_halfedge(v0, h2, g);
  face_descriptor f = add_face(g);
  set_face(h0,f,g);
  set_face(h1,f,g);
  set_face(h2,f,g);
  set_halfedge(f,h0,g);
  h0 = opposite(h0,g);
  h1 = opposite(h1,g);
  h2 = opposite(h2,g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  halfedge_descriptor h3 = halfedge(add_edge(g),g);
  halfedge_descriptor h4 = halfedge(add_edge(g),g);
  halfedge_descriptor h5 = halfedge(add_edge(g),g);
  set_target(h3, v3, g);
  set_target(h4, v3, g);
  set_target(h5, v3, g);
  set_halfedge(v3, h3, g);
  
  set_next(h0, h3, g);
  set_next(h1, h4, g);
  set_next(h2, h5, g);

  set_next(h3, opposite(h4,g), g);
  set_next(h4, opposite(h5,g), g);
  set_next(h5, opposite(h3,g), g);
  set_next(opposite(h4,g), h0, g);
  set_next(opposite(h5,g), h1, g);
  set_next(opposite(h3,g), h2, g);

  set_target(opposite(h3,g), v0, g);
  set_target(opposite(h4,g), v1, g);
  set_target(opposite(h5,g), v2, g);

  f = add_face(g);
  set_halfedge(f,h0,g);
  set_face(h0, f, g);
  set_face(h3, f, g);
  set_face(opposite(h4,g), f, g);
  f = add_face(g);
  set_halfedge(f,h1,g);
  set_face(h1, f, g);
  set_face(h4, f, g);
  set_face(opposite(h5,g), f, g);
  f = add_face(g);
  set_halfedge(f,h2,g);
  set_face(h2, f, g);
  set_face(h5, f, g);
  set_face(opposite(h3,g), f, g);
  
  return opposite(h2,g);
}


} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_HELPERS_H
