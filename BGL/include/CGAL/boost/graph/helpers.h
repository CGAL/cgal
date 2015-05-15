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

namespace CGAL {



template <typename FaceGraph>
bool is_border(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return face(hd,g) == boost::graph_traits<FaceGraph>::null_face();
}

template <typename FaceGraph>
bool is_border_edge(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return is_border(hd, g) || is_border(opposite(hd,g), g);
}

template <typename FaceGraph>
bool is_border(typename boost::graph_traits<FaceGraph>::edge_descriptor ed, const FaceGraph& g)
{
  return is_border_edge(halfedge(ed,g), g);
}

template <typename Graph>
boost::optional<typename boost::graph_traits<Graph>::halfedge_descriptor>
is_border(typename boost::graph_traits<Graph>::vertex_descriptor v,
          const Graph& g)
{
  CGAL::Halfedge_around_target_iterator<Graph> havib, havie;
  for(boost::tie(havib, havie) = halfedges_around_target(halfedge(v, g), g); havib != havie; ++havib) {
    if(is_border(*havib,g)) {
      typename boost::graph_traits<Graph>::halfedge_descriptor h = *havib;
      return h;
    }
  }
  // empty
  return boost::optional<typename boost::graph_traits<Graph>::halfedge_descriptor>();
}


 /*!
    returns `true` if there are no 
    border edges. 
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


template <typename FaceGraph>
bool is_bivalent(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == opposite(next(opposite(next(hd,g),g),g),g);
}

  /*!
    returns `true` if all 
    vertices have exactly two incident edges. 
  */ 
template <typename FaceGraph>
  bool is_pure_bivalent(const FaceGraph& g)  
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

template <typename FaceGraph>
bool is_trivalent(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
  return hd == opposite(next(opposite(next(opposite(next(hd,g),g),g),g),g),g);
}
	
  /*!
    returns `true` if all 
    vertices have exactly three incident edges. 
  */ 
template <typename FaceGraph>
  bool is_pure_trivalent(const FaceGraph& g)  
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(vertex_descriptor vd, vertices(g)){
    halfedge_descriptor hd = halfedge(vd,g);
    if((hd == boost::graph_traits<FaceGraph>::null_halfedge()) ||
       (! is_trivalent(hd,g))){
      return false;
    }
  }
  return true;
}

 /*!
    returns `true` iff the connected component denoted by `h` is a triangle. 
    \pre `g` must be valid.
  */ 
template <typename FaceGraph>
  bool is_triangle(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)  
{ 
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor beg = hd;
  if(is_border(hd,g)) return false;
  for(int i=0; i<3;i++){
    if(! is_border(opposite(hd,g),g)) return false;
    hd = next(hd,g);
  }
  return next(hd,g)== beg;
}

 /*!
    returns `true` iff the face is a triangle, that is it has three incident halfedges. 
 */
template <typename FaceGraph>
bool is_triangle(typename boost::graph_traits<FaceGraph>::face_descriptor fd, const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor hd = halfedge(fd,g);
  return hd == next(next(next(hd,g),g),g);
}

  /*!
    returns `true` if all faces are triangles. 
  */ 
template <typename FaceGraph>
  bool is_pure_triangle(const FaceGraph& g)  
{
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(g)){
    if(! is_triangle(fd,g)){
      return false;
    }
  }
  return true;
}


/*!
    returns `true` iff the connected component denoted by `h` is a quadrilateral. 
  */
template <typename FaceGraph>
bool is_quad(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph& g)
{
 typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor beg = hd;
  if(is_border(hd,g)) return false;
  for(int i=0; i<4;i++){
    if(! is_border(opposite(hd,g),g)) return false;
    hd = next(hd,g);
  }
  return next(hd,g)== beg;
}


 /*!
    returns `true` iff the face is a quad, that is it has four incident halfedges. 
 */
template <typename FaceGraph>
bool is_quad(typename boost::graph_traits<FaceGraph>::face_descriptor fd, const FaceGraph& g)
{
 typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor hd = halfedge(fd,g);
  return hd == next(next(next(next(hd,g),g),g),g);
}

  /*!
    returns `true` if all faces are quadrilaterals. 
  */ 
template <typename FaceGraph>
  bool is_pure_quad(const FaceGraph& g)  
{
    typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(g)){
    if(! is_quad(fd,g)){
      return false;
    }
  }
  return true;
}
 
  /*!
    returns `true` iff the connected component denoted by `h` is a tetrahedron. 
  */ 
template <typename FaceGraph>
bool is_tetrahedron( typename boost::graph_traits<FaceGraph>::halfedge_descriptor h1, const FaceGraph& g)   
{
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
bool is_valid( typename boost::graph_traits<FaceGraph>::halfedge_descriptor h, const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
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
bool is_valid( typename boost::graph_traits<FaceGraph>::vertex_descriptor v, const FaceGraph& g)
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
bool is_valid( typename boost::graph_traits<FaceGraph>::face_descriptor f, const FaceGraph& g)
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
    if(! is_valid(v,g)){
      return false;
    }
  }
  BOOST_FOREACH(halfedge_descriptor h, halfedges(g)){
    if(! is_valid(h,g)){
      return false;
    }
  }
  BOOST_FOREACH(face_descriptor f, faces(g)){
    if(! is_valid(f,g)){
      return false;
    }
  }
  return true;
}


} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_HELPERS_H
