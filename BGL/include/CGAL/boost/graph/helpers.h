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
// SPDX-License-Identifier: LGPL-3.0+
// 
// Author(s) : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_HELPERS_H
#define CGAL_BOOST_GRAPH_HELPERS_H


#include <boost/foreach.hpp>
#include <boost/range/empty.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/internal/Has_member_clear.h>
#include <CGAL/function_objects.h>
#include <boost/unordered_set.hpp>
#include <CGAL/IO/Verbose_ostream.h>


namespace CGAL {

  namespace Euler {

    template< typename Graph>
    void fill_hole(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                   Graph& g);

    template<typename Graph , typename VertexRange >
    typename boost::graph_traits<Graph>::face_descriptor add_face(const VertexRange& vr,
                                                         Graph& g);
  }//Euler

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

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of `g`.
 *
 * `g` is valid if it follows the rules of the `HalfedgeListGraph` concept,
 * and all of its associations are reciprocal.
 * For example, `prev(next(h, g), g)` must be `h`,
 * and `next(prev(h, g), g)` must be `h`.
 * \param g the `Graph` to test.
 * \param verb : if `true`, the details of the check will be written in the standard output.
 *
 * \tparam `Graph` a model of `HalfedgeListGraph`
 * \return `true` if `g` is valid, `false` otherwise.
 *
 */
template<typename Graph>
bool is_valid_halfedge_graph(const Graph& g, bool verb = false)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertices_size_type    vertex_size_type;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedges_size_type   halfedges_size_type;

  Verbose_ostream verr(verb);
  std::size_t num_v(std::distance(boost::begin(vertices(g)), boost::end(vertices(g)))),
              num_e(std::distance(boost::begin(edges(g)), boost::end(edges(g)))),
              num_h(std::distance(boost::begin(halfedges(g)), boost::end(halfedges(g))));

  bool valid = (1 != (num_h&1) && (2*num_e == num_h));
  if(!valid)
  {
    verr << "number of halfedges is odd." << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  // All halfedges.
  halfedges_size_type n = 0;
  BOOST_FOREACH(halfedge_descriptor begin, halfedges(g))
  {
    // Pointer integrity.
    valid = (next(begin, g) != boost::graph_traits<Graph>::null_halfedge());
    valid = valid && (opposite(begin, g) != boost::graph_traits<Graph>::null_halfedge());
    if(!valid)
    {
      verr << "halfedge " << n << " next / opposite halfedges are null." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // edge integrity
    valid = (halfedge(edge(begin, g), g) == begin);

    // opposite integrity.
    valid = valid && (opposite(begin, g) != begin);
    valid = valid && (opposite(opposite(begin, g), g) == begin);
    if(!valid)
    {
      verr << "halfedge " << n << " invalid halfedge opposite()." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // previous integrity.
    valid = (prev(next(begin, g), g) == begin);
    valid = valid && (next(prev(begin, g), g) == begin);
    if(!valid)
    {
      verr << "halfedge " << n << " prev(next(hd)) != hd OR next(prev(hd)) != hd" << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // vertex integrity.
    valid = (target(begin, g) != boost::graph_traits<Graph>::null_vertex());
    if(!valid)
    {
      verr << "halfedge " << n << " target of halfedge is the null vertex." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    valid = (target(begin, g) == target(opposite(next(begin, g), g), g));
    if(!valid)
    {
      verr << "halfedge " << n << " target(hd) != source(next(hd))." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    ++n;
  }

  valid = (n == num_h);
  if(!valid)
  {
    verr << "counting halfedges failed." << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  // All vertices.
  vertex_size_type v = 0;
  n = 0;
  BOOST_FOREACH(vertex_descriptor vbegin, vertices(g))
  {
    // Pointer integrity.
    if(halfedge(vbegin, g) != boost::graph_traits<Graph>::null_halfedge())
      valid = (target(halfedge(vbegin, g), g) == vbegin);
    else
      valid = false;

    if(!valid)
    {
      verr << "vertex " << v << " halfedge incident to vertex is the null halfedge." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // cycle-around-vertex test.
    halfedge_descriptor h = halfedge(vbegin, g);
    if(h != boost::graph_traits<Graph>::null_halfedge())
    {
      halfedge_descriptor ge = h;
      do
      {
        ++n;
        h = opposite(next(h, g), g);
        valid = (n <= num_h && n != 0);
        if(!valid)
        {
          verr << "vertex " << v << " too many halfedges around vertex." << std::endl;
          verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
          return false;
        }
      }
      while(h != ge);
    }

    ++v;
  }

  valid = (v == num_v);
  if(!valid)
  {
    verr << "counting vertices failed." << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  valid = (n == num_h);
  if(!valid)
  {
    verr << "counting halfedges via vertices failed." << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  // All halfedges.
  n = 0;
  BOOST_FOREACH(halfedge_descriptor i, halfedges(g))
  {
    // At least triangular facets and distinct geometry.
    valid = (next(i, g) != i) && (target(i, g) != target(opposite(i, g), g));
    if(!valid)
    {
      verr << "halfedge " << n << " pointer validity corrupted." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    ++n;
  }

  valid = (n == num_h);
  if(!valid)
    verr << "counting halfedges failed." << std::endl;

  verr << "Halfedge Graph Structure is " << (valid ? "valid." : "NOT VALID.") << std::endl;

  return valid;
}

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of `g`.
 *
 * `g` is valid if it is a valid `HalfedgeListGraph`, if it follows the rules
 * of the `FaceListGraph` concept, and all of its associations are reciprocal.
 * For example, `face(halfedge(f,g),g)` must be `f`.
 * calls `is_valid_halfedge_graph()`
 * \param g the `Graph` to test.
 * \param verb : if `true`, the details of the check will be written in the standard output.
 *
 * \tparam `Graph` a model of `FaceListGraph`
 * \return `true` if `g` is valid, `false` otherwise.
 *
 * \see `is_valid_halfedge_graph()`
 */
template<typename Graph>
bool is_valid_face_graph(const Graph& g, bool verb = false)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedges_size_type   halfedges_size_type;
  typedef typename boost::graph_traits<Graph>::face_descriptor       face_descriptor;
  typedef typename boost::graph_traits<Graph>::faces_size_type       faces_size_type;

  Verbose_ostream verr(verb);

  std::size_t num_f(std::distance(boost::begin(faces(g)), boost::end(faces(g)))),
              num_h(std::distance(boost::begin(halfedges(g)), boost::end(halfedges(g))));

  faces_size_type f = 0;
  std::size_t n = 0;
  std::size_t hn = 0;
  halfedges_size_type nb = 0;

  //is valid halfedge_graph ?
  bool valid = is_valid_halfedge_graph(g, verb);
  if(!valid)
    return false;

  // All faces.
  BOOST_FOREACH(face_descriptor fbegin, faces(g))
  {
    // Pointer integrity.
    if(halfedge(fbegin, g) != boost::graph_traits<Graph>::null_halfedge())
      valid = (face(halfedge(fbegin, g), g) == fbegin);
    else
      valid = false;

    if(!valid)
    {
      verr << "face " << f << " halfedge incident to face is the null halfedge." << std::endl;
      verr << "Face Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // cycle-around-face test.
    halfedge_descriptor h = halfedge( fbegin, g);
    if(h != boost::graph_traits<Graph>::null_halfedge())
    {
      halfedge_descriptor ge = h;
      do
      {
        ++n;
        h = next(h, g);
        valid = (n <= num_h && n != 0);
        if(!valid)
        {
          verr << "face " << f << " too many halfedges around face." << std::endl;
          verr << "Face Graph Structure is NOT VALID." << std::endl;
          return false;
        }
      }
      while(h != ge);
    }

    ++f;
  }

  valid = (f == num_f);
  if(!valid)
  {
    verr << "counting faces failed." << std::endl;
    verr << "Face Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  BOOST_FOREACH(halfedge_descriptor i, halfedges(g))
  {
    ++hn;

    //counting borders
    if(is_border(i, g))
      ++nb;

    // face integrity.
    valid = (face(i, g) == face(next(i, g), g));
    if(!valid)
    {
      verr << "halfedge " << hn << " face(hd) != face(next(hd))." << std::endl;
      verr << "Face Graph Structure is NOT VALID." << std::endl;
      return false;
    }
  }

  valid = (n + nb == num_h);
  if(!valid)
  {
    verr << "sum border halfedges (2*nb) = " << 2 * nb << std::endl;
    verr << "counting halfedges via faces failed." << std::endl;
    verr << "Face Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  valid = (f == num_f);
  if(!valid)
    verr << "counting faces failed." << std::endl;

  verr << "Face Graph Structure is " << (valid ? "valid." : "NOT VALID.") << std::endl;

  return valid;
}

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of `g`.
 *
 * `g` is valid if it is a valid `FaceListGraph` and it has distinct faces on each side of an edge.
 * calls `is_valid_face_graph()`.
 *
 * \param g the `Mesh` to test.
 * \param verb : if `true`, the details of the check will be written in the standard output.
 *
 * \tparam Mesh a model of `FaceListGraph` and `HalfedgeListGraph`, and follows
 * the definition of a \ref PMPDef "PolygonMesh"
 * \return `true` if `g` is valid, `false` otherwise.
 *
 * \see `is_valid_face_graph()`
 * \see `is_valid_halfedge_graph()`
 *
 */
template <typename Mesh>
bool is_valid_polygon_mesh(const Mesh& g, bool verb = false)
{
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;

  Verbose_ostream verr(verb);
  bool valid = is_valid_face_graph(g, verb);
  if(!valid)
    return false;

  // test for 2-manifoldness
  // Distinct facets on each side of an halfedge.
  BOOST_FOREACH(halfedge_descriptor i, halfedges(g))
  {
    valid = (face(i, g) != face(opposite(i, g), g));
    if(!valid)
    {
      verr << "both incident facets are equal." << std::endl;
      verr << "Polygon Mesh Structure is NOT VALID." << std::endl;
      return false;
    }

    valid = (next(next(i, g), g) != i);
    valid = valid && (target(i, g) != target(next(i, g), g));
    valid = valid && (target(i, g) != target(next(next(i, g), g), g));
    if(!valid)
    {
      verr << "incident facet is not at least a triangle." << std::endl;
      verr << "Polygon Mesh Structure is NOT VALID." << std::endl;
      return false;
    }
  }

  verr << "Polygon Mesh Structure is valid." << std::endl;
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
 * \brief Creates an isolated triangle
 * with its vertices initialized to `p0`, `p1` and `p2`, and adds it to the graph `g`.
 * \returns the non-border halfedge that has the target vertex associated with `p0`.
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

namespace internal {

template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_quad(typename boost::graph_traits<Graph>::vertex_descriptor v0,
          typename boost::graph_traits<Graph>::vertex_descriptor v1, 
          typename boost::graph_traits<Graph>::vertex_descriptor v2,
          typename boost::graph_traits<Graph>::vertex_descriptor v3, Graph& g)
{ 
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor face_descriptor;
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

} // namespace internal

/** 
 * \ingroup PkgBGLHelperFct
 * \brief Creates an isolated quad with
 * its vertices initialized to `p0`, `p1`, `p2`, and `p3`, and adds it to the graph `g`.
 * \returns the non-border halfedge that has the target vertex associated with `p0`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_quad(const P& p0, const P& p1, const P& p2, const P& p3, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
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
  return internal::make_quad(v0, v1, v2, v3, g);
}

/** 
 * \ingroup PkgBGLHelperFct
 * \brief Creates an isolated hexahedron
 * with its vertices initialized to `p0`, `p1`, ...\ , and `p7`, and adds it to the graph `g`.
 * \returns the halfedge that has the target vertex associated with `p0`, in the face with the vertices with the points `p0`, `p1`, `p2`, and `p3`.
 **/ 
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_hexahedron(const P& p0, const P& p1, const P& p2, const P& p3,
                const P& p4, const P& p5, const P& p6, const P& p7, Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;

  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  Point_property_map ppmap = get(CGAL::vertex_point, g);
  vertex_descriptor v0, v1, v2, v3, v4, v5, v6, v7;
  v0 = add_vertex(g);
  v1 = add_vertex(g);
  v2 = add_vertex(g);
  v3 = add_vertex(g);
  v4 = add_vertex(g);
  v5 = add_vertex(g);
  v6 = add_vertex(g);
  v7 = add_vertex(g);
  ppmap[v0] = p0;
  ppmap[v1] = p1;
  ppmap[v2] = p2;
  ppmap[v3] = p3;
  ppmap[v4] = p4;
  ppmap[v5] = p5;
  ppmap[v6] = p6;
  ppmap[v7] = p7;

  halfedge_descriptor ht = internal::make_quad(v7, v4, v5, v6, g);
  halfedge_descriptor hb = prev(internal::make_quad(v1, v0, v3, v2, g),g);
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
  return next(next(hb,g),g);
}
/** 
 * \ingroup PkgBGLHelperFct
 * \brief Creates an isolated tetrahedron
 * with its vertices initialized to `p0`, `p1`, `p2`, and `p3`, and adds it to the graph `g`.
 * \returns the halfedge that has the target vertex associated with `p0`, in the face with the vertices with the points `p0`, `p1`, and `p2`.
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
  v2 = add_vertex(g); // this and the next line are switched to keep points in order
  v1 = add_vertex(g);
  v3 = add_vertex(g);

  ppmap[v0] = p0;
  ppmap[v1] = p2;// this and the next line are switched to reorient the surface
  ppmap[v2] = p1;
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
  set_next(h0, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
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

/// \cond SKIP_IN_DOC
template <class Traits, class TriangleMesh, class VertexPointMap>
bool is_degenerate_triangle_face(
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
  TriangleMesh& tmesh,
  const VertexPointMap& vpmap,
  const Traits& traits)
{
  CGAL_assertion(!is_border(hd, tmesh));

  const typename Traits::Point_3& p1 = get(vpmap, target( hd, tmesh) );
  const typename Traits::Point_3& p2 = get(vpmap, target(next(hd, tmesh), tmesh) );
  const typename Traits::Point_3& p3 = get(vpmap, source( hd, tmesh) );
  return traits.collinear_3_object()(p1, p2, p3);
}

template <class Traits, class TriangleMesh, class VertexPointMap>
bool is_degenerate_triangle_face(
  typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
  TriangleMesh& tmesh,
  const VertexPointMap& vpmap,
  const Traits& traits)
{
  return is_degenerate_triangle_face(halfedge(fd,tmesh), tmesh, vpmap, traits);
}
/// \endcond

/**
 * \ingroup PkgBGLHelperFct
 * \brief Creates a triangulated regular prism, outward oriented,
 * having `nb_vertices` vertices in each of its bases and adds it to the graph `g`.
 * If `center` is (0, 0, 0), then the first point of the prism is (`radius`, `height`, 0)
 * \param nb_vertices the number of vertices per base. It must be greater than or equal to 3.
 * \param g the graph in which the regular prism will be created.
 * \param base_center the center of the circle in which the lower base is inscribed.
 * \param height the distance between the two bases.
 * \param radius the radius of the circles in which the bases are inscribed.
 * \param is_closed determines if the bases must be created or not. If `is_closed` is `true`, `center` is a vertex.
 * \returns the halfedge that has the target vertex associated with the first point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_regular_prism(
    typename boost::graph_traits<Graph>::vertices_size_type nb_vertices,
    Graph& g,
    const P& base_center = P(0,0,0),
    typename CGAL::Kernel_traits<P>::Kernel::FT height = 1.0,
    typename CGAL::Kernel_traits<P>::Kernel::FT radius = 1.0,
    bool is_closed = true)
{
  CGAL_assertion(nb_vertices >= 3);
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename CGAL::Kernel_traits<P>::Kernel::FT FT;

  const FT to_rad = CGAL_PI / 180.0;
  const FT precision = 360.0/nb_vertices;
  const FT diameter = 2*radius;
  Point_property_map vpmap = get(CGAL::vertex_point, g);
  std::vector<vertex_descriptor> vertices;
  vertices.resize(nb_vertices*2);
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices*2; ++i)
    vertices[i] = add_vertex(g);

  //fill vertices
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i < nb_vertices; ++i)
  {
    put(vpmap,
        vertices[i],
        P(0.5*diameter*cos(i*precision*to_rad)+base_center.x(),
          height+base_center.y(),
          -0.5*diameter*sin(i*precision*to_rad) + base_center.z()));

    put(vpmap,
        vertices[i+nb_vertices],
        P(0.5*diameter*cos(i*precision*to_rad)+base_center.x(),
          base_center.y(),
          -0.5*diameter*sin(i*precision*to_rad)+base_center.z()));
  }
  std::vector<vertex_descriptor> face;
  face.resize(3);
  //fill faces
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
  {
    face[0] = vertices[(i+1)%(nb_vertices)];
    face[1] = vertices[i];
    face[2] = vertices[(i+1)%(nb_vertices) + nb_vertices];
    Euler::add_face(face, g);

    face[0] = vertices[(i+1)%(nb_vertices) + nb_vertices];
    face[1] = vertices[i];
    face[2] = vertices[i + nb_vertices];
    Euler::add_face(face, g);
  }

  //close
  if(is_closed)
  {
    //add the base_center of the fans
    vertex_descriptor top = add_vertex(g);
    vertex_descriptor bot = add_vertex(g);
    put(vpmap, top, P(base_center.x(),height+base_center.y(),base_center.z()));
    put(vpmap, bot, P(base_center.x(),base_center.y(),base_center.z()));

    //add the faces
    for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
    {
      face[0] = vertices[i];
      face[1] = vertices[(i+1)%(nb_vertices)];
      face[2] = top;
      Euler::add_face(face, g);

      face[0] = bot;
      face[1] = vertices[(i+1)%(nb_vertices) + nb_vertices];
      face[2] = vertices[i + nb_vertices];
      Euler::add_face(face, g);
    }
  }
  return halfedge(vertices[0], vertices[1], g).first;
}

/**
 * \ingroup PkgBGLHelperFct
 * \brief Creates a pyramid, outward oriented, having `nb_vertices` vertices in its base and adds it to the graph `g`.
 *
 * If `center` is (0, 0, 0), then the first point of the base is (`radius`, 0`, 0)
 * \param nb_vertices the number of vertices in the base. It must be greater than or equal to 3.
 * \param g the graph in which the pyramid will be created
 * \param base_center the center of the circle in which the base is inscribed.
 * \param height the distance between the base and the apex.
 * \param radius the radius of the circle in which the base is inscribed.
 * \param is_closed determines if the base must be created or not. If `is_closed` is `true`, `center` is a vertex.
 * \returns the halfedge that has the target vertex associated with the apex point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_pyramid(
    typename boost::graph_traits<Graph>::vertices_size_type nb_vertices,
    Graph& g,
    const P& base_center = P(0,0,0),
    typename CGAL::Kernel_traits<P>::Kernel::FT height = 1.0,
    typename CGAL::Kernel_traits<P>::Kernel::FT radius = 1.0,
    bool is_closed = true)
{
  CGAL_assertion(nb_vertices >= 3);
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename CGAL::Kernel_traits<P>::Kernel::FT FT;
  const FT to_rad = CGAL_PI / 180.0;
  const FT precision = 360.0/nb_vertices;
  const FT diameter = 2*radius;
  Point_property_map vpmap = get(CGAL::vertex_point, g);
  std::vector<vertex_descriptor> vertices;
  vertices.resize(nb_vertices);
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0;
      i<nb_vertices; ++i)
    vertices[i] = add_vertex(g);
  vertex_descriptor apex = add_vertex(g);

  //fill vertices
  put(vpmap,
      apex,
      P(base_center.x(),
        base_center.y() + height,
        base_center.z()));
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0;
      i < nb_vertices; ++i)
  {

    put(vpmap,
        vertices[i],
        P(0.5*diameter*cos(i*precision*to_rad)+base_center.x(),
          base_center.y(),
          -0.5*diameter*sin(i*precision*to_rad)+base_center.z()));
  }
  std::vector<vertex_descriptor> face;
  face.resize(3);
  //fill faces
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0;
      i<nb_vertices; ++i)
  {
    face[0] = apex;
    face[1] = vertices[i];
    face[2] = vertices[(i+1)%(nb_vertices)];
    Euler::add_face(face, g);
  }

  //close
  if(is_closed)
  {
    //add the center of the fan
    vertex_descriptor bot = add_vertex(g);
    put(vpmap, bot, P(base_center.x(),base_center.y(),base_center.z()));

    //add the faces
    for(typename boost::graph_traits<Graph>::vertices_size_type i=0;
        i<nb_vertices; ++i)
    {
      face[0] = bot;
      face[1] = vertices[(i+1)%(nb_vertices)];
      face[2] = vertices[i];
      Euler::add_face(face, g);
    }
  }
  return halfedge(vertices[0], apex, g).first;
}

/**
 * \ingroup PkgBGLHelperFct
 * \brief Creates an icosahedron, outward oriented, centered in `center` and adds it to the graph `g`.
 * \param g the graph in which the icosahedron will be created.
 * \param center the center of the sphere in which the icosahedron is inscribed.
 * \param radius the radius of the sphere in which the icosahedron is inscribed.
 * \returns the halfedge that has the target vertex associated with the first point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_icosahedron(
    Graph& g,
    const P& center = P(0,0,0),
    typename CGAL::Kernel_traits<P>::Kernel::FT radius = 1.0)
{
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  Point_property_map vpmap = get(CGAL::vertex_point, g);
  // create the initial icosahedron
  std::vector<vertex_descriptor> v_vertices;
  v_vertices.resize(12);
  for(int i=0; i<12; ++i)
    v_vertices[i] = add_vertex(g);
  typename CGAL::Kernel_traits<P>::Kernel::FT t =
      (radius + radius*CGAL::approximate_sqrt(5.0)) / 2.0;

  put(vpmap, v_vertices[0],P(-radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[1],P( radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[2],P(-radius + center.x(), -t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[3],P( radius + center.x(), -t + center.y(), 0.0 + center.z()));

  put(vpmap, v_vertices[4],P( 0.0 + center.x(), -radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[5],P( 0.0 + center.x(),  radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[6],P( 0.0 + center.x(), -radius + center.y(), -t + center.z()));
  put(vpmap, v_vertices[7],P( 0.0 + center.x(),  radius + center.y(), -t + center.z()));

  put(vpmap, v_vertices[8],P(  t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[9],P(  t + center.x(), 0.0 + center.y(),  radius + center.z()));
  put(vpmap, v_vertices[10],P(-t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[11],P(-t + center.x(), 0.0 + center.y(),  radius + center.z()));

  std::vector<vertex_descriptor> face;
  face.resize(3);
  face[1] = v_vertices[0]; face[0] = v_vertices[5]; face[2] = v_vertices[11];
  Euler::add_face(face, g);
  face[1] = v_vertices[0]; face[0] = v_vertices[1]; face[2] = v_vertices[5];
  Euler::add_face(face, g);
  face[1] = v_vertices[0]; face[0] = v_vertices[7]; face[2] = v_vertices[1];
  Euler::add_face(face, g);
  face[1] = v_vertices[0]; face[0] = v_vertices[10]; face[2] = v_vertices[7];
  Euler::add_face(face, g);
  face[1] = v_vertices[0]; face[0] = v_vertices[11]; face[2] = v_vertices[10];
  Euler::add_face(face, g);

  face[1] = v_vertices[1] ; face[0] = v_vertices[9] ; face[2] = v_vertices[5];
  Euler::add_face(face, g);
  face[1] = v_vertices[5] ; face[0] = v_vertices[4]; face[2] = v_vertices[11];
  Euler::add_face(face, g);
  face[1] = v_vertices[11]; face[0] = v_vertices[2]; face[2] = v_vertices[10];
  Euler::add_face(face, g);
  face[1] = v_vertices[10]; face[0] = v_vertices[6] ; face[2] = v_vertices[7];
  Euler::add_face(face, g);
  face[1] = v_vertices[7] ; face[0] = v_vertices[8] ; face[2] = v_vertices[1];
  Euler::add_face(face, g);

  face[1] = v_vertices[3] ; face[0] = v_vertices[4] ; face[2] = v_vertices[9];
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[3] ; face[0] = v_vertices[2] ; face[2] = v_vertices[4];
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[3] ; face[0] = v_vertices[6] ; face[2] = v_vertices[2];
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[3] ; face[0] = v_vertices[8] ; face[2] = v_vertices[6];
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[3] ; face[0] = v_vertices[9] ; face[2] = v_vertices[8];
  Euler::add_face(face, g);                                                  
                                                                             
  face[1] = v_vertices[4] ; face[0] = v_vertices[5] ; face[2] = v_vertices[9] ;
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[2] ; face[0] = v_vertices[11] ; face[2] = v_vertices[4];
  Euler::add_face(face, g);                                                   
  face[1] = v_vertices[6] ; face[0] = v_vertices[10] ; face[2] = v_vertices[2];
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[8] ; face[0] = v_vertices[7] ; face[2] = v_vertices[6] ;
  Euler::add_face(face, g);                                                  
  face[1] = v_vertices[9] ; face[0] = v_vertices[1] ; face[2] = v_vertices[8] ;
  Euler::add_face(face, g);

  return halfedge(v_vertices[1], v_vertices[0], g).first;
}


/*!
 * \ingroup PkgBGLHelperFct
 *
 * \brief Creates a row major ordered grid with `i` cells along the width and `j` cells
 * along the height and adds it to the graph `g`.
 *
 * \param i the number of cells along the width.
 * \param j the number of cells along the height.
 * \param g the graph in which the grid will be created.
 * \param calculator the functor that will assign coordinates to the grid vertices.
 * \param triangulated decides if a cell is composed of one quad or two triangles.
 * If `triangulated` is `true`, the diagonal of each cell is oriented from (0,0) to (1,1)
 * in the cell coordinates.
 *
 * \tparam CoordinateFunctor that takes two `boost::graph_traits<Graph>::%vertices_size_type`
 * and outputs a `boost::property_traits<boost::property_map<Graph,CGAL::vertex_point_t>::%type>::%value_type`.
 * <p>%Default: a point with positive integer coordinates (`w`, `h`, 0), with `w` in [0..`i`] and `h` in [0..`j`]
 * \returns the non-border non-diagonal halfedge that has the target vertex associated with the first point of the grid (default is (0,0,0) ).
 */
#ifndef DOXYGEN_RUNNING
template<class Graph, class CoordinateFunctor>
#else
template<class Graph, class CoordinateFunctor = CGAL::Creator_uniform_3<
           typename boost::graph_traits<Graph>::vertices_size_type,
           typename boost::property_traits<typename boost::property_map<Graph, vertex_point_t>::type>::value_type> >
#endif
typename boost::graph_traits<Graph>::halfedge_descriptor
make_grid(typename boost::graph_traits<Graph>::vertices_size_type i,
          typename boost::graph_traits<Graph>::vertices_size_type j,
          Graph& g,
          const CoordinateFunctor& calculator,
          bool triangulated = false)
{
  typedef typename boost::property_map<Graph,vertex_point_t>::type Point_property_map;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typename boost::graph_traits<Graph>::vertices_size_type w(i+1), h(j+1);
  Point_property_map vpmap = get(CGAL::vertex_point, g);
  // create the initial icosahedron
  //create the vertices
  std::vector<vertex_descriptor> v_vertices;
  v_vertices.resize(static_cast<std::size_t>(w*h));
  for(std::size_t k = 0; k < v_vertices.size(); ++k)
    v_vertices[k] = add_vertex(g);
  //assign the coordinates
  for(typename boost::graph_traits<Graph>::vertices_size_type a = 0; a<w; ++a)
  {
    for(typename boost::graph_traits<Graph>::vertices_size_type b=0; b<h; ++b)
    {
      put(vpmap, v_vertices[a+w*b], calculator(a,b,0));
    }
  }

  //create the faces
  std::vector<vertex_descriptor> face;
  if(triangulated)
    face.resize(3);
  else
    face.resize(4);
  for(typename boost::graph_traits<Graph>::vertices_size_type a = 0; a<w-1; ++a)
  {
    for(typename boost::graph_traits<Graph>::vertices_size_type b = 0; b<h-1; ++b)
    {
      if(triangulated)
      {
        face[0] = v_vertices[w*b+a];
        face[1] = v_vertices[w*b+a+1];
        face[2] = v_vertices[w*(b+1)+a];
        Euler::add_face(face, g);
        face[0] = v_vertices[w*b+a+1];
        face[1] = v_vertices[w*(b+1)+a+1];
        face[2] = v_vertices[w*(b+1)+a];
        Euler::add_face(face, g);
      }
      else
      {
        face[0] = v_vertices[w*b+ a];
        face[1] = v_vertices[w*b+ a+1];
        face[2] = v_vertices[w*(b+1)+ a+1];
        face[3] = v_vertices[w*(b+1)+ a];
        Euler::add_face(face, g);
      }
    }
  }
  return halfedge(v_vertices[1], v_vertices[0], g).first;
}

//default Functor
template<class Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_grid(typename boost::graph_traits<Graph>::vertices_size_type w,
          typename boost::graph_traits<Graph>::vertices_size_type h,
          Graph& g,
          bool triangulated = false)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type Size_type;
  typedef typename boost::property_traits<typename boost::property_map<Graph, vertex_point_t>::type>::value_type Point;
  return make_grid(w, h, g, CGAL::Creator_uniform_3<Size_type, Point>(), triangulated);
}

namespace internal {

template<typename FaceGraph>
inline
typename boost::enable_if<Has_member_clear<FaceGraph>, void>::type
clear_impl(FaceGraph& g)
{ g.clear(); }

template<typename FaceGraph>
inline
typename boost::disable_if<Has_member_clear<FaceGraph>, void>::type
clear_impl(FaceGraph& g)
{
  while(boost::begin(edges(g))!=boost::end(edges(g)))
    remove_edge(*boost::begin(edges(g)), g);
  while(boost::begin(faces(g))!=boost::end(faces(g)))
    remove_face(*boost::begin(faces(g)), g);
  while(boost::begin(vertices(g))!=boost::end(vertices(g)))
    remove_vertex(*boost::begin(vertices(g)), g);
}

template <class FaceGraph>
void swap_vertices(
  typename boost::graph_traits<FaceGraph>::vertex_descriptor& p,
  typename boost::graph_traits<FaceGraph>::vertex_descriptor& q,
  FaceGraph& g)
{
 typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor hq=halfedge(q, g);
  halfedge_descriptor hp=halfedge(p, g);
  BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(hq, g))
    set_target(h, p, g);
  BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(hp, g))
    set_target(h, q, g);
  set_halfedge(p, hq, g);
  set_halfedge(q, hp, g);
}

template <class FaceGraph>
void swap_edges(
  const typename boost::graph_traits<FaceGraph>::halfedge_descriptor& h1,
  const typename boost::graph_traits<FaceGraph>::halfedge_descriptor& h2,
  FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  const halfedge_descriptor oh1 = opposite(h1, g), oh2 = opposite(h2, g);

  // backup vertex pointers
  vertex_descriptor s1 = target(oh1, g), s2 = target(oh2, g);
  vertex_descriptor t1 = target(h1, g), t2 = target(h2, g);

  // backup face pointers
  face_descriptor f1 = face(h1, g), f2 = face(h2, g);
  face_descriptor fo1 = face(oh1, g), fo2 = face(oh2, g);

  // backup next prev pointers
  halfedge_descriptor nh1 = next(h1, g), nh2 = next(h2, g);
  halfedge_descriptor ph1 = prev(h1, g), ph2 = prev(h2, g);
  halfedge_descriptor noh1 = next(oh1, g), noh2 = next(oh2, g);
  halfedge_descriptor poh1 = prev(oh1, g), poh2 = prev(oh2, g);

  // handle particular cases where next/prev are halfedges to be swapt
  if (nh1 == oh2) nh1 = oh1;
  if (nh1 == h2) nh1 = h1;
  if (nh2 == oh1) nh2 = oh2;
  if (nh2 == h1) nh2 = h2;
  if (ph1 == oh2) ph1 = oh1;
  if (ph1 == h2) ph1 = h1;
  if (ph2 == oh1) ph2 = oh2;
  if (ph2 == h1) ph2 = h2;
  if (noh1 == oh2) noh1 = oh1;
  if (noh1 == h2) noh1 = h1;
  if (noh2 == oh1) noh2 = oh2;
  if (noh2 == h1) noh2 = h2;
  if (poh1 == oh2) poh1 = oh1;
  if (poh1 == h2) poh1 = h1;
  if (poh2 == oh1) poh2 = oh2;
  if (poh2 == h1) poh2 = h2;

  // (1) exchange next pointers
  set_next(h1, nh2, g);
  set_next(h2, nh1, g);
  set_next(ph1, h2, g);
  set_next(ph2, h1, g);
  set_next(oh1, noh2, g);
  set_next(oh2, noh1, g);
  set_next(poh1, oh2, g);
  set_next(poh2, oh1, g);

  // (2) exchange vertex-halfedge pointers
  set_target(h1, t2, g);
  set_target(h2, t1, g);
  set_target(oh1, s2, g);
  set_target(oh2, s1, g);
  if (halfedge(t1, g)==h1) set_halfedge(t1, h2, g);
  if (halfedge(t2, g)==h2) set_halfedge(t2, h1, g);
  if (halfedge(s1, g)==oh1) set_halfedge(s1, oh2, g);
  if (halfedge(s2, g)==oh2) set_halfedge(s2, oh1, g);

  // (3) exchange face-halfedge pointers
  set_face(h1, f2, g);
  set_face(h2, f1, g);
  set_face(oh1, fo2, g);
  set_face(oh2, fo1, g);

  face_descriptor nf = boost::graph_traits<FaceGraph>::null_face();
  if (f1 != nf && halfedge(f1, g)==h1) set_halfedge(f1, h2, g);
  if (f2 != nf && halfedge(f2, g)==h2) set_halfedge(f2, h1, g);
  if (fo1 != nf && halfedge(fo1, g)==oh1) set_halfedge(fo1, oh2, g);
  if (fo2 != nf && halfedge(fo2, g)==oh2) set_halfedge(fo2, oh1, g);
}


} //end of internal namespace

/**
 * \ingroup PkgBGLHelperFct
 *
 * removes all vertices, faces and halfedges from a graph. Calls
 * `remove_edge()`, `remove_vertex()`, and `remove_face()` for each
 * edge, vertex or face.
 *
 * If the graph has a member function `clear()`, it will be called
 * instead.
 * 
 * @tparam FaceGraph model of `MutableHalfedgeGraph` and `MutableFaceGraph`
 *
 * @param g the graph to clear
 *
 **/
template<typename FaceGraph>
void clear(FaceGraph& g)
{ 
  internal::clear_impl(g);
  CGAL_postcondition(std::distance(boost::begin(edges(g)),boost::end(edges(g))) == 0);
  CGAL_postcondition(std::distance(boost::begin(vertices(g)),boost::end(vertices(g))) == 0);
  CGAL_postcondition(std::distance(boost::begin(faces(g)),boost::end(faces(g))) == 0);
}

/**
* \ingroup PkgBGLHelperFct
*
* checks whether the graph is empty, by checking that it does not contain any vertex.
*
* @tparam FaceGraph model of `FaceGraph`
*
* @param g the graph to test
*
**/
template<typename FaceGraph>
bool is_empty(const FaceGraph& g)
{
  return boost::empty(vertices(g));
}

} // namespace CGAL

// Include "Euler_operations.h" at the end, because its implementation
// requires this header.
#include <CGAL/boost/graph/Euler_operations.h>

#endif // CGAL_BOOST_GRAPH_HELPERS_H
