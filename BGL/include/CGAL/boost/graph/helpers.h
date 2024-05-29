// Copyright (c) 2014 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_HELPERS_H
#define CGAL_BOOST_GRAPH_HELPERS_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/internal/Has_member_clear.h>
#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/function_objects.h>
#include <CGAL/IO/Verbose_ostream.h>

#include <boost/range/empty.hpp>

#include <type_traits>

namespace CGAL {

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
std::optional<typename boost::graph_traits<FaceGraph>::halfedge_descriptor>
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
  return std::optional<typename boost::graph_traits<FaceGraph>::halfedge_descriptor>();
}

namespace BGL {

template <typename Graph>
bool is_valid_vertex_descriptor(typename boost::graph_traits<Graph>::vertex_descriptor v,
                                const Graph& g,
                                const bool verb = false)
{
  Verbose_ostream verr(verb);
  bool valid = true;

  // null vertex
  valid = (v != boost::graph_traits<Graph>::null_vertex());
  if(!valid)
  {
    verr << "vertex is null." << std::endl;
    return false;
  }

  if(!CGAL::internal::is_isolated(v, g))
  {
    // Incident halfedge integrity
    valid = (target(halfedge(v, g), g) == v);
    if(!valid)
    {
      verr << "vertex has invalid halfedge()." << std::endl;
      return false;
    }
  }

  return true;
}

template <typename Graph>
bool is_valid_halfedge_descriptor(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                                  const Graph& g,
                                  const bool verb = false)
{
  Verbose_ostream verr(verb);
  bool valid = true;

  // null halfedge
  valid = (h != boost::graph_traits<Graph>::null_halfedge());
  if(!valid)
  {
    verr << "halfedge is null." << std::endl;
    return false;
  }

  // Pointer integrity.
  valid = (prev(h, g) != boost::graph_traits<Graph>::null_halfedge());
  valid = valid && (next(h, g) != boost::graph_traits<Graph>::null_halfedge());
  valid = valid && (opposite(h, g) != boost::graph_traits<Graph>::null_halfedge());
  if(!valid)
  {
    verr << "halfedge's prev / next / opposite halfedges are null." << std::endl;
    return false;
  }

  // degeneracies
  valid = (next(h, g) != h);
  valid = valid && (prev(h, g) != h);
  valid = valid && (opposite(h, g) != h);
  valid = valid && (target(h, g) != target(opposite(h, g), g));
  if(!valid)
  {
    verr << "combinatorial degeneracies." << std::endl;
    return false;
  }

  // edge integrity
  valid = (halfedge(edge(h, g), g) == h);
  if(!valid)
  {
    verr << "halfedge has an invalid edge." << std::endl;
    return false;
  }

  // opposite integrity.
  valid = (opposite(h, g) != h);
  valid = valid && (opposite(opposite(h, g), g) == h);
  if(!valid)
  {
    verr << "halfedge has invalid opposite()." << std::endl;
    return false;
  }

  // previous integrity.
  valid = (prev(next(h, g), g) == h);
  valid = valid && (next(prev(h, g), g) == h);
  if(!valid)
  {
    verr << "prev(next(hd)) != hd OR next(prev(hd)) != hd" << std::endl;
    return false;
  }

  // vertex integrity.
  valid = (target(h, g) != boost::graph_traits<Graph>::null_vertex());
  if(!valid)
  {
    verr << "target of halfedge is the null vertex." << std::endl;
    return false;
  }

  valid = (target(h, g) == target(opposite(next(h, g), g), g));
  valid = valid && (target(opposite(h, g), g) == target(prev(h, g), g));
  if(!valid)
  {
    verr << "vertex inconsistencies with prev/next." << std::endl;
    return false;
  }

  return true;
}

template <typename FaceGraph>
bool is_valid_edge_descriptor(typename boost::graph_traits<FaceGraph>::edge_descriptor e,
                              const FaceGraph& g,
                              const bool verb = false)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  Verbose_ostream verr(verb);
  bool valid = true;

  // there is no null_edge() in the Graph concepts

  // Pointer integrity.
  const halfedge_descriptor h = halfedge(e, g);
  valid = (h != boost::graph_traits<FaceGraph>::null_halfedge());
  if(!valid)
  {
    verr << "halfedge incident to edge is the null halfedge." << std::endl;
    return false;
  }

  // halfedge integrity
  valid = (edge(h, g) == e);
  if(!valid)
  {
    verr << "edge has an invalid halfedge()." << std::endl;
    return false;
  }

  return true;
}

template <typename FaceGraph>
bool is_valid_face_descriptor(typename boost::graph_traits<FaceGraph>::face_descriptor f,
                              const FaceGraph& g,
                              const bool verb = false)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  Verbose_ostream verr(verb);
  bool valid = true;

  // null face
  valid = (f != boost::graph_traits<FaceGraph>::null_face());
  if(!valid)
  {
    verr << "face is null." << std::endl;
    return false;
  }

  // Pointer integrity.
  const halfedge_descriptor h = halfedge(f, g);
  valid = (h != boost::graph_traits<FaceGraph>::null_halfedge());
  if(!valid)
  {
    verr << "halfedge incident to face is the null halfedge." << std::endl;
    return false;
  }

  valid = (face(h, g) == f);
  if(!valid)
  {
    verr << "face has an invalid halfedge()." << std::endl;
    return false;
  }

  // face integrity.
  valid = (face(h, g) == face(next(h, g), g));
  valid = valid && (face(h, g) == face(prev(h, g), g));
  if(!valid)
  {
    verr << "different face incident to face halfedges." << std::endl;
    return false;
  }

  return true;
}

} // namespace BGL

// These empty functions simply calling the BGL versions (just above) are done such that
// a specific graph type (e.g. Surface_mesh) can overload those and still call the BGL versions
// without duplicating code
template <typename Graph>
bool is_valid_vertex_descriptor(typename boost::graph_traits<Graph>::vertex_descriptor v,
                                const Graph& g,
                                const bool verb = false)
{
  return BGL::is_valid_vertex_descriptor(v, g, verb);
}

template <typename Graph>
bool is_valid_halfedge_descriptor(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                                  const Graph& g,
                                  const bool verb = false)
{
  return BGL::is_valid_halfedge_descriptor(h, g, verb);
}

template <typename Graph>
bool is_valid_edge_descriptor(typename boost::graph_traits<Graph>::edge_descriptor e,
                              const Graph& g,
                              const bool verb = false)
{
  return BGL::is_valid_edge_descriptor(e, g, verb);
}

template <typename Graph>
bool is_valid_face_descriptor(typename boost::graph_traits<Graph>::face_descriptor f,
                              const Graph& g,
                              const bool verb = false)
{
  return BGL::is_valid_face_descriptor(f, g, verb);
}

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of the graph `g`.
 *
 * The graph `g` is valid if it follows the rules of the `HalfedgeListGraph` concept
 * and all of its associations are reciprocal (for example, `prev(next(h, g), g)` must be `h`,
 * and `next(prev(h, g), g)` must be `h`).
 *
 * \param g the graph to test
 * \param verb if `true`, the details of the check will be written in the standard output.
 *
 * \tparam Graph a model of `HalfedgeListGraph`
 *
 * \return `true` if `g` is valid, `false` otherwise.
 *
 */
template<typename Graph>
bool is_valid_halfedge_graph(const Graph& g, bool verb = false)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor   halfedge_descriptor;

  Verbose_ostream verr(verb);

  std::size_t num_v = CGAL::internal::exact_num_vertices(g),
              num_e = CGAL::internal::exact_num_edges(g),
              num_h = CGAL::internal::exact_num_halfedges(g);

  bool valid = (1 != (num_h&1) && (2*num_e == num_h));
  if(!valid)
  {
    verr << "number of halfedges is odd." << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  // All halfedges.
  std::size_t hc = 0;
  for(halfedge_descriptor h : halfedges(g))
  {
    if(!is_valid_halfedge_descriptor(h, g, verb))
    {
      verr << "halfedge " << hc << " is invalid." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    ++hc;
  }

  valid = (hc == num_h);
  if(!valid)
  {
    verr << "counting halfedges failed: " << hc << " vs " << num_h << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  // All vertices.
  std::size_t vc = 0;
  hc = 0;
  for(vertex_descriptor v : vertices(g))
  {
    if(!is_valid_vertex_descriptor(v, g, verb))
    {
      verr << "vertex " << vc << " is invalid." << std::endl;
      verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // cycle-around-vertex test.
    if(!CGAL::internal::is_isolated(v, g))
    {
      halfedge_descriptor h = halfedge(v, g), done = h;
      do
      {
        ++hc;
        h = opposite(next(h, g), g);
        valid = (hc <= num_h);
        if(!valid)
        {
          verr << "vertex " << vc << " too many halfedges around vertex." << std::endl;
          verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
          return false;
        }
      }
      while(h != done);
    }

    ++vc;
  }

  valid = (vc == num_v);
  if(!valid)
  {
    verr << "counting vertices failed: " << vc << " vs " << num_v << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  valid = (hc == num_h);
  if(!valid)
  {
    verr << "counting halfedges via vertices failed: " << hc << " vs " << num_h << std::endl;
    verr << "Halfedge Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  verr << "Halfedge Graph Structure is valid" << std::endl;

  return valid;
}

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of the graph `g`.
 *
 * The graph `g` is a valid face graph if it is a valid halfedge graph, and if it follows the rules
 * of the `FaceListGraph` concept and all of its associations are reciprocal (for example,
 * `face(halfedge(f,g),g)` must be `f`).
 *
 * \param g the graph to test
 * \param verb if `true`, the details of the check will be written in the standard output
 *
 * \tparam FaceGraph a model of `FaceListGraph` and `HalfedgeListGraph`
 *
 * \return `true` if `g` is valid, `false` otherwise.
 *
 * \see `is_valid_halfedge_graph()`
 */
template<typename FaceGraph>
bool is_valid_face_graph(const FaceGraph& g, bool verb = false)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor       face_descriptor;

  Verbose_ostream verr(verb);

  std::size_t num_f = CGAL::internal::exact_num_faces(g),
              num_h = CGAL::internal::exact_num_halfedges(g);

  std::size_t fc = 0, hc = 0, nb = 0;

  bool valid = is_valid_halfedge_graph(g, verb);
  if(!valid)
    return false;

  // All faces.
  for(face_descriptor f : faces(g))
  {
    if(!is_valid_face_descriptor(f, g, verb))
    {
      verr << "face " << fc << " is invalid." << std::endl;
      verr << "Face Graph Structure is NOT VALID." << std::endl;
      return false;
    }

    // cycle-around-face test.
    halfedge_descriptor h = halfedge(f, g), done(h);
    do
    {
      ++hc;
      valid = (hc <= num_h);
      if(!valid)
      {
        verr << "face " << fc << " too many halfedges around face." << std::endl;
        verr << "Face Graph Structure is NOT VALID." << std::endl;
        return false;
      }
      h = next(h, g);
    }
    while(h != done);

    ++fc;
  }

  valid = (fc == num_f);
  if(!valid)
  {
    verr << "counting faces failed: " << fc << " vs " << num_f << std::endl;
    verr << "Face Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  for(halfedge_descriptor h : halfedges(g))
  {
    //counting borders
    if(is_border(h, g))
      ++nb;
  }

  valid = (hc + nb == num_h);
  if(!valid)
  {
    verr << "counting halfedges via faces failed." << std::endl;
    verr << "sum border halfedges (2*nb) = " << 2 * nb << " vs " << num_h << std::endl;
    verr << "Face Graph Structure is NOT VALID." << std::endl;
    return false;
  }

  verr << "Face Graph Structure is valid" << std::endl;

  return valid;
}

/*!
  \ingroup PkgBGLHelperFct
 * \brief checks the integrity of the mesh `g`.
 *
 * The mesh `g` is a valid polygon mesh if it is a valid face graph and if it follows the rules
 * defined in \ref PMPDef "PolygonMesh".
 *
 * \param g the `Mesh` to test
 * \param verb if `true`, the details of the check will be written in the standard output
 *
 * \tparam Mesh a model of `FaceListGraph` and `HalfedgeListGraph`
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
  for(halfedge_descriptor i : halfedges(g))
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
    returns `true` if there are no border edges.
  */
template <typename FaceGraph>
bool is_closed(const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  for(halfedge_descriptor hd : halfedges(g)){
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
  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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
  for(vertex_descriptor vd : vertices(g)){
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
  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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
  for(vertex_descriptor vd : vertices(g)){
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
bool is_isolated_triangle(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                          const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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
bool is_triangle(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                 const FaceGraph& g)
{
  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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
  for(face_descriptor fd : faces(g)){
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
bool is_isolated_quad(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                      const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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
bool is_quad(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
             const FaceGraph& g)
{
  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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

  for(face_descriptor fd : faces(g)){
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
bool is_tetrahedron(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                    const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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

  /*!
   \ingroup PkgBGLHelperFct
    returns `true` iff the connected component denoted by `hd` is a hexahedron.
  */
template <typename FaceGraph>
bool is_hexahedron(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                   const FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(hd, g));

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

namespace internal {

template <class FaceGraph>
void swap_vertices(typename boost::graph_traits<FaceGraph>::vertex_descriptor& p,
                   typename boost::graph_traits<FaceGraph>::vertex_descriptor& q,
                   FaceGraph& g)
{
 typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_vertex_descriptor(p, g) && is_valid_vertex_descriptor(q, g));

  halfedge_descriptor hq=halfedge(q, g);
  halfedge_descriptor hp=halfedge(p, g);
  for(halfedge_descriptor h : halfedges_around_target(hq, g))
    set_target(h, p, g);
  for(halfedge_descriptor h : halfedges_around_target(hp, g))
    set_target(h, q, g);
  set_halfedge(p, hq, g);
  set_halfedge(q, hp, g);
}

template <class FaceGraph>
void swap_edges(const typename boost::graph_traits<FaceGraph>::halfedge_descriptor& h1,
                const typename boost::graph_traits<FaceGraph>::halfedge_descriptor& h2,
                FaceGraph& g)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(h1, g) && is_valid_halfedge_descriptor(h2, g));

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

template <typename Graph>
void collect_garbage(Graph&)
{
  // nothing by default
}

} //end of internal namespace

/**
 * \ingroup PkgBGLHelperFct
 *
 * removes all vertices, faces and halfedges from a graph. Calls
 * \link MutableHalfedgeGraph `remove_vertex()`\endlink,
 * \link MutableHalfedgeGraph `remove_edge()`\endlink, and
 * \link MutableFaceGraph `remove_face()`\endlink, for each vertex, edge, and face.
 *
 * Note that some graphs have a specialized version of this function to improve
 * complexity.
 *
 * @warning This function does not perform anything more than what is advertised above. It is
 * up to the user to e.g. clean garbage or remove internal property maps (if relevant, and desired).
 *
 * @tparam FaceGraph model of `MutableHalfedgeGraph` and `MutableFaceGraph`
 *
 * @param g the graph whose elements will be removed
 *
 * @sa `CGAL::clear()`
 **/
template<typename FaceGraph>
void remove_all_elements(FaceGraph& g)
{
  while(std::begin(edges(g)) != std::end(edges(g)))
    remove_edge(*std::begin(edges(g)), g);
  while(std::begin(faces(g)) != std::end(faces(g)))
    remove_face(*std::begin(faces(g)), g);
  while(std::begin(vertices(g)) != std::end(vertices(g)))
    remove_vertex(*std::begin(vertices(g)), g);

  CGAL_postcondition(std::distance(std::cbegin(vertices(g)), std::cend(vertices(g))) == 0);
  CGAL_postcondition(std::distance(std::cbegin(edges(g)), std::cend(edges(g))) == 0);
  CGAL_postcondition(std::distance(std::cbegin(faces(g)), std::cend(faces(g))) == 0);
}

namespace internal {

template<typename FaceGraph>
inline
std::enable_if_t<Has_member_clear<FaceGraph>::value, void>
clear_impl(FaceGraph& g)
{
  g.clear();
}

template<typename FaceGraph>
inline
std::enable_if_t<!Has_member_clear<FaceGraph>::value, void>
clear_impl(FaceGraph& g)
{
  remove_all_elements(g);
}

} // namespace internal

/**
 * \ingroup PkgBGLHelperFct
 *
 * removes all vertices, faces and halfedges from a graph. Calls
 * \link MutableHalfedgeGraph `remove_vertex()`\endlink,
 * \link MutableHalfedgeGraph `remove_edge()`\endlink, and
 * \link MutableFaceGraph `remove_face()`\endlink, for each vertex, edge, and face.
 *
 * If the graph has a member function `clear()`, it will be called
 * instead.
 *
 * @warning If it exists, the `clear()` function of a graph might do more than
 * simply remove elements. For example, `CGAL::Surface_mesh::clear()` collects garbage
 * and removes *all* property maps added by a call to `CGAL::Surface_mesh::add_property_map()` for all simplex types.
 *
 * @tparam FaceGraph model of `MutableHalfedgeGraph` and `MutableFaceGraph`
 *
 * @param g the graph to clear
 *
 * @sa `CGAL::remove_all_elements()`
 **/
template<typename FaceGraph>
void clear(FaceGraph& g)
{
  internal::clear_impl(g);

  CGAL_postcondition(std::distance(std::cbegin(vertices(g)), std::cend(vertices(g))) == 0);
  CGAL_postcondition(std::distance(std::cbegin(edges(g)), std::cend(edges(g))) == 0);
  CGAL_postcondition(std::distance(std::cbegin(faces(g)), std::cend(faces(g))) == 0);
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

/// \ingroup PkgBGLHelperFct
///
/// \brief returns the number of calls to `next()` one has to apply to the halfedge `hd`
///        for `source(hd, mesh) == vd` to be true, starting from `hd = halfedge(fd, tm)`.
///
/// \tparam Graph a model of `FaceGraph`
///
/// \param vd a vertex of `g` whose index is sought
/// \param fd a face of `g` in which the index of `vd` is sought
/// \param g a mesh of type `Graph`
///
/// \pre `vd` is a vertex of `fd`.
template <typename Graph>
int vertex_index_in_face(const typename boost::graph_traits<Graph>::vertex_descriptor vd,
                         const typename boost::graph_traits<Graph>::face_descriptor fd,
                         const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_vertex_descriptor(vd, g) && is_valid_face_descriptor(fd, g));

  halfedge_descriptor start = halfedge(fd, g);
  halfedge_descriptor current = start;
  int counter = 0;

  do
  {
    if(source(current, g) == vd)
      break;

    ++counter;
    current = next(current, g);
  }
  while(current != start);

  if(counter != 0 && current == start)
  {
    CGAL_assertion_msg(false, "Could not find vertex in face");
    return -1;
  }

  return counter;
}

/// \ingroup PkgBGLHelperFct
///
/// \brief returns the number of calls to `next(hd, tm)` one has to apply to `hd` for `hd == he`
///        to be true, starting from `hd = halfedge(face(he, tm), tm)`.
///
/// \tparam Graph a model of `FaceGraph`.
///
/// \param he a halfedge of `g` whose index in `face(he, tm)` is sought
/// \param g an object of type `Graph`
///
template <typename Graph>
int halfedge_index_in_face(typename boost::graph_traits<Graph>::halfedge_descriptor he,
                           const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor     face_descriptor;

  CGAL_precondition(is_valid_halfedge_descriptor(he, g));
  CGAL_precondition(!is_border(he, g));

  face_descriptor f = face(he, g);
  halfedge_descriptor start = halfedge(f, g);
  halfedge_descriptor current = start;
  int count = 0;

  while(current != he)
  {
    current = next(current, g);
    ++count;
  }

  return count;
}

} // namespace CGAL


#endif // CGAL_BOOST_GRAPH_HELPERS_H
