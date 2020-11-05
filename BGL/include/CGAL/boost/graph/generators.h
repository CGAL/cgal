// Copyright (c) 2014, 2017 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno,
//             Mael Rouxel-Labbé

#ifndef CGAL_BOOST_GRAPH_GENERATORS_H
#define CGAL_BOOST_GRAPH_GENERATORS_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Random.h>
#include <CGAL/function_objects.h>

namespace CGAL {
namespace Euler {

// Some forward declaration to break the helpers.h > generators.h > Euler_operations.h cycle
template< typename Graph>
void fill_hole(typename boost::graph_traits<Graph>::halfedge_descriptor h,
               Graph& g);

template<typename Graph , typename VertexRange >
typename boost::graph_traits<Graph>::face_descriptor add_face(const VertexRange& vr,
                                                              Graph& g);

} // namespace Euler

namespace internal {

template <class FaceGraph>
void swap_vertices(typename boost::graph_traits<FaceGraph>::vertex_descriptor& p,
                   typename boost::graph_traits<FaceGraph>::vertex_descriptor& q,
                   FaceGraph& g);

template<typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
random_entity_in_range(InputIterator first, InputIterator beyond,
                       CGAL::Random& rnd = get_default_random())
{
  typedef typename std::iterator_traits<InputIterator>::difference_type  size_type;

  size_type zero = 0, ne = std::distance(first, beyond);
  std::advance(first, rnd.uniform_int(zero, ne - 1));

  return *first;
}

template<typename InputIterator>
typename std::iterator_traits<InputIterator>::value_type
random_entity_in_range(const CGAL::Iterator_range<InputIterator>& range,
                       CGAL::Random& rnd = get_default_random())
{
  return random_entity_in_range(range.begin(), range.end(), rnd);
}

// \brief returns a random non-null vertex incident to the face `fd` of the polygon mesh `g`.
// \tparam Graph a model of `HalfedgeGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::vertex_descriptor
random_vertex_in_face(typename boost::graph_traits<Graph>::face_descriptor fd,
                      const Graph& g,
                      CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(vertices_around_face(halfedge(fd, g), g), rnd);
}

// \brief returns a random non-null halfedge incident to the face `fd` of the polygon mesh `g`.
// \tparam Graph a model of `HalfedgeGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
random_halfedge_in_face(typename boost::graph_traits<Graph>::face_descriptor fd,
                        const Graph& g,
                        CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(halfedges_around_face(halfedge(fd, g), g), rnd);
}

// \brief returns a random non-null vertex of the polygon mesh `g`.
// \tparam Graph a model of `VertexListGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::vertex_descriptor
random_vertex_in_mesh(const Graph& g, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(vertices(g), rnd);
}

// \brief returns a random non-null halfedge of the polygon mesh `g`.
// \tparam Graph a model of `HalfedgeListGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
random_halfedge_in_mesh(const Graph& g, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(halfedges(g), rnd);
}

// \brief returns a random non-null edge of the polygon mesh `g`.
// \tparam Graph a model of `EdgeListGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::edge_descriptor
random_edge_in_mesh(const Graph& g, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(edges(g), rnd);
}

// \brief returns a random non-null face of the polygon mesh `g`.
// \tparam Graph a model of `FaceListGraph`
template<typename Graph>
typename boost::graph_traits<Graph>::face_descriptor
random_face_in_mesh(const Graph& g, CGAL::Random& rnd = get_default_random())
{
  return internal::random_entity_in_range(faces(g), rnd);
}

} // namespace internal

/**
 * \ingroup PkgBGLHelperFct
 *
 * \brief Creates an isolated triangle
 * with its vertices initialized to `p0`, `p1` and `p2`, and adds it to the graph `g`.
 *
 * \returns the non-border halfedge that has the target vertex associated with `p0`.
 **/
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_triangle(const P& p0, const P& p1, const P& p2, Graph& g)
{
  typedef typename boost::graph_traits<Graph>                       Traits;
  typedef typename Traits::vertex_descriptor                        vertex_descriptor;
  typedef typename Traits::halfedge_descriptor                      halfedge_descriptor;

  typedef typename Traits::face_descriptor                          face_descriptor;
  typedef typename boost::property_map<Graph,vertex_point_t>::type  Point_property_map;
  Point_property_map ppmap = get(CGAL::vertex_point, g);
  vertex_descriptor v0, v1, v2;
  v0 = add_vertex(g);
  v1 = add_vertex(g);
  v2 = add_vertex(g);

  ppmap[v0] = p0;
  ppmap[v1] = p1;
  ppmap[v2] = p2;
  halfedge_descriptor h0 = halfedge(add_edge(g), g);
  halfedge_descriptor h1 = halfedge(add_edge(g), g);
  halfedge_descriptor h2 = halfedge(add_edge(g), g);
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
  set_face(h0, f, g);
  set_face(h1, f, g);
  set_face(h2, f, g);
  set_halfedge(f, h0, g);
  h0 = opposite(h0, g);
  h1 = opposite(h1, g);
  h2 = opposite(h2, g);
  set_next(h0, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  set_face(h0, boost::graph_traits<Graph>::null_face(), g);
  set_face(h1, boost::graph_traits<Graph>::null_face(), g);
  set_face(h2, boost::graph_traits<Graph>::null_face(), g);

  return opposite(h2, g);
}

namespace internal {

template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_quad(typename boost::graph_traits<Graph>::vertex_descriptor v0,
          typename boost::graph_traits<Graph>::vertex_descriptor v1,
          typename boost::graph_traits<Graph>::vertex_descriptor v2,
          typename boost::graph_traits<Graph>::vertex_descriptor v3,
          Graph& g)
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor face_descriptor;
  halfedge_descriptor h0 = halfedge(add_edge(g), g);
  halfedge_descriptor h1 = halfedge(add_edge(g), g);
  halfedge_descriptor h2 = halfedge(add_edge(g), g);
  halfedge_descriptor h3 = halfedge(add_edge(g), g);
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
  set_face(h0, f, g);
  set_face(h1, f, g);
  set_face(h2, f, g);
  set_face(h3, f, g);
  set_halfedge(f, h0, g);
  h0 = opposite(h0, g);
  h1 = opposite(h1, g);
  h2 = opposite(h2, g);
  h3 = opposite(h3, g);
  set_next(h0, h3, g);
  set_next(h3, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  set_target(h3, v3, g);
  set_face(h0, boost::graph_traits<Graph>::null_face(), g);
  set_face(h1, boost::graph_traits<Graph>::null_face(), g);
  set_face(h2, boost::graph_traits<Graph>::null_face(), g);
  set_face(h3, boost::graph_traits<Graph>::null_face(), g);
  return opposite(h3, g);
}

// default Functor for make_grid
template<typename Size_type, typename Point>
struct Default_grid_maker
    : public CGAL::Creator_uniform_3<Size_type, Point>
{
  Point operator()(const Size_type& i, const Size_type& j) const {
    return CGAL::Creator_uniform_3<Size_type, Point>::operator ()(i,j,0);
  }
};

} // namespace internal

/**
 * \ingroup PkgBGLHelperFct
 *
 * \brief Creates an isolated quad with
 * its vertices initialized to `p0`, `p1`, `p2`, and `p3`, and adds it to the graph `g`.
 *
 * \returns the non-border halfedge that has the target vertex associated with `p0`.
 **/
template<typename Graph, typename P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_quad(const P& p0, const P& p1, const P& p2, const P& p3, Graph& g)
{
  typedef typename boost::graph_traits<Graph>                       Traits;
  typedef typename Traits::vertex_descriptor                        vertex_descriptor;
  typedef typename boost::property_map<Graph,vertex_point_t>::type  Point_property_map;

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
 * \image html hexahedron.png
 * \image latex hexahedron.png
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

  halfedge_descriptor ht = internal::make_quad(v4, v5, v6, v7, g);
  halfedge_descriptor hb = prev(internal::make_quad(v0, v3, v2, v1, g), g);
  for(int i=0; i <4; ++i)
  {
    halfedge_descriptor h = halfedge(add_edge(g), g);
    set_target(h,target(hb, g), g);
    set_next(h, opposite(hb, g), g);
    set_next(opposite(prev(ht, g), g), h, g);
    h = opposite(h, g);
    set_target(h, source(prev(ht, g), g), g);
    set_next(h, opposite(next(next(ht, g), g), g), g);
    set_next(opposite(next(hb, g), g), h, g);
    hb = next(hb, g);
    ht = prev(ht, g);
  }
  for(int i=0; i <4; ++i)
  {
    Euler::fill_hole(opposite(hb, g), g);
    hb = next(hb, g);
  }

  return next(next(hb, g), g);
}

/**
 * \ingroup PkgBGLHelperFct
 * \brief Creates an isolated tetrahedron
 * with its vertices initialized to `p0`, `p1`, `p2`, and `p3`, and adds it to the graph `g`.
 * \image html tetrahedron.png
 * \image latex tetrahedron.png
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
  halfedge_descriptor h0 = halfedge(add_edge(g), g);
  halfedge_descriptor h1 = halfedge(add_edge(g), g);
  halfedge_descriptor h2 = halfedge(add_edge(g), g);
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
  set_face(h0, f, g);
  set_face(h1, f, g);
  set_face(h2, f, g);
  set_halfedge(f, h0, g);
  h0 = opposite(h0, g);
  h1 = opposite(h1, g);
  h2 = opposite(h2, g);
  set_next(h0, h2, g);
  set_next(h2, h1, g);
  set_next(h1, h0, g);
  set_target(h0, v0, g);
  set_target(h1, v1, g);
  set_target(h2, v2, g);
  halfedge_descriptor h3 = halfedge(add_edge(g), g);
  halfedge_descriptor h4 = halfedge(add_edge(g), g);
  halfedge_descriptor h5 = halfedge(add_edge(g), g);
  set_target(h3, v3, g);
  set_target(h4, v3, g);
  set_target(h5, v3, g);
  set_halfedge(v3, h3, g);

  set_next(h0, h3, g);
  set_next(h1, h4, g);
  set_next(h2, h5, g);

  set_next(h3, opposite(h4, g), g);
  set_next(h4, opposite(h5, g), g);
  set_next(h5, opposite(h3, g), g);
  set_next(opposite(h4, g), h0, g);
  set_next(opposite(h5, g), h1, g);
  set_next(opposite(h3, g), h2, g);

  set_target(opposite(h3, g), v0, g);
  set_target(opposite(h4, g), v1, g);
  set_target(opposite(h5, g), v2, g);

  f = add_face(g);
  set_halfedge(f, h0, g);
  set_face(h0, f, g);
  set_face(h3, f, g);
  set_face(opposite(h4, g), f, g);
  f = add_face(g);
  set_halfedge(f, h1, g);
  set_face(h1, f, g);
  set_face(h4, f, g);
  set_face(opposite(h5, g), f, g);
  f = add_face(g);
  set_halfedge(f, h2, g);
  set_face(h2, f, g);
  set_face(h5, f, g);
  set_face(opposite(h3, g), f, g);

  return opposite(h2, g);
}

/**
 * \ingroup PkgBGLHelperFct
 *
 * \brief Creates a triangulated regular prism, outward oriented,
 * having `nb_vertices` vertices in each of its bases and adds it to the graph `g`.
 * If `center` is (0, 0, 0), then the first point of the prism is (`radius`, `height`, 0)
 *
 * \param nb_vertices the number of vertices per base. It must be greater than or equal to 3.
 * \param g the graph in which the regular prism will be created.
 * \param base_center the center of the circle in which the lower base is inscribed.
 * \param height the distance between the two bases.
 * \param radius the radius of the circles in which the bases are inscribed.
 * \param is_closed determines if the bases must be created or not. If `is_closed` is `true`, `center` is a vertex.
 *
 * \returns the halfedge that has the target vertex associated with the first point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_regular_prism(typename boost::graph_traits<Graph>::vertices_size_type nb_vertices,
                   Graph& g,
                   const P& base_center = P(0,0,0),
                   typename CGAL::Kernel_traits<P>::Kernel::FT height = 1.0,
                   typename CGAL::Kernel_traits<P>::Kernel::FT radius = 1.0,
                   bool is_closed = true)
{
  CGAL_assertion(nb_vertices >= 3);

  typedef typename boost::graph_traits<Graph>::vertex_descriptor     vertex_descriptor;
  typedef typename CGAL::Kernel_traits<P>::Kernel::FT                FT;

  typedef typename boost::property_map<Graph, vertex_point_t>::type  Point_property_map;
  Point_property_map vpmap = get(CGAL::vertex_point, g);

  const FT to_rad = CGAL_PI / 180.0;
  const FT precision = 360.0 / nb_vertices;
  const FT diameter = 2 * radius;

  std::vector<vertex_descriptor> vertices;
  vertices.resize(nb_vertices*2);
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices*2; ++i)
    vertices[i] = add_vertex(g);

  //fill vertices
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i < nb_vertices; ++i)
  {
    put(vpmap, vertices[i],
        P(0.5*diameter * cos(i*precision*to_rad) + base_center.x(),
          height+base_center.y(),
          -0.5*diameter * sin(i*precision*to_rad) + base_center.z()));

    put(vpmap,
        vertices[i+nb_vertices],
        P(0.5*diameter * cos(i*precision*to_rad) + base_center.x(),
          base_center.y(),
          -0.5*diameter * sin(i*precision*to_rad) + base_center.z()));
  }

  //fill faces
  std::vector<vertex_descriptor> face;
  face.resize(3);
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
    put(vpmap, top, P(base_center.x(), height+base_center.y(),base_center.z()));
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
 * If `center` is `(0, 0, 0)`, then the first point of the base is `(radius, 0, 0)`
 *
 * \param nb_vertices the number of vertices in the base. It must be greater than or equal to 3.
 * \param g the graph in which the pyramid will be created
 * \param base_center the center of the circle in which the base is inscribed.
 * \param height the distance between the base and the apex.
 * \param radius the radius of the circle in which the base is inscribed.
 * \param is_closed determines if the base must be created or not. If `is_closed` is `true`, `center` is a vertex.
 *
 * \returns the halfedge that has the target vertex associated with the apex point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_pyramid(typename boost::graph_traits<Graph>::vertices_size_type nb_vertices,
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
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
    vertices[i] = add_vertex(g);

  vertex_descriptor apex = add_vertex(g);

  //fill vertices
  put(vpmap, apex,
      P(base_center.x(),
        base_center.y() + height,
        base_center.z()));

  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
  {
    put(vpmap, vertices[i],
        P(0.5*diameter*cos(i*precision*to_rad)+base_center.x(),
          base_center.y(),
          -0.5*diameter*sin(i*precision*to_rad)+base_center.z()));
  }

  //fill faces
  std::vector<vertex_descriptor> face;
  face.resize(3);
  for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
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
    for(typename boost::graph_traits<Graph>::vertices_size_type i=0; i<nb_vertices; ++i)
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
 *
 * \brief Creates an icosahedron, outward oriented, centered in `center` and adds it to the graph `g`.
 *
 * \param g the graph in which the icosahedron will be created.
 * \param center the center of the sphere in which the icosahedron is inscribed.
 * \param radius the radius of the sphere in which the icosahedron is inscribed.
 *
 * \returns the halfedge that has the target vertex associated with the first point in the first face.
 */
template<class Graph, class P>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_icosahedron(Graph& g,
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

  typename CGAL::Kernel_traits<P>::Kernel::FT t = radius * (1.0 + CGAL::approximate_sqrt(5.0)) / 2.0;

  put(vpmap, v_vertices[0], P(-radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[1], P( radius + center.x(),  t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[2], P(-radius + center.x(), -t + center.y(), 0.0 + center.z()));
  put(vpmap, v_vertices[3], P( radius + center.x(), -t + center.y(), 0.0 + center.z()));

  put(vpmap, v_vertices[4], P( 0.0 + center.x(), -radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[5], P( 0.0 + center.x(),  radius + center.y(),  t + center.z()));
  put(vpmap, v_vertices[6], P( 0.0 + center.x(), -radius + center.y(), -t + center.z()));
  put(vpmap, v_vertices[7], P( 0.0 + center.x(),  radius + center.y(), -t + center.z()));

  put(vpmap, v_vertices[8], P(  t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[9], P(  t + center.x(), 0.0 + center.y(),  radius + center.z()));
  put(vpmap, v_vertices[10], P(-t + center.x(), 0.0 + center.y(), -radius + center.z()));
  put(vpmap, v_vertices[11], P(-t + center.x(), 0.0 + center.y(),  radius + center.z()));

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

  face[1] = v_vertices[1]; face[0] = v_vertices[9]; face[2] = v_vertices[5];
  Euler::add_face(face, g);
  face[1] = v_vertices[5]; face[0] = v_vertices[4]; face[2] = v_vertices[11];
  Euler::add_face(face, g);
  face[1] = v_vertices[11]; face[0] = v_vertices[2]; face[2] = v_vertices[10];
  Euler::add_face(face, g);
  face[1] = v_vertices[10]; face[0] = v_vertices[6]; face[2] = v_vertices[7];
  Euler::add_face(face, g);
  face[1] = v_vertices[7]; face[0] = v_vertices[8]; face[2] = v_vertices[1];
  Euler::add_face(face, g);

  face[1] = v_vertices[3]; face[0] = v_vertices[4]; face[2] = v_vertices[9];
  Euler::add_face(face, g);
  face[1] = v_vertices[3]; face[0] = v_vertices[2]; face[2] = v_vertices[4];
  Euler::add_face(face, g);
  face[1] = v_vertices[3]; face[0] = v_vertices[6]; face[2] = v_vertices[2];
  Euler::add_face(face, g);
  face[1] = v_vertices[3]; face[0] = v_vertices[8]; face[2] = v_vertices[6];
  Euler::add_face(face, g);
  face[1] = v_vertices[3]; face[0] = v_vertices[9]; face[2] = v_vertices[8];
  Euler::add_face(face, g);

  face[1] = v_vertices[4]; face[0] = v_vertices[5]; face[2] = v_vertices[9];
  Euler::add_face(face, g);
  face[1] = v_vertices[2]; face[0] = v_vertices[11]; face[2] = v_vertices[4];
  Euler::add_face(face, g);
  face[1] = v_vertices[6]; face[0] = v_vertices[10]; face[2] = v_vertices[2];
  Euler::add_face(face, g);
  face[1] = v_vertices[8]; face[0] = v_vertices[7]; face[2] = v_vertices[6];
  Euler::add_face(face, g);
  face[1] = v_vertices[9]; face[0] = v_vertices[1]; face[2] = v_vertices[8];
  Euler::add_face(face, g);

  return halfedge(v_vertices[1], v_vertices[0], g).first;
}

/*!
 * \ingroup PkgBGLHelperFct
 *
 * \brief Creates a row major ordered grid with `i` cells along the width and `j` cells
 * along the height and adds it to the graph `g`.
 * An internal property map for `CGAL::vertex_point_t` must be available in `Graph`.
 *
 * \param i the number of cells along the width.
 * \param j the number of cells along the height.
 * \param g the graph in which the grid will be created.
 * \param calculator the functor that will assign coordinates to the grid vertices.
 * \param triangulated decides if a cell is composed of one quad or two triangles.
 * If `triangulated` is `true`, the diagonal of each cell is oriented from (0,0) to (1,1)
 * in the cell coordinates.
 *
 *\tparam CoordinateFunctor a function object providing:
 * `%Point_3 operator()(size_type I, size_type J)`, with `%Point_3` being the value_type
 * of the internal property_map for `CGAL::vertex_point_t` and outputs an object of type
 * `boost::property_traits<boost::property_map<Graph,CGAL::vertex_point_t>::%type>::%value_type`.
 *  It will be called with arguments (`w`, `h`), with `w` in [0..`i`] and `h` in [0..`j`].<br>
 * %Default: a point with positive integer coordinates (`w`, `h`, 0), with `w` in [0..`i`] and `h` in [0..`j`]
 *
 * \returns the non-border non-diagonal halfedge that has the target vertex associated with the first point of the grid (default is (0,0,0) ).
 */
template<class Graph, class CoordinateFunctor>
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
  //create the vertices
  std::vector<vertex_descriptor> v_vertices;
  v_vertices.resize(static_cast<std::size_t>(w*h));
  for(std::size_t k = 0; k < v_vertices.size(); ++k)
    v_vertices[k] = add_vertex(g);
  //assign the coordinates
  for(typename boost::graph_traits<Graph>::vertices_size_type a=0; a<w; ++a)
    for(typename boost::graph_traits<Graph>::vertices_size_type b=0; b<h; ++b)
      put(vpmap, v_vertices[a+w*b], calculator(a,b));

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

template<class Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
make_grid(typename boost::graph_traits<Graph>::vertices_size_type w,
          typename boost::graph_traits<Graph>::vertices_size_type h,
          Graph& g,
          bool triangulated = false)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type Size_type;
  typedef typename boost::property_traits<typename boost::property_map<Graph, vertex_point_t>::type>::value_type Point;

  return make_grid(w, h, g, internal::Default_grid_maker<Size_type, Point>(), triangulated);
}

} // namespace CGAL

// Here at the bottom because helpers.h must include generators (for backward compatibility reasons),
// and Euler_operations.h needs helpers.h
#include <CGAL/boost/graph/Euler_operations.h>

#endif // CGAL_BOOST_GRAPH_GENERATORS_H
