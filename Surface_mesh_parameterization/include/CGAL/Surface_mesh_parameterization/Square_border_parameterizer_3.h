// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_SQUARE_BORDER_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_SQUARE_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/circulator.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/foreach.hpp>

#include <cfloat>
#include <climits>
#include <vector>

/// \file Square_border_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

//
// Class Square_border_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a square.
///
/// `Square_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated. It implements most of the algorithm. Subclasses only
/// have to implement the function `compute_edge_length()` to compute a segment's
/// length.
///
/// The user can provide four vertices on the border of the mesh, which will be
/// mapped to the four corners of the square.
///
/// Implementation note:
/// To simplify the implementation, the border parameterizer knows only the
/// `TriangleMesh` class and does not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_uniform_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_arc_length_parameterizer_3<TriangleMesh>`
///
template <class TriangleMesh_>
class Square_border_parameterizer_3
{
// Public types
public:
  typedef TriangleMesh_                                   TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  typedef Halfedge_around_face_iterator<TriangleMesh>  halfedge_around_face_iterator;

// Protected types
protected:
  // Traits subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::PPM       PPM;
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel    Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Vector_3                                 Vector_3;

  typedef typename std::vector<NT>                                  Offset_map;

// Private members
private:
  vertex_descriptor v0, v1, v2, v3;
  bool vertices_given;

// Protected operations
protected:
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target) const = 0;

// Private operations
private:
  // Compute the total length of the border.
  double compute_border_length(const TriangleMesh& mesh,
                               halfedge_descriptor bhd) const
  {
    double len = 0.0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)) {
      len += compute_edge_length(mesh, source(hd, mesh), target(hd, mesh));
    }
    return len;
  }

  // Utility method for parameterize().
  // Compute the mesh iterator whose offset is closest to 'value'.
  halfedge_around_face_iterator closest_iterator(const TriangleMesh& mesh,
                                                 halfedge_descriptor bhd,
                                                 Offset_map& offset,
                                                 double value) const
  {
    halfedge_around_face_iterator b, e, best;
    double min = DBL_MAX; // distance for 'best'

    std::size_t id = 0, min_id = -1;
    for(boost::tie(b,e) = halfedges_around_face(bhd, mesh); b!=e; ++b, ++id) {
      double d = CGAL::abs(offset[id] - value);
      if(d < min) {
        best = b;
        min = d;
        min_id = id;
      }
    }

    // snap the offset value to the corner
    offset[min_id] = value;

    return best;
  }

  // Set the corners by splitting the border of the mesh in four
  // approximately equal segments.
  template<typename VertexParameterizedMap>
  halfedge_descriptor compute_offsets_without_given_vertices(const TriangleMesh& mesh,
                                                             halfedge_descriptor bhd,
                                                             VertexParameterizedMap vpmap,
                                                             Offset_map& offset) const
  {
    // map to [0,4[
    double len = 0.0; // current position on square in [0, total_len[
    double total_len = compute_border_length(mesh, bhd);

    halfedge_around_face_iterator b, e;
    boost::tie(b,e) =  halfedges_around_face(bhd, mesh);
    for(halfedge_around_face_iterator it = b; it!= e; ++it) {
      vertex_descriptor vs = source(*it, mesh);
      vertex_descriptor vt = target(*it, mesh);

      offset.push_back(4.0f * len / total_len);
                              // current position on square in [0,4[

      len += compute_edge_length(mesh, vs, vt);
    }

    // First square corner is mapped to first vertex.
    // Then find closest points for three other corners.
    halfedge_around_face_iterator it0 = b;
    offset[0] = 0; // snap the vertex to the first corner

    halfedge_around_face_iterator it1 = closest_iterator(mesh, bhd, offset, 1.0);
    halfedge_around_face_iterator it2 = closest_iterator(mesh, bhd, offset, 2.0);
    halfedge_around_face_iterator it3 = closest_iterator(mesh, bhd, offset, 3.0);

    vertex_descriptor vd0 = source(*it0, mesh);
    vertex_descriptor vd1 = source(*it1, mesh);
    vertex_descriptor vd2 = source(*it2, mesh);
    vertex_descriptor vd3 = source(*it3, mesh);

    put(vpmap, vd0, true);
    put(vpmap, vd1, true);
    put(vpmap, vd2, true);
    put(vpmap, vd3, true);

    return bhd;
  }

  // Compute the offset values for all the vertices of the border of
  // the mesh. The vertices between two given vertices vi and vj are
  // sent to the same side of the square.
  template<typename VertexParameterizedMap>
  halfedge_descriptor compute_offsets(const TriangleMesh& mesh,
                                      halfedge_descriptor bhd,
                                      VertexParameterizedMap vpmap,
                                      Offset_map& offset)
  {
    assert(offset.empty());

    put(vpmap, v0, true);
    put(vpmap, v1, true);
    put(vpmap, v2, true);
    put(vpmap, v3, true);

    // Move till the border halfedge has a given vertex as source
    halfedge_descriptor start_hd = bhd;
    while(!get(vpmap, source(start_hd, mesh)))
      start_hd = next(start_hd, mesh);

    // Initialize the offset for each side
    double len = 0.0;
    std::size_t index_of_previous_corner = 0, current_index = 0;
    double corner_offset = 0.0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(start_hd, mesh)) {
      vertex_descriptor vs = source(hd, mesh);
      vertex_descriptor vt = target(hd, mesh);

      if(get(vpmap, vs)) { // if the source is a corner
        index_of_previous_corner = current_index;
        offset.push_back(corner_offset);
      }
      else
        offset.push_back(len);

      len += compute_edge_length(mesh, vs, vt);

      // If the target is a corner vertex, we have the complete length of a side in 'len'
      // and we must "normalize" the previous entries
      if(get(vpmap, vt)) {
        // If both extremeties of a segment are corners, offsets are already correct
        if(!get(vpmap, vs)) {
          assert(len != 0.0);
          double ld = 1.0 / len;
          for(std::size_t i=index_of_previous_corner+1; i<=current_index; ++i) {
            // ld * offset[i] is in [0;1[
            offset[i] = corner_offset + ld * offset[i];
          }
        }
        len = 0.0; // reset the length of the side
        corner_offset += 1.0; // offset for the next corner
      }

      current_index++;
    }

    return start_hd;
  }

// Public operations
public:
  // Default constructor, copy constructor and operator =() are fine

  /// Assign to the vertices of the border of the mesh a 2D position
  /// (i.e.\ a (u,v) pair) on the border's shape. Mark them as <i>parameterized</i>.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (vimap must be initialized).
  ///
  template<typename VertexUVMap,
           typename VertexIndexMap,
           typename VertexParameterizedMap>
  Error_code parameterize(const TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          VertexIndexMap /* vimap */,
                          VertexParameterizedMap vpmap)
  {
    // Nothing to do if no border
    if (bhd == halfedge_descriptor())
      return ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(mesh, bhd);
    if (total_len == 0)
        return ERROR_BORDER_TOO_SHORT;

    // check the number of border edges
    std::size_t size_of_border = halfedges_around_face(bhd, mesh).size();
    if(size_of_border < 4) {
      return ERROR_BORDER_TOO_SHORT;
    }

    // make sure that the given vertices all belong to the border defined by 'bhd'
    unsigned int v_counter = 0;
    if(vertices_given) {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)) {
        vertex_descriptor vd = source(hd, mesh);
        if(vd == v0 || vd == v1 || vd == v2 || vd == v3)
          v_counter++;
      }

      if(v_counter != 4) {
        std::cerr << "Error: Fixed vertices must belong to the same border";
        std::cerr << " (defined by 'bhd')." << std::endl;
        return ERROR_NON_CONVEX_BORDER;
      }
    }

    halfedge_descriptor start_hd; // border halfedge whose source is the first corner

    // The offset is a vector of double between 0 and 4 that gives the position
    // of the vertex through its distance to v0 on the "unrolled" border
    Offset_map offset;
    if(vertices_given)
      start_hd = compute_offsets(mesh, bhd, vpmap, offset);
    else
      start_hd = compute_offsets_without_given_vertices(mesh, bhd, vpmap, offset);

    unsigned int corners_encountered = 0;
    std::size_t counter = 0;

    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(start_hd, mesh)) {
      vertex_descriptor vd = source(hd, mesh);
      Point_2 uv;
      assert(counter < offset.size());

      if(corners_encountered == 0)
        uv = Point_2(offset[counter++], 0.0);
      else if(corners_encountered == 1)
        uv = Point_2(1.0, offset[counter++] - 1.0);
      else if(corners_encountered == 2)
        uv = Point_2(3.0 - offset[counter++], 1.0);
      else if(corners_encountered == 3)
        uv = Point_2(0.0, 4.0 - offset[counter++]);

      put(uvmap, vd, uv);
      put(vpmap, vd, true);

      if(get(vpmap, target(hd, mesh)))
        ++corners_encountered;
    }

    return OK;
  }

  /// Indicate if the border's shape is convex.
  bool is_border_convex() const { return true; }

public:
  virtual ~Square_border_parameterizer_3() { }

  /// Constructor.
  Square_border_parameterizer_3()
    :
      vertices_given(false)
  { }

  /// Constructor with user-defined corners: the user provides four vertices
  /// of the border of the mesh, which will be parameterized to the corners
  /// of the square.
  ///
  /// @pre The given vertices must be on the border.
  Square_border_parameterizer_3(vertex_descriptor v0, vertex_descriptor v1,
                                vertex_descriptor v2, vertex_descriptor v3)
    :
      v0(v0), v1(v1), v2(v2), v3(v3),
      vertices_given(true)
  { }
};

//
// Class Square_border_uniform_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a square
/// in a uniform manner: points are equally spaced.
///
/// Square_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \cgalModels `Parameterizer_3`
///
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_arc_length_parameterizer_3<TriangleMesh>`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
template<class TriangleMesh_>
class Square_border_uniform_parameterizer_3
  : public Square_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef TriangleMesh_ TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  /// @endcond

// Private types
private:
  typedef Square_border_parameterizer_3<TriangleMesh_>          Base;
  typedef typename Base::NT                                     NT;

public:
  virtual ~Square_border_uniform_parameterizer_3() { }

  /// Constructor.
  Square_border_uniform_parameterizer_3() : Base() { }

  /// Constructor with user-defined corners: the user provides four vertices
  /// of the border of the mesh, which will be parameterized to the corners
  /// of the square.
  ///
  /// @pre The given vertices must be on the border.
  Square_border_uniform_parameterizer_3(vertex_descriptor v0, vertex_descriptor v1,
                                        vertex_descriptor v2, vertex_descriptor v3)
    :
      Base(v0, v1, v2, v3)
  { }

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual NT compute_edge_length(const TriangleMesh& /* mesh */,
                                 vertex_descriptor /* source */,
                                 vertex_descriptor /* target */) const
  {
    /// Uniform border parameterization: points are equally spaced.
    return 1.;
  }
};

//
// Class Square_border_arc_length_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a square,
/// with an arc-length parameterization: `(u,v)` values are
/// proportional to the length of border edges.
///
/// Square_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \cgalModels `Parameterizer_3`
///
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Square_border_uniform_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Square_border_arc_length_parameterizer_3
  : public Square_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  /// @cond SKIP_IN_MANUAL
  typedef TriangleMesh_                                         TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  /// @endcond

// Private types
private:
  typedef Square_border_parameterizer_3<TriangleMesh_>          Base;
  typedef typename Base::PPM                                    PPM;
  typedef typename Base::NT                                     NT;
  typedef typename Base::Vector_3                               Vector_3;

public:
  virtual ~Square_border_arc_length_parameterizer_3() { }

  /// Constructor.
  Square_border_arc_length_parameterizer_3() : Base() { }

  /// Constructor with user-defined corners: the user provides four vertices
  /// of the border of the mesh, which will be parameterized to the corners
  /// of the square.
  ///
  /// @pre The given vertices must be on the border.
  Square_border_arc_length_parameterizer_3(vertex_descriptor v0, vertex_descriptor v1,
                                           vertex_descriptor v2, vertex_descriptor v3)
    :
      Base(v0, v1, v2, v3)
  { }

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual NT compute_edge_length(const TriangleMesh& mesh,
                                 vertex_descriptor source,
                                 vertex_descriptor target) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    /// In this arc-length border parameterization: the `(u,v)` values are
    /// proportional to the length of the border edges.
    Vector_3 v = get(ppmap, target) - get(ppmap, source);
    return std::sqrt(v * v);
  }
};

} // Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_SQUARE_BORDER_PARAMETERIZER_3_H
