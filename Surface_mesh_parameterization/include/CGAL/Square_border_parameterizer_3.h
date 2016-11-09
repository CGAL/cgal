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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SQUARE_BORDER_PARAMETERIZER_3_H
#define CGAL_SQUARE_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/Parameterizer_traits_3.h>

#include <CGAL/circulator.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/foreach.hpp>

#include <cfloat>
#include <climits>
#include <fstream>
#include <string>
#include <vector>

/// \file Square_border_parameterizer_3.h

namespace CGAL {

//
// Class Square_border_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a square.
/// `Square_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated.
///
/// It implements most of the algorithm. Subclasses only have to implement
/// the function `compute_edge_length()` to compute a segment's length.
///
/// A CGAL selection file can be used to provide four chosen vertices that
/// will be the four corners of the square.
///
/// Implementation note:
/// To simplify the implementation, `BorderParameterizer_3` models know only the
/// `TriangleMesh` class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
///
/// \sa `CGAL::Square_border_uniform_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Square_border_arc_length_parameterizer_3<TriangleMesh>`
///
template <class TriangleMesh_>
class Square_border_parameterizer_3
{
// Public types
public:
  typedef TriangleMesh_  TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  typedef Halfedge_around_face_iterator<TriangleMesh>  halfedge_around_face_iterator;

// Private types
private:
    // Mesh_TriangleMesh_3 subtypes:
  typedef Parameterizer_traits_3<TriangleMesh>          Traits;
  typedef typename Traits::VPM                          VPM;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Point_3                      Point_3;
  typedef typename Traits::Vector_3                     Vector_3;
  typedef typename Traits::Error_code                   Error_code;

  typedef typename std::vector<double>                  Offset_map;


// Protected operations
protected:
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target) const = 0;

// Private operations
private:
  /// Compute the total length of the border.
  double compute_border_length(const TriangleMesh& mesh,
                               halfedge_descriptor bhd) const
  {
    double len = 0.0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)) {
      len += compute_edge_length(mesh, source(hd, mesh), target(hd, mesh));
    }
    return len;
  }

  /// Utility method for parameterize().
  /// Compute the mesh iterator whose offset is closest to 'value'.
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

  /// Set the corners by splitting the border of the mesh in four
  /// approximately equal segments.
  template<typename VertexParameterizedMap>
  Error_code set_default_corners(const TriangleMesh& mesh,
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
    offset[0] = 0; // snap the vertex to the corner

    halfedge_around_face_iterator it1 = closest_iterator(mesh, bhd, offset, 1.0);
    halfedge_around_face_iterator it2 = closest_iterator(mesh, bhd, offset, 2.0);
    halfedge_around_face_iterator it3 = closest_iterator(mesh, bhd, offset, 3.0);

    // We may get into trouble if the border is too short
    if (it0 == it1 || it1 == it2 || it2 == it3 || it3 == it0)
      return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    vertex_descriptor vd0 = source(*it0, mesh);
    vertex_descriptor vd1 = source(*it1, mesh);
    vertex_descriptor vd2 = source(*it2, mesh);
    vertex_descriptor vd3 = source(*it3, mesh);

    put(vpmap, vd0, true);
    put(vpmap, vd1, true);
    put(vpmap, vd2, true);
    put(vpmap, vd3, true);

    return Traits::OK;
  }

  /// Set the corners. If a selection file is provided and has four valid vertices,
  /// use those. Otherwise, split the border in four approximately equal segments
  /// and initialize the corners from the extremities of these segments.
  template<typename VertexParameterizedMap>
  Error_code set_corners(const TriangleMesh& mesh,
                         halfedge_descriptor bhd,
                         VertexParameterizedMap vpmap,
                         const char* filename,
                         Offset_map& offset) const
  {
    assert(offset.empty());

    if(filename == NULL) // no file given in input
      return set_default_corners(mesh, bhd, vpmap, offset);

    std::string str = filename;
    if(str.substr(str.length() - 14) != ".selection.txt") {
      std::cout << "Error: fixed vertices must be given by a *.selection.txt file" << std::endl;
      return set_default_corners(mesh, bhd, vpmap, offset);
    }

    std::ifstream in(filename);
    std::string line;
    if(!std::getline(in, line)) {
      std::cout << "Error: could not read input file: " << filename << std::endl;
      return set_default_corners(mesh, bhd, vpmap, offset);
    }

    // The selection file is a list of integers, so we need to build a correspondence
    // between vertices and the integers.
    std::vector<vertex_descriptor> vds;
    vds.reserve(num_vertices(mesh));
    typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator  vertex_iterator;
    vertex_iterator vi = vertices(mesh).begin(),
                    vi_end = vertices(mesh).end();
    CGAL_For_all(vi, vi_end) {
      vds.push_back(*vi);
    }

    // Get the first line and read the fixed vertex indices
    unsigned int counter = 0;
    std::istringstream point_line(line);
    std::size_t s;
    while(point_line >> s) {
      assert(s < vds.size());
      vertex_descriptor vd = vds[s];
      assert(is_border(halfedge(vd, mesh), mesh)); // must be on the border
      put(vpmap, vd, true);
      ++counter;
    }

    if(counter != 4) {
      std::cout << "Error: exactly four vertices must form the corners" << std::endl;
      return set_default_corners(mesh, bhd, vpmap, offset);
    }

    // All clear, initialize the offset vectors with values on [0;4[ with the
    // given corners

    // Move till the border halfedge has a fixed point as source
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
              offset[i] = corner_offset + ld * offset[i];
          }
        }
        len = 0.0; // reset the length of the side
        corner_offset += 1.0; // offset for the next corner
      }

      current_index++;
    }
    return Traits::OK;
  }

// Public operations
public:
  // Default constructor, copy constructor and operator =() are fine

  /// Assign to the vertices of the border of the mesh a 2D position
  /// (i.e.\ a (u,v) pair) on the border's shape. Mark them as <i>parameterized</i>.
  template<typename VertexUVMap, typename VertexParameterizedMap>
  Error_code parameterize(const TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          VertexParameterizedMap vpmap,
                          const char* filename = "../data/square_corners.selection.txt")
  {
  #ifdef DEBUG_TRACE
    std::cerr << "Map to a square" << std::endl;
  #endif

    // Nothing to do if no border
    if (bhd == halfedge_descriptor())
      return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(mesh, bhd);
    if (total_len == 0)
        return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    Offset_map offset;
    set_corners(mesh, bhd, vpmap, filename, offset);

    halfedge_descriptor start_hd = bhd;
    while(!get(vpmap, source(start_hd, mesh)))
      start_hd = next(start_hd, mesh);

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

    return Traits::OK;
  }

  /// Indicate if the border's shape is convex.
  bool is_border_convex() const { return true; }
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
/// algorithm. This class implements only compute_edge_length() to compute a
/// segment's length.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Square_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Square_border_arc_length_parameterizer_3<TriangleMesh>`
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

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual double compute_edge_length(const TriangleMesh& /* mesh */,
                                     vertex_descriptor /* source */,
                                     vertex_descriptor /* target */) const
  {
    /// Uniform border parameterization: points are equally spaced.
    return 1;
  }
};

//
// Class Square_border_arc_length_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a square,
/// with an arc-length parameterization: (u,v) values are
/// proportional to the length of border edges.
///
/// Square_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only compute_edge_length() to compute a
/// segment's length.
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Square_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Square_border_uniform_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Square_border_arc_length_parameterizer_3
  : public Square_border_parameterizer_3<TriangleMesh_>
{
public:
  /// @cond SKIP_IN_MANUAL
  typedef TriangleMesh_             TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  /// @endcond

private:
  typedef Parameterizer_traits_3<TriangleMesh>          Traits;
  typedef typename Traits::Vector_3                     Vector_3;
  typedef typename Traits::VPM                          VPM;

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target) const
  {
    VPM ppmap = get(vertex_point, mesh);

    /// Arc-length border parameterization: (u,v) values are proportional
    /// to the length of border edges.
    Vector_3 v = get(ppmap, target) - get(ppmap, source);
    return std::sqrt(v * v);
  }
};

} // namespace CGAL

#endif // CGAL_SQUARE_BORDER_PARAMETERIZER_3_H
