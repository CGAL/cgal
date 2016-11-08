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

#ifndef CGAL_CIRCULAR_BORDER_PARAMETERIZER_3_H
#define CGAL_CIRCULAR_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/Parameterizer_traits_3.h>

#include <CGAL/boost/graph/iterator.h>

#include <boost/foreach.hpp>

/// \file Circular_border_parameterizer_3.h

namespace CGAL {

//
// Class Circular_border_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a circle.
/// The class `Circular_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated.
///
/// It implements most of the algorithm. Subclasses only have to implement
/// the function `compute_edge_length()` to compute a segment's length.
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
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Circular_border_parameterizer_3
{
// Public types
public:
  typedef TriangleMesh_ TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  // Private types
private:
  typedef Parameterizer_traits_3<TriangleMesh>  Traits;
  typedef typename Traits::VPM                  VPM;
  typedef typename Traits::Point_3              Point_3;
  typedef typename Traits::Vector_3             Vector_3;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Error_code           Error_code;

// Protected operations
protected:
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target) = 0;

// Private operations
private:
  /// Compute the total length of the border
  double compute_border_length(const TriangleMesh& mesh, halfedge_descriptor bhd)
  {
    double len = 0.0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)) {
      len += compute_edge_length(mesh, source(hd, mesh), target(hd, mesh));
    }
    return len;
  }

// Public operations
public:
  // Default constructor, copy constructor and operator =() are fine

  /// Assign to the mesh's border vertices a 2D position (i.e.\ a (u,v) pair)
  /// on the circle. Mark them as <i>parameterized</i>.
  template <typename VertexUVmap, typename VertexParameterizedMap>
  Error_code
  parameterize_border(const TriangleMesh& mesh,
                      halfedge_descriptor bhd,
                      VertexUVmap uvmap,
                      VertexParameterizedMap vpmap)
  {
    // Nothing to do if no border
    if (bhd == halfedge_descriptor())
      return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(mesh, bhd);
    if (total_len == 0)
      return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    const double tmp = 2 * CGAL_PI / total_len;
    double len = 0.0;           // current position on circle in [0, total_len]

    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)) {
      vertex_descriptor vd = source(hd, mesh);
      double angle = len * tmp; // current position on the circle in radians

      // Map vertex on unit circle
      Point_2 uv(0.5 + 0.5*std::cos(-angle), 0.5 + 0.5*std::sin(-angle));
      put(uvmap, vd, uv);

      // Mark vertex as parameterized
      put(vpmap, vd, true);

      len += compute_edge_length(mesh, vd, target(hd, mesh));
    }

    return Parameterizer_traits_3<TriangleMesh>::OK;
  }

  /// Indicate if border's shape is convex.
  bool is_border_convex() { return true; }
};

//
// Class Circular_border_uniform_parameterizer_3
//

/// \ingroup PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle
/// in a uniform manner: points are equally spaced.
///
/// Circular_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only compute_edge_length() to compute a
/// segment's length.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
template<class TriangleMesh_>
class Circular_border_uniform_parameterizer_3
  : public Circular_border_parameterizer_3<TriangleMesh_>
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
                                     vertex_descriptor /* target */)
  {
    /// Uniform border parameterization: points are equally spaced.
    return 1.;
  }
};

//
// Class Circular_border_arc_length_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle,
/// with an arc-length parameterization: the (u,v) values are proportional
/// to the length of border edges.
/// The class `Circular_border_parameterizer_3` implements most of the border
/// parameterization algorithm.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \tparam TriangleMesh must be a model of `FaceGraph`.
///
/// \sa `CGAL::Circular_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Circular_border_arc_length_parameterizer_3
  : public Circular_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef TriangleMesh_          TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  /// @endcond

// Private types
private:
  typedef Parameterizer_traits_3<TriangleMesh>          Traits;
  typedef typename Traits::Vector_3                     Vector_3;
  typedef typename Traits::VPM                          VPM;

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target)
  {
    VPM ppmap = get(vertex_point, mesh);

    /// Arc-length border parameterization: (u,v) values are proportional
    /// to the length of border edges.
    Vector_3 v = get(ppmap, target) - get(ppmap, source);
    return std::sqrt(v * v);
  }
};

} // namespace CGAL

#endif // CGAL_CIRCULAR_BORDER_PARAMETERIZER_3_H
