// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_CIRCULAR_BORDER_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_CIRCULAR_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/number_utils.h>


#include <cmath>

/// \file Circular_border_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

//
// Class Circular_border_parameterizer_3
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
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
/// To simplify the implementation, the border parameterizer knows only the
/// `TriangleMesh` class and does not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_uniform_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
///
template< typename TriangleMesh_ >
class Circular_border_parameterizer_3
{
// Public types
public:
  /// Triangle mesh type
  typedef TriangleMesh_ Triangle_mesh;

  /// Mesh vertex type
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor   vertex_descriptor;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

// Protected types
protected:
  typedef typename internal::Kernel_traits<Triangle_mesh>::PPM       PPM;
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel    Kernel;
  typedef typename Kernel::FT                                        NT;
  typedef typename Kernel::Point_2                                   Point_2;
  typedef typename Kernel::Vector_3                                  Vector_3;

// Protected operations
protected:
  virtual NT compute_edge_length(const Triangle_mesh& mesh,
                                 vertex_descriptor source,
                                 vertex_descriptor target) const = 0;

// Private operations
private:
  // Compute the total length of the border
  NT compute_border_length(const Triangle_mesh& mesh,
                           halfedge_descriptor bhd) const
  {
    NT len = 0.0;
    for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
      len += compute_edge_length(mesh, source(hd, mesh), target(hd, mesh));
    }
    return len;
  }

// Public operations
public:
  /// assigns to the mesh's border vertices a 2D position (i.e.\ a `(u,v)` pair)
  /// on the circle. Mark them as <i>parameterized</i>.
  ///
  /// The distribution of vertices over the circle depends on the function
  /// `compute_edge_length()`.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Triangle_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
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
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(const Triangle_mesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap /* vimap */,
                          VertexParameterizedMap vpmap)
  {
    // Nothing to do if no border
    if (bhd == halfedge_descriptor())
      return ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    NT total_len = compute_border_length(mesh, bhd);
    if (total_len == 0)
      return ERROR_BORDER_TOO_SHORT;

    const NT tmp = 2 * CGAL_PI / total_len;
    NT len = 0.0; // current position on circle in [0, total_len]

    for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
      vertex_descriptor vd = source(hd, mesh);
      NT angle = len * tmp; // current position on the circle in radians

      // Map vertex on unit circle
      Point_2 uv(0.5 + 0.5*std::cos(-angle), 0.5 + 0.5*std::sin(-angle));
      put(uvmap, vd, uv);

      // Mark vertex as parameterized
      put(vpmap, vd, true);

      len += compute_edge_length(mesh, vd, target(hd, mesh));
    }

    return OK;
  }

  /// indicates if border's shape is convex.
  bool is_border_convex() const { return true; }

  virtual ~Circular_border_parameterizer_3() { }
};

//
// Class Circular_border_uniform_parameterizer_3
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle
/// in a uniform manner: points are equally spaced.
///
/// Circular_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \cgalModels `Parameterizer_3`
///
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
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
  typedef TriangleMesh_                                                  Triangle_mesh;
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;

// Private types
private:
  typedef Circular_border_parameterizer_3<Triangle_mesh>    Base;
  typedef typename Base::NT                                NT;

// Protected operations
protected:
  /// computes the length of an edge.
  virtual NT compute_edge_length(const Triangle_mesh& /* mesh */,
                                 vertex_descriptor /* source */,
                                 vertex_descriptor /* target */) const
  {
    /// Uniform border parameterization: points are equally spaced.
    return 1.;
  }

public:
  virtual ~Circular_border_uniform_parameterizer_3() { }
};

//
// Class Circular_border_arc_length_parameterizer_3
//

/// \ingroup  PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle,
/// with an arc-length parameterization: the `(u,v)` values are proportional
/// to the length of border edges.
/// The class `Circular_border_parameterizer_3` implements most of the border
/// parameterization algorithm.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Circular_border_uniform_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Circular_border_arc_length_parameterizer_3
  : public Circular_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  // We have to repeat the types exported by superclass
  typedef TriangleMesh_                                                    Triangle_mesh;
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

// Private types
private:
  typedef Circular_border_parameterizer_3<Triangle_mesh>    Base;

  typedef typename Base::PPM                                PPM;
  typedef typename Base::NT                                 NT;
  typedef typename Base::Vector_3                           Vector_3;

// Protected operations
protected:
  /// computes the length of an edge.
  virtual NT compute_edge_length(const Triangle_mesh& mesh,
                                 vertex_descriptor source,
                                 vertex_descriptor target) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    /// Arc-length border parameterization: `(u,v)` values are proportional
    /// to the length of border edges.
    Vector_3 v = get(ppmap, target) - get(ppmap, source);
    return CGAL::sqrt(v * v);
  }

public:
  virtual ~Circular_border_arc_length_parameterizer_3() { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_CIRCULAR_BORDER_PARAMETERIZER_3_H
