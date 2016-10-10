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

/// \file Circular_border_parameterizer_3.h

namespace CGAL {

//
// Class Circular_border_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a circle.
/// `Circular_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated.
/// It implements most of the algorithm.
///

/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<TriangleMesh>`

template<class TriangleMesh_> ///< a model of `FaceGraph`
class Circular_border_parameterizer_3
{
// Public types
public:
  typedef TriangleMesh_ TriangleMesh;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  // Private types
private:
  typedef Parameterizer_traits_3<TriangleMesh>  Traits;
  typedef typename Traits::VPM                  VPM;
  typedef typename Traits::Point_3              Point_3;
  typedef typename Traits::Vector_3             Vector_3;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Error_code           Error_code;

// Public operations
public:
  /// Destructor of base class should be virtual.
  virtual ~Circular_border_parameterizer_3() { }

  // Default constructor, copy constructor and operator =() are fine

  /// Assign to mesh's border vertices a 2D position (i.e.\ a (u,v) pair)
  /// on border's shape. Mark them as <i>parameterized</i>.
  template <typename VertexUVmap, typename VertexParameterizedMap>
  Error_code
  parameterize_border(const TriangleMesh& mesh,
                      halfedge_descriptor bhd,
                      VertexUVmap uvmap,
                      VertexParameterizedMap vpmap)
  {
    VPM vpm = get(vertex_point, mesh);
    // TODO  Nothing to do if no border
    //if (! is_border(bhd,tmesh)){
    // return Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;
    //}

    // Compute the total border length
    double total_len = compute_border_length(mesh,bhd);
    if (total_len == 0)
      return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    const double tmp = 2*CGAL_PI/total_len;
    double len = 0.0;           // current position on circle in [0, total_len]

    Halfedge_around_face_circulator<TriangleMesh> circ(bhd,mesh), done(circ);
    do {
      halfedge_descriptor hd = *circ;
      --circ;
      vertex_descriptor vd = target(opposite(hd,mesh),mesh);

      double angle = len*tmp; // current position on the circle in radians

      // map vertex on unit circle
      Point_2 uv(0.5+0.5*std::cos(-angle),0.5+0.5*std::sin(-angle));
      put(uvmap, vd, uv);
      put(vpmap, vd, true);

      len += CGAL::sqrt(squared_distance(get(vpm, target(hd,mesh)), get(vpm,vd)));
    } while(circ != done);

    return Parameterizer_traits_3<TriangleMesh>::OK;
  }


/// Indicate if border's shape is convex.
  bool is_border_convex () { return true; }


private:
  /// Compute the total length of the border
  double compute_border_length(const TriangleMesh& tmesh, halfedge_descriptor bhd)
  {
    VPM vpm = get(CGAL::vertex_point,tmesh);
    double len = 0.0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd,tmesh)){
      len += CGAL::sqrt(squared_distance(get(vpm, source(hd,tmesh)),
                                         get(vpm, target(hd,tmesh))));
    }
    return len;
  }
};

//
// Class Circular_border_arc_length_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle,
/// with an arc-length parameterization: (u,v) values are
/// proportional to the length of border edges.
/// `Circular_border_parameterizer_3` implements most of the border parameterization
/// algorithm.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<TriangleMesh>`

template<class TriangleMesh_> ///< 3D surface
  class Circular_border_arc_length_parameterizer_3
    : public Circular_border_parameterizer_3<TriangleMesh_>
{
  // Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef TriangleMesh_          TriangleMesh;
  /// @endcond

  // Private types
private:
  typedef typename Parameterizer_traits_3<TriangleMesh>::Point_2 Point_2;
  typedef typename Parameterizer_traits_3<TriangleMesh>::Vector_3 Vector_3;


  // Public operations
public:
   typedef typename Circular_border_parameterizer_3<TriangleMesh_>::vertex_descriptor vertex_descriptor;
   typedef typename Circular_border_parameterizer_3<TriangleMesh_>::halfedge_descriptor halfedge_descriptor;
// Private types
  // Default constructor, copy constructor and operator =() are fine

  // Protected operations
protected:
  /// Compute the length of an edge.
  virtual double compute_edge_length(const TriangleMesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target)
  {
    typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(vertex_point, mesh);
    /// Arc-length border parameterization: (u,v) values are
    /// proportional to the length of border edges.
    Vector_3 v = get(ppmap, target) - get(ppmap,source);
    return std::sqrt(v*v);
  }
};

} // namespace CGAL

#endif // CGAL_CIRCULAR_BORDER_PARAMETERIZER_3_H
