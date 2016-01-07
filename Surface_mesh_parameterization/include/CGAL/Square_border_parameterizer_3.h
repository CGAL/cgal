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


#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Parameterizer_traits_3.h>

#include <cfloat>
#include <climits>
#include <vector>

/// \file Square_border_parameterizer_3.h

namespace CGAL {


//
// Class Square_border_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a square.
/// `Square_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated.
///
/// It implements most of the algorithm. Subclasses just
/// have to implement `compute_edge_length(`) to compute a segment's length.
///
/// Implementation note:
/// To simplify the implementation, `BorderParameterizer_3` models know only the
/// `ParameterizationMesh_3` class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `BorderParameterizer_3`
///

template<class ParameterizationMesh_3>      //< 3D surface
class Square_border_parameterizer_3
{
// Public types
public:
    /// Export ParameterizationMesh_3 template parameter.
    typedef ParameterizationMesh_3  TriangleMesh;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef Halfedge_around_face_iterator<TriangleMesh> halfedge_around_face_iterator;

// Private types
private:
    // Mesh_TriangleMesh_3 subtypes:
  typedef Parameterizer_traits_3<TriangleMesh> Traits;
  typedef typename Traits::VPM      VPM;
  typedef typename Traits::Point_3  Point_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::Point_2  Point_2;

    typedef typename std::vector<double>    Offset_map;

// Public operations
public:
    /// Destructor of base class should be virtual.
    virtual ~Square_border_parameterizer_3() {}

    // Default constructor, copy constructor and operator =() are fine

    /// Assign to mesh's border vertices a 2D position (i.e.\ a (u,v) pair)
    /// on border's shape. Mark them as <i>parameterized</i>.
    typename Parameterizer_traits_3<TriangleMesh>::Error_code
                                        parameterize_border(TriangleMesh& mesh);

    /// Indicate if border's shape is convex.
    bool  is_border_convex () { return true; }

 
// Private operations
private:
    /// Compute the total length of the border.
  double compute_border_length(const TriangleMesh& tmesh, halfedge_descriptor bhd);

    /// Get mesh iterator whose offset is closest to 'value'.
    halfedge_around_face_iterator closest_iterator(TriangleMesh& mesh,
                                                   const Offset_map& offsets,
                                                   double value);
};


// Compute the total length of the border.
template<class TriangleMesh>
inline
double Square_border_parameterizer_3<TriangleMesh>::compute_border_length(
                                                        const TriangleMesh& tmesh, halfedge_descriptor bhd)
{
  VPM vpm = get(CGAL::vertex_point,tmesh);
  double len = 0.0;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd,tmesh)){
    len += CGAL::sqrt(squared_distance(get(vpm, source(hd,tmesh)), get(vpm, target(hd,tmesh))));
  }
  return len;
}

// Assign to mesh's border vertices a 2D position (i.e. a (u,v) pair)
// on border's shape. Mark them as "parameterized".
template<class TriangleMesh>
inline
typename Parameterizer_traits_3<TriangleMesh>::Error_code
Square_border_parameterizer_3<TriangleMesh>::parameterize_border(TriangleMesh& tmesh,
                      halfedge_descriptor bhd,
                      VertexUVmap uvmap,
                      VertexParameterizedMap vpmap)
{
    VPM vpm = get(vertex_point, tmesh);
#ifdef DEBUG_TRACE
    std::cerr << "  map on a square" << std::endl;
#endif
    
    // Nothing to do if no border
    // if (mesh.main_border().empty())
    //    return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(tmesh,bhd);
    if (total_len == 0)
        return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;

    // map to [0,4[
    double len = 0.0;           // current position on square in [0, total_len[
    Offset_map offset;          // vertex index -> offset map
    offset.resize(mesh.count_mesh_vertices());

    halfedge_around_face_iterator b,e;
    BOOST_FOREACH(halfedge_descriptor hd,  halfedges_around_face(bhd,tmesh)){
    
      vertex_descriptor vs = source(hd,tmesh);
      vertex_descriptor vt = target(hd,tmesh);

        offset[mesh.get_vertex_index(vs)] = 4.0f*len/total_len;
                                // current position on square in [0,4[

        len += compute_edge_length(mesh, vs, vt);
    }

    // First square corner is mapped to first vertex.
    // Then find closest points for three other corners.
    halfedge_around_face_iterator it0 = b;
    halfedge_around_face_iterator it1 = closest_iterator(mesh, offset, 1.0);
    halfedge_around_face_iterator it2 = closest_iterator(mesh, offset, 2.0);
    halfedge_around_face_iterator it3 = closest_iterator(mesh, offset, 3.0);
    //
    // We may get into trouble if the border is too short
    if (it0 == it1 || it1 == it2 || it2 == it3 || it3 == it0)
        return Parameterizer_traits_3<TriangleMesh>::ERROR_BORDER_TOO_SHORT;
    //
    // Snap these vertices to corners
    offset[mesh.get_vertex_index(source(*it0,tmesh))] = 0.0;
    offset[mesh.get_vertex_index(source(*it1,tmesh))] = 1.0;
    offset[mesh.get_vertex_index(source(*it2,tmesh))] = 2.0;
    offset[mesh.get_vertex_index(source(*it3,tmesh))] = 3.0;

    // Set vertices along square's sides and mark them as "parameterized"
    for(halfedge_around_face_iterator it = it0; it != it1; it++) // 1st side
    {
      Point_2 uv(offset[mesh.get_vertex_index(source(*it,tmesh))], 0.0);
        mesh.set_vertex_uv(source(*it,tmesh), uv);
        mesh.set_vertex_parameterized(source(*it,tmesh), true);
    }
    for(halfedge_around_face_iterator it = it1; it != it2; it++) // 2nd side
    {
        Point_2 uv(1.0, offset[mesh.get_vertex_index(source(*it,tmesh))]-1);
        mesh.set_vertex_uv(source(*it,tmesh), uv);
        mesh.set_vertex_parameterized(source(*it,tmesh), true);
    }
    for(halfedge_around_face_iterator it = it2; it != it3; it++) // 3rd side
    {
        Point_2 uv(3-offset[mesh.get_vertex_index(source(*it,tmesh))], 1.0);
        mesh.set_vertex_uv(source(*it,tmesh), uv);
        mesh.set_vertex_parameterized(source(*it,tmesh), true);
    }
    for(halfedge_around_face_iterator it = it3; it != e; it++) // 4th side
    {
        Point_2 uv(0.0, 4-offset[mesh.get_vertex_index(source(*it,tmesh))]);
        mesh.set_vertex_uv(source(*it,tmesh), uv);
        mesh.set_vertex_parameterized(source(*it,tmesh), true);
    }

    return Parameterizer_traits_3<TriangleMesh>::OK;
}

// Utility method for parameterize_border().
// Compute mesh iterator whose offset is closest to 'value'.
template<class TriangleMesh>
inline
typename Square_border_parameterizer_3<TriangleMesh>::halfedge_around_face_iterator
Square_border_parameterizer_3<TriangleMesh>::closest_iterator(TriangleMesh& mesh,
                                                       const Offset_map& offset,
                                                       double value)
{
    const TriangleMesh& tmesh = mesh.get_adapted_mesh();
    halfedge_around_face_iterator b, e, best;
    double min = DBL_MAX;           // distance for 'best'

    for(boost::tie(b,e) = halfedges_around_face(mesh.main_border(), tmesh);
        b!=e;
        ++b)
    {
      double d = CGAL::abs(offset[mesh.get_vertex_index(source(*b,tmesh))] - value);
        if (d < min)
        {
            best = b;
            min = d;
        }
    }

    return best;
 }


//
// Class Square_border_uniform_parameterizer_3
//


/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
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

template<class ParameterizationMesh_3>      //< 3D surface
class Square_border_uniform_parameterizer_3
    : public Square_border_parameterizer_3<ParameterizationMesh_3>
{
// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef ParameterizationMesh_3          TriangleMesh; 
    typedef typename TriangleMesh::Polyhedron TriangleMesh;

    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

    /// @endcond

// Private types
private:
    // Mesh_TriangleMesh_3 subtypes:
    typedef typename TriangleMesh::Point_2       Point_2;
    typedef typename TriangleMesh::Vector_3      Vector_3;


// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute the length of an edge.
    virtual double compute_edge_length(const TriangleMesh& /* mesh */,
                                       vertex_descriptor /* source */,
                                       vertex_descriptor /* target */)
    {
        /// Uniform border parameterization: points are equally spaced.
        return 1;
    }
};


//
// Class Square_border_arc_length_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a square,
/// with an arc-length parameterization: (u,v) values are
/// proportional to the length of border edges.
///
/// Square_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only compute_edge_length() to compute a
/// segment's length.
///
/// \cgalModels `BorderParameterizer_3`
///

template<class ParameterizationMesh_3>      //< 3D surface
class Square_border_arc_length_parameterizer_3
    : public Square_border_parameterizer_3<ParameterizationMesh_3>
{
// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef ParameterizationMesh_3          TriangleMesh;
   /// @endcond

// Private types
private:
    // Mesh_TriangleMesh_3 subtypes:
    typedef typename TriangleMesh::Point_2       Point_2;
    typedef typename TriangleMesh::Vector_3      Vector_3;

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute the length of an edge.
    virtual double compute_edge_length(const TriangleMesh& mesh,
                                       vertex_descriptor source,
                                       vertex_descriptor target)
    { 
      typedef typename boost::property_map<typename TriangleMesh::Polyhedron, boost::vertex_point_t>::const_type PPmap;
      PPmap ppmap = get(vertex_point, mesh.get_adapted_mesh());
        /// Arc-length border parameterization: (u,v) values are
        /// proportional to the length of border edges.
      Vector_3 v = get(ppmap, target)
                   - get(ppmap,source);
        return std::sqrt(v*v);
    }
};


} //namespace CGAL

#endif //CGAL_SQUARE_BORDER_PARAMETERIZER_3_H
