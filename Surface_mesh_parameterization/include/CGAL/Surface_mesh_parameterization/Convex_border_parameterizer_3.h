// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.1/Surface_mesh_parameterization/include/CGAL/Surface_mesh_parameterization/Convex_border_parameterizer_3.h $
// $Id: Convex_border_parameterizer_3.h 9c1ad66 2022-08-24T08:35:13+01:00 Andreas Fabri
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sirio Bolaños Puchet

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_CONVEX_BORDER_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_CONVEX_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/assertions.h>

#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/circulator.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/convexity_check_2.h>


#include <cfloat>
#include <climits>
#include <vector>

#define CGAL_SMP_CBP_DEBUG_L0 0

/// \file Convex_border_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

//
// Class Convex_border_parameterizer_3
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a convex polygon. If the input point set is
/// not strongly convex, the convex hull is computed and used instead.
///
/// `Convex_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated. It implements most of the algorithm. Subclasses only
/// have to implement the function `compute_edge_length()` to compute a segment's
/// length.
///
/// The user can provide vertices on the border of the mesh, which will be
/// mapped to the vertices of the convex polygon. If the input point set is
/// not strongly convex, the number of vertices provided must match the
/// (possibly unknown) number of vertices in the convex hull.
///
/// \attention The convex border parameterizer may create degenerate faces in the parameterization:
///            if an input border vertex has valence `1` and if it is mapped to the same edge of the
///            convex polygon as its two adjacent (border) vertices, for example.
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
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_uniform_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_arc_length_parameterizer_3<TriangleMesh>`
///
template <class TriangleMesh_>
class Convex_border_parameterizer_3
{
// Public types
public:
  /// Triangle mesh type
  typedef TriangleMesh_                                   Triangle_mesh;

  typedef TriangleMesh_                                   TriangleMesh;

  /// Mesh vertex type
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;

// Protected types
protected:
  typedef Halfedge_around_face_iterator<Triangle_mesh>              halfedge_around_face_iterator;

  // Traits subtypes:
  typedef typename internal::Kernel_traits<Triangle_mesh>::PPM      PPM;
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel   Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Vector_3                                 Vector_3;
  typedef typename Kernel::Segment_2                                Segment_2;

  typedef typename std::vector<NT>                                  Offset_map;

// Private members
private:
  int N;
  double perimeter;
  std::vector<Segment_2> seg;
  std::vector<vertex_descriptor> vx;

// Protected operations
protected:
  virtual double compute_edge_length(const Triangle_mesh& mesh,
                                     vertex_descriptor source,
                                     vertex_descriptor target) const = 0;

// Private operations
private:
  // Compute the total length of the border.
  double compute_border_length(const Triangle_mesh& mesh,
                               halfedge_descriptor bhd) const
  {
    double len = 0.0;
    for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
      len += compute_edge_length(mesh, source(hd, mesh), target(hd, mesh));
    }
    return len;
  }

  // Utility method for parameterize().
  // Compute the mesh iterator whose offset is closest to 'value'.
  halfedge_around_face_iterator closest_iterator(const Triangle_mesh& mesh,
                                                 halfedge_descriptor bhd,
                                                 Offset_map& offset,
                                                 double value) const
  {
    halfedge_around_face_iterator b, e, best;
    double min = DBL_MAX; // distance for 'best'

    std::size_t id = 0, min_id = (std::numeric_limits<std::size_t>::max)();
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

  // Set the corners by splitting the border of the mesh in N
  // segments with approximately same relative lengths as the
  // edges of the convex polygon.
  std::vector<vertex_descriptor> compute_isodistant_corners(const Triangle_mesh& mesh,
                                                            halfedge_descriptor bhd,
                                                            Offset_map& offset) const
  {
    // map to [0,perimeter[
    double len = 0.0; // current position on convex hull in [0, total_len[
    double total_len = compute_border_length(mesh, bhd);

    halfedge_around_face_iterator b, e;
    boost::tie(b,e) =  halfedges_around_face(bhd, mesh);
    for(halfedge_around_face_iterator it = b; it!= e; ++it) {
      vertex_descriptor vs = source(*it, mesh);
      vertex_descriptor vt = target(*it, mesh);

      offset.push_back(perimeter * len / total_len);
                              // current position on convex polygon in [0,perimeter[

      len += compute_edge_length(mesh, vs, vt);
    }

    std::vector<vertex_descriptor> vx_new;
    // First convex polygon corner is mapped to first vertex.
    // Then find closest points for all other corners.
    halfedge_around_face_iterator it = b;
    vertex_descriptor vd = source(*it, mesh);
    double runlen = sqrt(seg[0].squared_length());
    vx_new.push_back(vd); // add corner
    offset[0] = 0; // snap the vertex to the first corner

    for(int i = 1; i < N; ++i) {
        it = closest_iterator(mesh, bhd, offset, runlen);
        vd = source(*it, mesh);
        runlen += sqrt(seg[i].squared_length());
        vx_new.push_back(vd); // add corner
    }

    offset.clear();

    return vx_new;
  }

  // Compute the offset values for all the vertices of the border of
  // the mesh. The vertices between two given vertices vi and vj are
  // sent to the same edge of the convex polygon.
  template<typename VertexParameterizedMap>
  halfedge_descriptor compute_offsets(const Triangle_mesh& mesh,
                                      halfedge_descriptor bhd,
                                      VertexParameterizedMap vpmap,
                                      Offset_map& offset)
  {
    CGAL_assertion(offset.empty());

    for(int i = 0; i < N; ++i)
        put(vpmap, vx[i], true);

    // Move till the border halfedge has a given vertex as source
    halfedge_descriptor start_hd = bhd;
    while(!get(vpmap, source(start_hd, mesh)))
      start_hd = next(start_hd, mesh);

    // Initialize the offset for each side
    double len = 0.0;
    std::size_t index_of_previous_corner = 0, current_index = 0;
    double corner_offset = 0.0;
    for(halfedge_descriptor hd : halfedges_around_face(start_hd, mesh)) {
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
        // If both extremities of a segment are corners, offsets are already correct
        if(!get(vpmap, vs)) {
          CGAL_assertion(len != 0.0);
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

  /// assigns to the vertices of the border of the mesh a 2D position
  /// (i.e.\ a (u,v) pair) on the border's shape. Mark them as <i>parameterized</i>.
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
  /// \param uvmap an instantiation of the class `VertexUVmap`.
  /// \param vpmap an instantiation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (vimap must be initialized).
  ///
  template<typename VertexUVMap,
           typename VertexIndexMap,
           typename VertexParameterizedMap>
  Error_code parameterize(const Triangle_mesh& mesh,
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
    if(size_of_border < 3) { // minimal shape allowed is a triangle
      return ERROR_BORDER_TOO_SHORT;
    }

    halfedge_descriptor start_hd; // border halfedge whose source is the first corner

    // The offset is a vector of double between 0 and N that gives the position
    // of the vertex through its (edge-wise normalized) distance to v0 on the "unrolled" border
    Offset_map offset;

    if(vx.empty()) {
      vx = compute_isodistant_corners(mesh, bhd, offset);
    } else { // make sure that the provided vertices all belong to the border defined by 'bhd'
        unsigned int v_counter = 0;
        if(!vx.empty()) {
          for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
            vertex_descriptor vd = source(hd, mesh);
            for(int i = 0; i < N; ++i) { // FIXME: N*M lookup
              if(vd == vx[i]) {
                  v_counter++;
                  break;
              }
            }
          }

          if(v_counter != N) {
            std::cerr << "Error: Fixed vertices must belong to the same border";
            std::cerr << " (defined by 'bhd')." << std::endl;
            return ERROR_NON_CONVEX_BORDER;
          }
        }
    }
    // Here we obtain offsets "normalized" to each edge, i.e. every interval [n, n+1[
    // with n integer in 0 .. N-1, maps to one edge in the convex polygon
    start_hd = compute_offsets(mesh, bhd, vpmap, offset);

    if(CGAL_SMP_CBP_DEBUG_L0)
      for(int i = 0; i < N; ++i)
        std::cout << vx[i] << std::endl;

    unsigned int corners_encountered = 0;
    std::size_t counter = 0;

    for(halfedge_descriptor hd : halfedges_around_face(start_hd, mesh)) {
      vertex_descriptor vd = source(hd, mesh);
      Point_2 uv;
      Segment_2 edge;
      NT lambda;
      CGAL_assertion(counter < offset.size());

      edge = seg[corners_encountered];
      lambda = offset[counter++] - (corners_encountered * 1.0); // in [0,1[
      CGAL_assertion(lambda < 1);
      uv = edge.source() + lambda * edge.to_vector();

      put(uvmap, vd, uv);
      put(vpmap, vd, true); // Mark vertex as parameterized

      if(get(vpmap, target(hd, mesh))) // is a corner
        ++corners_encountered;
    }

    return OK;
  }

  /// indicates if the border's shape is convex.
  bool is_border_convex() const { return true; }


public:
  virtual ~Convex_border_parameterizer_3() { }

  /// Constructor with point set
  Convex_border_parameterizer_3(std::vector<std::pair<double,double>> pvec)
  {
    std::vector<Point_2> points, result;
    for(std::pair<double,double> p : pvec)
      points.push_back(Point_2(p.first,p.second));
    if(is_ccw_strongly_convex_2(points.begin(), points.end())) {
      if(CGAL_SMP_CBP_DEBUG_L0)
        std::cerr << "Info: Point set is already strongly convex" << std::endl;
      result = points;
    } else
      convex_hull_2(points.begin(), points.end(), std::back_inserter(result));
    N = result.size();
    perimeter = 0.0;
    for(int i = 0; i < N; ++i) {
      Point_2 src = result[i];
      Point_2 tgt = result[(i + 1) % N];
      Segment_2 s(src,tgt);
      perimeter += sqrt(s.squared_length());
      seg.push_back(s);
    }
  }

  /// Constructor with point set and user-defined corners
  Convex_border_parameterizer_3(std::vector<std::pair<double,double>> pvec,
                                std::vector<vertex_descriptor> vx)
    :
      vx(vx)
  {
    std::vector<Point_2> points, result;
    for(std::pair<double,double> p : pvec)
      points.push_back(Point_2(p.first,p.second));
    if(is_ccw_strongly_convex_2(points.begin(), points.end())) {
      std::cerr << "Info: Point set is already strongly convex" << std::endl;
      result = points;
    } else
      convex_hull_2(points.begin(), points.end(), std::back_inserter(result));
    N = result.size();
    for(int i = 0; i < N; ++i) {
      Point_2 src = result[i];
      Point_2 tgt = result[(i + 1) % N];
      seg.push_back(Segment_2(src,tgt));
    }
  }
};

//
// Class Convex_border_uniform_parameterizer_3
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a convex polygon
/// in a uniform manner: points are equally spaced.
///
/// Convex_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \attention The convex border parameterizer may create degenerate faces in the parameterization:
///            if an input border vertex has valence `1` and if it is mapped to the same edge of the
///            convex polygon as its two adjacent (border) vertices, for example.
///
/// \cgalModels `Parameterizer_3`
///
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_arc_length_parameterizer_3<TriangleMesh>`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
template<class TriangleMesh_>
class Convex_border_uniform_parameterizer_3
  : public Convex_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  // We have to repeat the types exported by superclass
  typedef TriangleMesh_                                                  TriangleMesh;
  typedef TriangleMesh_                                                  Triangle_mesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;

// Private types
private:
  typedef Convex_border_parameterizer_3<TriangleMesh_>                   Base;
  typedef typename Base::NT                                              NT;

public:
  virtual ~Convex_border_uniform_parameterizer_3() { }

  Convex_border_uniform_parameterizer_3(std::vector<std::pair<double,double>> pvec)
    :
      Base(pvec)
  { }

  Convex_border_uniform_parameterizer_3(std::vector<std::pair<double,double>> pvec,
                                        std::vector<vertex_descriptor> vx)
    :
      Base(pvec, vx)
  { }

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
};

//
// Class Convex_border_arc_length_parameterizer_3
//

/// \ingroup PkgSurfaceMeshParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a convex polygon,
/// with an arc-length parameterization: `(u,v)` values are
/// proportional to the length of border edges.
///
/// Convex_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \attention The convex border parameterizer may create degenerate faces in the parameterization:
///            if an input border vertex has valence `1` and if it is mapped to the same edge of the
///            convex polygon as its two adjacent (border) vertices, for example.
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \cgalModels `Parameterizer_3`
///
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_parameterizer_3<TriangleMesh>`
/// \sa `CGAL::Surface_mesh_parameterization::Convex_border_uniform_parameterizer_3<TriangleMesh>`
///
template<class TriangleMesh_>
class Convex_border_arc_length_parameterizer_3
  : public Convex_border_parameterizer_3<TriangleMesh_>
{
// Public types
public:
  typedef TriangleMesh_                                                  TriangleMesh;
  typedef TriangleMesh_                                                  Triangle_mesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;

// Private types
private:
  typedef Convex_border_parameterizer_3<Triangle_mesh>          Base;
  typedef typename Base::PPM                                    PPM;
  typedef typename Base::NT                                     NT;
  typedef typename Base::Vector_3                               Vector_3;

public:
  virtual ~Convex_border_arc_length_parameterizer_3() { }

  Convex_border_arc_length_parameterizer_3(std::vector<std::pair<double,double>> pset)
    :
      Base(pset)
  { }

  Convex_border_arc_length_parameterizer_3(std::vector<std::pair<double,double>> pset,
                                           std::vector<vertex_descriptor> vx)
    :
      Base(pset, vx)
  { }

// Protected operations
protected:
  /// Compute the length of an edge.
  virtual NT compute_edge_length(const Triangle_mesh& mesh,
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

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_CONVEX_BORDER_PARAMETERIZER_3_H
