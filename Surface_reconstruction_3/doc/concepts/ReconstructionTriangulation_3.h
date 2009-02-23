// Copyright (c) 2007  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent Saboret, Pierre Alliez


/// The ReconstructionTriangulation_3 concept defines
/// the interface of a 3D Delaunay triangulation
/// requested by the Poisson_reconstruction_function class.
/// The cell class must be a model of
/// ReconstructionCellBase_3 and the vertex class
/// must be a model of ReconstructionVertexBase_3.
///
/// It provides the interface requested by the Poisson_reconstruction_function class:
/// - Each vertex owns a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
/// or not, and has a unique index.
///
/// CAUTION: invalidate_bounds() must be called
/// after modifying the points.
///
/// @heading Has Models:
/// Reconstruction_triangulation_3<GeomTraits, TriangulationDataStructure_3>

class ReconstructionTriangulation_3 : public DelaunayTriangulation_3,
                                      public DefaultConstructible
{
// Public types
public:

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// The geometric traits class's Point_3 type is a model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.

  /// Iterator over all normals.
  typedef xxx Normal_iterator;

  /// Point type
  static const unsigned char INPUT = 0;
  static const unsigned char STEINER = 1;

  /// Iterator over input vertices.
  typedef xxx Input_vertices_iterator;

  /// Iterator over input points.
  typedef xxx Project_point<Vertex> >  Input_point_iterator;

// Public methods
public:

  /// Get first iterator over finite vertices normals.
  Normal_iterator normals_begin();
  /// Get past-the-end iterator over finite vertices normals.
  Normal_iterator normals_end();

  /// Get first iterator over input vertices.
  Input_vertices_iterator input_vertices_begin() const;
  /// Get past-the-end iterator over input vertices.
  Input_vertices_iterator input_vertices_end() const;

  /// Get first iterator over input points.
  Input_point_iterator input_points_begin() const;
  /// Get past-the-end iterator over input points.
  Input_point_iterator input_points_end() const;
  /// Get the bounding box of all points.
  Iso_cuboid bounding_box() const;


  /// Get the bounding box of input points.
  Iso_cuboid input_points_bounding_box() const;
  /// Get the bounding sphere of all points.
  Sphere bounding_sphere() const;


  /// Get the bounding sphere of input points.
  Sphere input_points_bounding_sphere() const;

  /// Get the barycenter of all points.
  Point barycenter() const;

  /// Get the standard deviation of the distance to barycenter (for all points).
  FT diameter_standard_deviation() const;

  /// Update barycenter, bounding box, bounding sphere and standard deviation.
  /// Owner is responsible to call this function after modifying the triangulation.
  void invalidate_bounds();

  /// Insert point in the triangulation.
  /// Default type is INPUT.
  Vertex_handle insert(const Point& p,
                       unsigned char type = INPUT /* INPUT or STEINER */,
                       Cell_handle start = Cell_handle());

  /// Insert points in the triangulation using a spatial sort.
  /// Default type is INPUT.
  ///
  /// Precondition: the value type of InputIterator must be 'Point'.
  ///
  /// @param first First point to add to pdt.
  /// @param beyond Past-the-end point to add to pdt.
  /// @return the number of inserted points.
  template < class InputIterator >
  int insert(InputIterator first, InputIterator beyond,
             unsigned char type = INPUT /* INPUT or STEINER */);

  /// Delaunay refinement callback:
  /// insert STEINER point in the triangulation.
  template <class CellIt>
  Vertex_handle
  insert_in_hole(const Point & p, CellIt cell_begin, CellIt cell_end,
	         Cell_handle begin, int i,
                 unsigned char type = STEINER /* INPUT or STEINER */);

  /// Index all (finite) vertices following the order of Finite_vertices_iterator.
  /// @return the number of (finite) vertices.
  unsigned int index_vertices();

  /// Index unconstraint vertices following the order of Finite_vertices_iterator.
  /// @return the number of unconstraint vertices.
  unsigned int index_unconstrained_vertices();
};

