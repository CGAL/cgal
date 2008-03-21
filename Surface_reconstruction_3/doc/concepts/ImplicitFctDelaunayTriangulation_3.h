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
//
// Author(s)     : Laurent Saboret, Pierre Alliez


/// The ImplicitFctDelaunayTriangulation_3 concept defines
/// the interface of a 3D Delaunay triangulation
/// requested by the Poisson_implicit_function class.
/// The cell base class must be a model of
/// ImplicitFctDelaunayTriangulationCellBase_3 and the vertex base class
/// must be a model of ImplicitFctDelaunayTriangulationVertexBase_3.
///
/// @heading Refines: DelaunayTriangulation_3,
///                   DefaultConstructible, CopyConstructible, Assignable.
///
/// @heading Has Models: 
/// Implicit_fct_delaunay_triangulation_3<GeomTraits, TriangulationDataStructure_3>

class ImplicitFctDelaunayTriangulation_3 : public DelaunayTriangulation_3
{
// Public types
public:

  // Geometric types
	typedef Geom_traits::FT FT;
	typedef Geom_traits::Vector_3 Vector;
	typedef Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
	typedef Geom_traits::Sphere_3 Sphere;
	
  /// The geometric traits class's Point_3 type is a model of PointWithNormal_3
	typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
	typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of OrientedNormal_3 concept.

  /// Iterator over normals
  typedef xxx Normal_iterator;

	/// Point type
	static const unsigned char INPUT = 0;
	static const unsigned char STEINER = 1;

// Public methods
public:

  /// Get first iterator over finite vertices normals.
  Normal_iterator normals_begin();
  /// Get past-the-end iterator over finite vertices normals.
  Normal_iterator normals_end();

  /// Get the bounding box.
	Iso_cuboid_3 bounding_box() const;

  /// Get bounding sphere.
	Sphere bounding_sphere() const;

	/// Get points barycenter.
	Point barycenter() const;

	/// Get the standard deviation of the distance to barycenter.
	FT standard_deviation() const;

  /// Update barycenter, bounding box, bounding sphere and standard deviation.
  /// Owner is responsible to call this function after modifying the triangulation.
	void invalidate_bounding_box();

  /// Insert point to the triangulation.
  Vertex_handle insert(const Point& p, 
                       unsigned char type /* INPUT or STEINER */, 
                       Cell_handle start = Cell_handle());

  /// Insert points to the triangulation using a spatial sort.
  ///
  /// Precondition: the value type of InputIterator must 'Point'.
  ///
  /// @param first First point to add to pdt.
  /// @param beyond Past-the-end point to add to pdt.
  /// @return the number of inserted points.
  template < class InputIterator >
  int insert(InputIterator first, InputIterator beyond,
             unsigned char type /* INPUT or STEINER */);

};

