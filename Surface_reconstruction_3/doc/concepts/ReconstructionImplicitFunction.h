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


/// ReconstructionImplicitFunction is the concept of
/// implicit function representing a 3D surface
/// used by the Surface_reconstruction_3 package.
///
/// @heading Design Pattern:
/// A model of ReconstructionImplicitFunction is a
/// Strategy [GHJV95]: it implements a strategy of surface mesh reconstruction.
///
/// @heading Has Models: Poisson_implicit_function<PoissonDelaunayTriangulation_3>

class ReconstructionImplicitFunction : public Surface_mesher::ImplicitFunction
{
// Public types
public:

  typedef xxx Geom_traits; ///< Kernel's geometric traits

	typedef Geom_traits::FT FT;
	typedef Geom_traits::Point_3 Point;
	typedef Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
	typedef Geom_traits::Sphere_3 Sphere;

  typedef typename Triangulation::Point_with_normal Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Triangulation::Normal Normal; ///< Model of OrientedNormal_3 concept.

// Public methods
public:

  /// Create an empty implicit function.
  ReconstructionImplicitFunction();

  ///// Create an implicit function from a point set.
  ///// Precondition: the value type of InputIterator must be 'Point' or Point_with_normal.
 // ///
 // /// @param first First point to add.
 // /// @param beyond Past-the-end point to add.
  //template < class InputIterator >
  //ReconstructionImplicitFunction(InputIterator first, InputIterator beyond);

  ///// Insert points.
  ///// Precondition: the value type of InputIterator must be 'Point' or Point_with_normal.
  ///// Return the number of inserted points.
  //template < class InputIterator >
  //int insert(InputIterator first, InputIterator beyond);

  /// Get the bounding box.
	Iso_cuboid_3 bounding_box() const;

  /// Get bounding sphere.
	Sphere bounding_sphere() const;

  /// Get the region of interest, ignoring the outliers.
  /// This method is used to define the OpenGL arcball sphere.
	Sphere region_of_interest() const;

  /// You should call compute_implicit_function() once when points insertion is over.
  /// Return false on error.
  /// TODO: add parameters to compute_implicit_function()?
	bool compute_implicit_function();

  /// Shift and orient the implicit function such that:
  /// - the implicit function = 0 for points / f() = contouring_value,
  /// - the implicit function < 0 inside the surface.
  ///
  /// Return the minimum value of the implicit function.
	FT set_contouring_value(FT contouring_value);

  /// [ImplicitFunction interface]
  ///
  /// Evaluate implicit function for any 3D point.
	FT operator() (Point p);

  /// Get point / the implicit function is minimum
	const Point& sink() const;

  /// Get average value of the implicit function over input vertices
	FT average_value_at_input_vertices() const;

  /// Get median value of the implicit function over input vertices
	FT median_value_at_input_vertices() const;

};


