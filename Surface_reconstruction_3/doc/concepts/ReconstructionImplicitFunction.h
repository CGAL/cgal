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
/// The surface is the zero-level of the implicit function.
///
/// The Surface Mesh Generation package makes copies of implicit functions,
/// thus such a class must be lightweight and stateless.
///
/// @heading Design Pattern:
/// A model of ReconstructionImplicitFunction is a
/// Strategy [GHJV95]: it implements a strategy of surface mesh reconstruction.
///
/// @heading Has Models:
/// - Poisson_implicit_function<PoissonDelaunayTriangulation_3>
/// - APSS_implicit_function<Gt>

class ReconstructionImplicitFunction : public Surface_mesher::ImplicitFunction
{
// Public types
public:

  typedef xxx Geom_traits; ///< Kernel's geometric traits

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  typedef xxx Point_with_normal;                     ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.
  typedef typename Geom_traits::Vector_3 Vector;

// Public methods
public:

  /// Create an implicit function from a point set.
  ///
  /// Precondition: the value type of InputIterator must be convertible to Point_with_normal.
  ///
  /// @param first First point to add.
  /// @param beyond Past-the-end point to add.
  template < class InputIterator >
  ReconstructionImplicitFunction(InputIterator first, InputIterator beyond);

  /// Get the surface's bounding sphere.
  Sphere bounding_sphere() const;

  /// Get the region of interest, ignoring the outliers.
  /// This method is used to define the OpenGL arcball sphere.
  Sphere region_of_interest() const;

  /// You should call compute_implicit_function() once when points insertion is over.
  /// Return false on error.
  bool compute_implicit_function();

  /// [ImplicitFunction interface]
  ///
  /// Evaluate implicit function for any 3D point.
  FT operator()(const Point& p) const;

  /// Get point inside the surface.
  Point get_inner_point() const;
};


