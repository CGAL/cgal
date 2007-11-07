// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://gankit@scm.gforge.inria.fr/svn/cgal/trunk/Principal_component_analysis/include/CGAL/tags.h $
// $Id: tags.h 37882 2007-04-03 15:15:30Z spion $
//
// Author(s)     : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_TAGS_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_TAGS_H

CGAL_BEGIN_NAMESPACE

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::These Tags are meant to represent the type of
//::geometry we are dealing with. The CGAL kernel
//::provides only for circles, triangles, spheres,
//: rectangles and so on. But these can be solid,
//::hollow, just outlined etc. These tags resolve
//::between these different types of the same geometry
//::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Fully filled geometry object, such as disks in 2d and balls in 3d
struct PCA_dimension_3_tag {};

// Empty geometry objects, such as triangles in 2d and hollow cubes with just the planar facets in 3d
struct PCA_dimension_2_tag {};

// Only relevant for 3d objects, outlines of 3d geometry objects formed by just lines joing appropriate vetrices, such as outline of a cuboid.
struct PCA_dimension_1_tag {};

// For the vertices of objects.
struct PCA_dimension_0_tag {};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Struct to denote dimension compile time decisions
//::::::::::::::::::::::::::::::::::::::::::::::::::::::

template <class Value_type /* type of fitted objects */,
          class K>
struct PCA_default_dimension{};

template < typename K >
struct PCA_default_dimension<typename K::Point_3> {
  typedef CGAL::PCA_dimension_0_tag Tag;
};

template < typename K >
struct PCA_default_dimension<typename K::Segment_3> {
  typedef CGAL::PCA_dimension_1_tag Tag;
};

template < typename K >
struct PCA_default_dimension<typename K::Triangle_3> {
  typedef CGAL::PCA_dimension_2_tag Tag;
};


CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_TAGS_H
