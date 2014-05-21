// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

/*!
  \file Barycentric_traits_2.h
*/

#ifndef CGAL_BARYCENTRIC_TRAITS_2_H
#define CGAL_BARYCENTRIC_TRAITS_2_H

// STL headers.
#include <vector>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Barycentric_traits_2` is a default traits class for all two-dimensional barycentric coordinates, and it provides all the necessary types and functions used internally in the coordinate classes.
 * This class is parameterized by the `Kernel` class.
*/

template <class Kernel>
	class Barycentric_traits_2 : public Kernel
{

public:

	/// \name Types
    /// @{

	/// Number type.
	typedef typename Kernel::FT           FT;

	/// Type of the used kernel.
	typedef Kernel                         K;

	/// Type of a general point.
	typedef typename Kernel::Point_2 Point_d;

	/// Type of 3D point.
	typedef typename Kernel::Point_3 Point_3;

	/// Type of 2D point.
	typedef typename Kernel::Point_2 Point_2; 

	/// @}    
  
    /// \name Projecting a point
    /// @{

	/// This function is required to convert a query point of the user-defined type to the type `CGAL::Point_2` used internally in all coordinate classes.
	inline const Point_2& project(const Point_d &query_point) const 
	{ 
		return query_point; 
	} 

	/// @}

    /// \name Projecting a sequence of points
    /// @{

	/// This function is required to convert a sequence of the polygon's vertices of the user-defined type to the sequence of points of the type `CGAL::Point_2`.
	/// This new sequence of points is stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> used internally in all coordinate classes.
	template<class InputIterator> 
		inline const std::vector<Point_2> project(const InputIterator &first_vertex, const InputIterator &last_vertex) const
	{
			std::vector<Point_2> vertex;
			vertex.resize(last_vertex - first_vertex);

			int i = 0;
			for(InputIterator it = first_vertex; it != last_vertex; ++it) {
				vertex[i] = project(*it);
				++i;
			}

			return vertex;
	}

	/// @}
};

} // namespace Barycentric_coordinates

} //namespace CGAL

#endif // CGAL_BARYCENTRIC_TRAITS_2_H
