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

// CGAL headers.
#include <CGAL/number_utils.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Finds a square root of the provided value of the type `Kernel::FT` by first converting it to the double type and then taking the square root using the `CGAL::sqrt()` function.
template<class Kernel> 
	class Sqrt
{
	typedef typename Kernel::FT Scalar;

public:
    Scalar operator()(const Scalar &value) const
    { 
    	return Scalar(CGAL::sqrt(CGAL::to_double(value)));
    }
};

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Barycentric_traits_2` is a traits class for all two-dimensional barycentric coordinates. It refines `Kernel` by adding the function `AlgebraicStructureTraits_::Sqrt`.
 * This class is parameterized by the `Kernel` class.
*/

template<class Kernel>
	class Barycentric_traits_2 : public Kernel
{

public:

   	/// \name Types
    /// @{

	typedef Sqrt<Kernel> Sqrt;

	/// @}

    /// \name Operations
    /// @{
 
	inline const Sqrt sqrt_object() const
	{
		return Sqrt();
	}

	/// @}
};

} // namespace Barycentric_coordinates

} //namespace CGAL

#endif // CGAL_BARYCENTRIC_TRAITS_2_H