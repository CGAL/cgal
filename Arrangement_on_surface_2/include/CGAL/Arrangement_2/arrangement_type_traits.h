// Copyright (c) 2005,2009,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Ophir Setter <ophirset@post.tau.ac.il>


/*!
  \file   arrangement_type_traits.h
  \brief  The file contains meta-function related to the arrangement pakcage.
  Specifically, it contains the meta-function is_arrangement_2 that determines
  whether a given type is an arrangement.
*/

#ifndef CGAL_ARRANGEMENT_TYPE_TRAITS_H
#define CGAL_ARRANGEMENT_TYPE_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <boost/type_traits/integral_constant.hpp>

namespace CGAL
{

// The current implementation is not general. The following section describes
// the general solution that needs volunteers that wants to learn more generic
// programming.
// We detect "arrangements" by checking if they define types as specified and
// functions as specified. (For example, if they contain a Geometry_traits_2 
// type.
//
// 1) Write a meta-function that detects whether a type T has a nested type.
//    Consult Modern C++ Design Chapter 2.7. should look something like:
// 
// template <class T>
// class does_contain_Geometry_traits_2
// {
//   typedef char        True;
//   class False { char dummy[2]; };
  
//   True Test(typename T::Geometry_traits_2);
//   False Test(...);
//   [Continue code]
// public:
  // [Continue code]
// };
//
// 2) Make this class a macro.
// 3) Can you detect functions in a class T (consult boost::type_traits.
// 4) Use boost::mpl to have is_arrangement_2 derive from the correct type 
//    (true_type or false_type) according to the meta-functions you wrote.
// 5) Optional: Check boost type_traits for platform-specific problems and 
//    solutions.


// In the meanwhile we use a default implementation.
template <class T>
class is_arrangement_2 : public boost::false_type
{};

//--------------------------------  Arrangement_2
// forward declaration first.
template <class GeomTraits_, class DCEL_>
class Arrangement_2;

// specialization
template <class GeomTraits_, class DCEL_>
class is_arrangement_2< 
  Arrangement_2<GeomTraits_, DCEL_>
> : public boost::false_type
{};


//--------------------------------  Arrangement_on_surface_2
// forward declaration first.
template <class GeomTraits_, class TopTraits_>
class Arrangement_on_surface_2;

// specialization
template <class GeomTraits_, class TopTraits_>
class is_arrangement_2< 
  Arrangement_on_surface_2<GeomTraits_, TopTraits_>
> : public boost::true_type
{};

} // namespace CGAL

#endif // CGAL_ARRANGEMENT_TYPE_TRAITS_H
