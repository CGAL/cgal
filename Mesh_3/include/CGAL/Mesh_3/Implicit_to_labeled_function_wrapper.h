// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implicit_to_labeling_function_wrapper and
// Implicit_vector_to_labeling_function_wrapper class declaration
// and implementation.
//
// See classes description to have more information.
//******************************************************************************

#ifndef CGAL_MESH_3_IMPLICIT_TO_LABELED_FUNCTION_WRAPPER_H
#define CGAL_MESH_3_IMPLICIT_TO_LABELED_FUNCTION_WRAPPER_H

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

#define CGAL_DEPRECATED_HEADER "<CGAL/Mesh_3/Implicit_to_labeled_function_wrapper.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Implicit_to_labeling_function_wrapper.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <vector>

namespace CGAL {

namespace Mesh_3 {

#include <CGAL/config.h>

/**
 * @class Implicit_to_labeled_function_wrapper
 *
 * This class is designed to wrap an implicit function which describes a domain
 * by [p is inside if f(p)<0] to a function which takes its values into {0,1}.
 * f(p)=0 means that p is outside the domain.
 */
template<class Function_, class BGT>
class Implicit_to_labeled_function_wrapper
{
public:
  // Types
  typedef int                     return_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Implicit_to_labeled_function_wrapper(const Function_& f)
    : r_f_(f) {}

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Implicit_to_labeled_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, const bool = true) const
  {
    return ( (r_f_(p)<0) ? 1 : 0 );
  }

private:
  /// Function to wrap
  const Function_& r_f_;

};  // end class Implicit_to_labeled_function_wrapper



/**
 * \deprecated
 *
 * @class Implicit_vector_to_labeled_function_wrapper
 *
 * Wraps a set of implicit function [f1,f2,...] to one function F which
 * takes its values into N.
 *
 * Let p be a point.
 * F(p) = 0b000000(f2(p)<0)(f1(p)<0)
 *
 * It can handle at most 8 functions.
 */
template<class Function_, class BGT>
class Implicit_vector_to_labeled_function_wrapper
{
public:
  // Types
  typedef int                       return_type;
  typedef std::vector<Function_*>   Function_vector;
  typedef typename BGT::Point_3     Point_3;

  /// Constructor
  Implicit_vector_to_labeled_function_wrapper(const std::vector<Function_*>& v)
    : function_vector_(v) {}

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Implicit_vector_to_labeled_function_wrapper() {}

  /// Operator ()
  return_type operator()(const Point_3& p, const bool = true) const
  {
    int nb_func = static_cast<int>(function_vector_.size());
    if ( nb_func > 8 )
    {
      CGAL_error_msg("We support at most 8 functions !");
    }

    char bits = 0;
    for ( int i = 0 ; i < nb_func ; ++i )
    {
      // Insert value into bits : we compute fi(p) and insert result at
      // bit i of bits
      bits |= ( ((*function_vector_[i])(p) < 0) << i );
    }

    return ( static_cast<return_type>(bits) );
  }

private:
  /// Functions to wrap
  const Function_vector function_vector_;

};  // end class Implicit_to_labeled_function_wrapper

}  // end namespace Mesh_3

}  // end namespace CGAL



#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_MESH_3_IMPLICIT_TO_LABELED_FUNCTION_WRAPPER_H
