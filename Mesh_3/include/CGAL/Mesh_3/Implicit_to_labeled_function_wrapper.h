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

#define CGAL_DEPRECATED_HEADER "<CGAL/Mesh_3/Implicit_to_labeled_function_wrapper.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Implicit_to_labeling_function_wrapper.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>

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
  : public CGAL::Implicit_to_labeling_function_wrapper<Function_, BGT>
{
  typedef CGAL::Implicit_to_labeling_function_wrapper<Function_, BGT> Base;
public:
  /// Constructor
  Implicit_to_labeled_function_wrapper(const Function_& f)
    : Base(f)
  {}
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
  : public CGAL::Implicit_vector_to_labeling_function_wrapper<Function_, BGT>
{
  typedef CGAL::Implicit_vector_to_labeling_function_wrapper<Function_,
                                                             BGT> Base;
public:
  /// Constructor
  Implicit_vector_to_labeled_function_wrapper(const std::vector<Function_*>& v)
    : Base(v) {}

};  // end class Implicit_to_labeled_function_wrapper

}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMPLICIT_TO_LABELED_FUNCTION_WRAPPER_H
