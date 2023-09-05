// Copyright (c) 2013 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel   <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POINT_LOCATION_RESULT_H
#define CGAL_ARR_POINT_LOCATION_RESULT_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/assertions.h>

#include <optional>
#include <variant>

namespace CGAL {

template <typename Arrangement_>
struct Arr_point_location_result {
  typedef Arrangement_                                   Arrangement_2;

  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle      Face_const_handle;

  typedef typename std::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>     Type;
  typedef Type                                           type;

  // This function returns a result_type constructor
  // to generate return values.
  // In theory a one parameter variant could be returned, but this _could_
  // lead to conversion overhead, and so we rather go for the real type.
  // Overloads for empty returns are also provided.
  template <typename T>
  static
  inline Type make_result(T t) { return Type(t); }

  static
  inline std::optional<Type> empty_optional_result()
  { return std::optional<Type>(); }

  template <typename T>
  static
  inline const T* assign(const Type* obj) { return std::get_if<T>(obj); }

  //this one is only to remove warnings in functions
  static
  inline Type default_result(){
    CGAL_error_msg("This functions should have never been called!");
    return Type();
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
