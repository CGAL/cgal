// Copyright (c) 2013 Technical University Braunschweig (Germany).
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
// Author(s):  Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>

#ifndef CGAL_TEST_MODEL_METHODS_H
#define CGAL_TEST_MODEL_METHODS_H

#include <CGAL/basic.h>
#include <CGAL/test_utils.h>
#include <cassert>

namespace CGAL {

template <class _Visibility_2>
bool test_is_attached(_Visibility_2 visibility) {
    return visibility.is_attached();
}

template <class _Visibility_2>
void test_model_methods(_Visibility_2 &visibility, 
                      const typename _Visibility_2::Input_arrangement_2 &arr) {

  // Check concept obediance
	typedef _Visibility_2 Visibility_2;
	typedef typename Visibility_2::Input_arrangement_2     Input_arrangement_2;
	typedef typename Visibility_2::Output_arrangement_2    Output_arrangement_2;
  typedef typename Visibility_2::Regularization_tag      Regularization_tag;
  typedef typename Visibility_2::Supports_general_polygon_tag
                                                  Supports_general_polygon_tag;
  typedef typename Visibility_2::Supports_simple_polygon_tag 
                                                  Supports_simple_polygon_tag;

  assert(false == visibility.is_attached());
  visibility.attach(arr);
  assert(true == visibility.is_attached());
  visibility.detach();
  assert(false == visibility.is_attached());
  visibility.attach(arr);
  assert(true == visibility.is_attached());
  assert(true == (CGAL::test_are_equal<Input_arrangement_2>(arr, 
                                                            visibility.arr())));
}

} // end CGAL namespace
#endif