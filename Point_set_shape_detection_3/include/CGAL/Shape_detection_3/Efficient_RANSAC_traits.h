// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_3_EFFICIENT_RANSAC_TRAITS_H
#define CGAL_SHAPE_DETECTION_3_EFFICIENT_RANSAC_TRAITS_H


namespace CGAL {
  namespace Shape_detection_3 {
    /*!
      \ingroup PkgPointSetShapeDetection3
      \brief %Default traits class to use the shape detection class `Efficient_RANSAC`.
      \cgalModels `EfficientRANSACTraits`

      \tparam Gt Geometric traits class. It must provide `Gt::FT`, `Gt::Point_3` and `Gt::Vector_3`.
             `Gt::FT` must be a floating point number type like `double` or `float`.

      \tparam InputIt is a model of RandomAccessIterator

      \tparam Ppmap is a model of `ReadablePropertyMap`
              `key_type = InputIt` and `value_type = Gt::Point_3`.

      \tparam Npmap is a model of `ReadablePropertyMap`
              `key_type = InputIt` and `value_type = Gt::Vector_3`.
    */
  template <class Gt,
            class InputIt,
            class Ppmap,
            class Npmap>
  struct Efficient_RANSAC_traits {
    ///
    typedef Gt Geom_traits;
    ///
    typedef InputIt Input_iterator;
    ///
    typedef Ppmap Point_pmap;
    ///
    typedef Npmap Normal_pmap;
  };

} } // end of namespace CGAL::Shape_detection_3

#endif // CGAL_SHAPE_DETECTION_3_EFFICIENT_RANSAC_TRAITS_H
