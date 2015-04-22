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

      \tparam Gt A \cgal Kernel. Geometric traits class. It must provide `Gt::FT`,
	         `Gt::Point_3` and `Gt::Vector_3`. The type `Gt::FT` must be a floating
			 point number type like `double` or `float`. Additionally,
             `Gt::Line_3`, `Gt::Plane_3`, `Gt::Sphere_3`, `Gt::Circle_2`,
			 `Gt::Point_2` are required depending on the shapes registered
			 for detection.

      \tparam InputRange is a model of `Range` with random access iterators, 
              providing input points and normals through the following two property maps.

      \tparam InputPointMap is a model of `ReadablePropertyMap` with `std::iterator_traits<Input_range::iterator>::%value_type` as key type and `Geom_traits::Point_3` as value type.


      \tparam InputNormalMap is a model of `ReadablePropertyMap` with `std::iterator_traits<Input_range::iterator>::%value_type` as key type and `Geom_traits::Vector_3` as value type.
    */
  template <class Gt,
            class InputRange,
            class InputPointMap,
            class InputNormalMap>
  struct Efficient_RANSAC_traits {
    ///
    typedef typename Gt::FT FT;
    ///
    typedef typename Gt::Vector_3 Vector_3;
    ///
    typedef typename Gt::Point_3 Point_3;
    ///
    typedef typename Gt::Sphere_3 Sphere_3;
    ///
    typedef typename Gt::Line_3 Line_3;
    ///
    typedef typename Gt::Circle_2 Circle_2;
    ///
    typedef typename Gt::Plane_3 Plane_3;
    ///
    typedef typename Gt::Point_2 Point_2;
    ///
    typedef InputRange Input_range;
    ///
    typedef InputPointMap Point_map;
    ///
    typedef InputNormalMap Normal_map;
    ///
    typedef typename CGAL::Search_traits_3<typename Gt> Search_traits;
  };

} } // end of namespace CGAL::Shape_detection_3

#endif // CGAL_SHAPE_DETECTION_3_EFFICIENT_RANSAC_TRAITS_H
