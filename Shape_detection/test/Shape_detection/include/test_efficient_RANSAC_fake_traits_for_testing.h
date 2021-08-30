// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_TRAITS_COMPILATION_TEST_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_TRAITS_COMPILATION_TEST_H

#include <CGAL/Origin.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL {

  struct Point__3{
    Point__3(){}
    Point__3(Origin const&){}
    Point__3(double, double, double){}
  };

  struct Vector__3 {
    Vector__3(){}
    Vector__3(Null_vector const&){}
    Vector__3(double, double, double){}
  };

  namespace Shape_detection {
  template <class Gt,
            class InputRange,
            class InputPointMap,
            class InputNormalMap>
  struct Efficient_RANSAC_traits {

    typedef Efficient_RANSAC_traits<Gt, InputRange, InputPointMap, InputNormalMap> Self;

    typedef double FT;
    typedef Point__3 Point_3;
    typedef Vector__3 Vector_3;

    struct Sphere_3 {};

    struct Line_3 {};

    struct Circle_2 {};

    struct Vector_2 {};

    struct Plane_3 {};

    struct Point_2 {};

    typedef InputRange Input_range;

    typedef InputPointMap Point_map;

    typedef InputNormalMap Normal_map;

    class Search_traits {

    public:

      typedef Dimension_tag<3> Dimension;
      typedef typename Self::Point_3 Point_d;

      struct Cartesian_const_iterator_d
      {
        typedef Cartesian_const_iterator_d      self;

        typedef std::random_access_iterator_tag iterator_category;
        typedef FT                              value_type;
        typedef std::ptrdiff_t                  difference_type;
        typedef const value_type*               pointer;
        typedef const value_type&               reference;

        self& operator++() { return *this; }
        self  operator++(int) { return *this; }
        self& operator--() { return *this; }
        self  operator--(int) { return *this; }

        self& operator+=(difference_type) { return *this; }
        self& operator-=(difference_type) { return *this; }
        self operator+(difference_type) const { return *this; }
        self operator-(difference_type) const { return *this; }

        difference_type operator-(self ) const { return 0; }

        value_type operator*() const { return value_type(); }
        value_type operator[](difference_type ) const { return value_type(); }

        bool operator==(const self& ) const { return true; }
        bool operator!=(const self& ) const { return true; }
        bool operator<(const self& ) const  { return true; }
      };

      struct Construct_cartesian_const_iterator_d{
        typedef Point_3 Point_d;
        typedef Cartesian_const_iterator_d result_type;

        Cartesian_const_iterator_d operator()(Point_d const& ) const { return Cartesian_const_iterator_d(); }
        Cartesian_const_iterator_d operator()(Point_d const& , int) const { return Cartesian_const_iterator_d(); }
      };
      struct Iso_box_d{};
      struct Sphere_d{};
      struct Construct_iso_box_d{};
      struct Construct_min_vertex_d{};
      struct Construct_max_vertex_d{};
      struct Construct_center_d{};
      struct Compute_squared_radius_d{};
      typedef typename Self::FT FT;

      Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const {
         return Construct_cartesian_const_iterator_d();
      }
    };

    struct Construct_point_3
    {
      Point_3 operator()(FT, FT, FT) { return Point_3(); }
    };
    Construct_point_3 construct_point_3_object() const
    { return Construct_point_3(); }

    struct Construct_point_2
    {
      Point_2 operator()(FT, FT) { return Point_2(); }
    };
    Construct_point_2 construct_point_2_object() const
    { return Construct_point_2(); }

    struct Construct_vector_3
    {
      Vector_3 operator()(Null_vector const&) { return NULL_VECTOR; }
      Vector_3 operator()(Point_3 const&, Point_3 const&) { return Vector_3(); }
      Vector_3 operator()(Origin const&, Point_3 const&) { return Vector_3(); }
      Vector_3 operator()(Line_3 const&) { return Vector_3(); }
    };
    Construct_vector_3 construct_vector_3_object() const
    { return Construct_vector_3(); }

    struct Construct_vector_2
    {
      Vector_2 operator()(Null_vector const&) { return NULL_VECTOR; }
      Vector_2 operator()(Point_2 const&, Point_2 const&) { return Vector_2(); }
    };
    Construct_vector_2 construct_vector_2_object() const
    { return Construct_vector_2(); }

    struct Construct_sphere_3
    {
      Sphere_3 operator()(Point_3 const&, FT) { return Sphere_3(); }
    };
    Construct_sphere_3 construct_sphere_3_object() const
    { return Construct_sphere_3(); }

    struct Construct_circle_2
    {
      Circle_2 operator()(Point_2 const&, Point_2 const&, Point_2 const&) { return Circle_2(); }
    };
    Construct_circle_2 construct_circle_2_object() const
    { return Construct_circle_2(); }

    struct Construct_line_3
    {
      Line_3 operator()(Point_3 const&, Vector_3 const&) { return Line_3(); }
    };
    Construct_line_3 construct_line_3_object() const
    { return Construct_line_3(); }

    struct Construct_point_on_3
    {
      Point_3 operator()(Line_3 const&, int) { return Point_3(); }
    };
    Construct_point_on_3 construct_point_on_3_object() const
    { return Construct_point_on_3(); }

    struct Compute_x_2
    {
      FT operator()(Point_2 const&) const { return 0; }
      FT operator()(Vector_2 const&) const { return 0; }
    };
    Compute_x_2 compute_x_2_object() const
    { return Compute_x_2(); }

    struct Compute_y_2
    {
      FT operator()(Point_2 const&) const { return 0; }
      FT operator()(Vector_2 const&) const { return 0; }
    };
    Compute_y_2 compute_y_2_object() const
    { return Compute_y_2(); }

    struct Compute_x_3
    {
      FT operator()(Point_3 const&) const { return 0; }
      FT operator()(Vector_3 const&) const { return 0; }
    };
    Compute_x_3 compute_x_3_object() const
    { return Compute_x_3(); }

    struct Compute_y_3
    {
      FT operator()(Point_3 const&) const { return 0; }
      FT operator()(Vector_3 const&) const { return 0; }
    };
    Compute_y_3 compute_y_3_object() const
    { return Compute_y_3(); }

    struct Compute_z_3
    {
      FT operator()(Point_3 const&) const { return 0; }
      FT operator()(Vector_3 const&) const { return 0; }
    };

    Compute_z_3 compute_z_3_object() const
    { return Compute_z_3(); }

    struct Compute_squared_length_3
    { FT operator()(Vector_3 const&) const { return 0; } };
    Compute_squared_length_3 compute_squared_length_3_object() const
    { return Compute_squared_length_3(); }

    struct Compute_squared_length_2
    { FT operator()(Vector_2 const&) const { return 0; } };
    Compute_squared_length_2 compute_squared_length_2_object() const
    { return Compute_squared_length_2(); }

    struct Construct_scaled_vector_3
    {
      Vector_3 operator()(Vector_3 const&, FT) { return Vector_3(); }
    };
    Construct_scaled_vector_3 construct_scaled_vector_3_object() const
    { return Construct_scaled_vector_3(); }

    struct Construct_sum_of_vectors_3
    {
      Vector_3 operator()(Vector_3 const&, Vector_3 const&) { return Vector_3(); }
    };
    Construct_sum_of_vectors_3 construct_sum_of_vectors_3_object() const
    { return Construct_sum_of_vectors_3(); }

    struct Construct_translated_point_3
    { Point_3 operator()(Point_3 const&, Vector_3 const&) const { return Point_3(); } };
    Construct_translated_point_3 construct_translated_point_3_object() const
    { return Construct_translated_point_3(); }

    struct Compute_scalar_product_3
    { FT operator()(Vector_3 const&, Vector_3 const&) const { return 0; } };
    Compute_scalar_product_3 compute_scalar_product_3_object() const
    { return Compute_scalar_product_3(); }

    struct Construct_cross_product_vector_3
    { Vector_3 operator()(Vector_3 const&, Vector_3 const&) const { return Vector_3(); } };
    Construct_cross_product_vector_3 construct_cross_product_vector_3_object() const
    { return Construct_cross_product_vector_3(); }

    struct Construct_center_3
    { Point_3 operator()(Sphere_3 const&) const { return Point_3(); } };
    Construct_center_3 construct_center_3_object() const
    { return Construct_center_3(); }

    struct Compute_squared_radius_3
    { FT operator()(Sphere_3 const&) const { return 0; } };
    Compute_squared_radius_3 compute_squared_radius_3_object() const
    { return Compute_squared_radius_3(); }

    struct Compute_squared_radius_2
    { FT operator()(Circle_2 const&) const { return 0; } };
    Compute_squared_radius_2 compute_squared_radius_2_object() const
    { return Compute_squared_radius_2(); }

    struct Construct_center_2
    { Point_2 operator()(Circle_2 const&) const { return Point_2(); } };
    Construct_center_2 construct_center_2_object() const
    { return Construct_center_2(); }

    struct Collinear_2
    {
      bool operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
      { return false; }
    };
    Collinear_2 collinear_2_object() const
    { return Collinear_2(); }

    /* struct Compute_squared_distance_3
    { FT operator()(Point_3 const&, Point_3 const&) const { return 0; } };
    Compute_squared_distance_3 compute_squared_distance_3_object() const
    { return Compute_squared_distance_3(); } */
    ///
    Efficient_RANSAC_traits(const Gt& gt =  Gt())
      : m_gt(gt) {}

  private:
    Gt m_gt;
  };
}

}

#endif // CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_TRAITS_COMPILATION_TEST_H
