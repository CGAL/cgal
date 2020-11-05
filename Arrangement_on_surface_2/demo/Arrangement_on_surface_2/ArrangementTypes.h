// Copyright (c) 2005, 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Alex Tsui <alextsui05@gmail.com>
//            Saurabh Singh <ssingh@cs.iitr.ac.in>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_TYPES_H
#define ARRANGEMENT_DEMO_TYPES_H

#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

#ifdef CGAL_USE_CORE
  #include <CGAL/CORE_algebraic_number_traits.h>
  #include <CGAL/Arr_conic_traits_2.h>
  #include <CGAL/Arr_algebraic_segment_traits_2.h>
  #include <CGAL/Arr_Bezier_curve_traits_2.h>
  #include <CGAL/Arr_rational_function_traits_2.h>
  #include <CGAL/Algebraic_kernel_d_1.h>
  #include <CGAL/CORE_BigRat.h>
#endif

#include "RationalTypes.h"

#include <QColor>

// avoid polluting global namespace
// caused multiple bugs before!
namespace demo_types
{
struct DemoTypes
{
class Face_with_color : public CGAL::Arr_face_base
{
  QColor    m_color;
  bool      m_visited;

public:
  Face_with_color() :
      CGAL::Arr_face_base(), m_color(::Qt::white), m_visited(false)
  {
  }

  QColor color() const { return m_color; }
  void set_color(const QColor& c) { m_color = c; }
  bool visited() const{ return m_visited; }
  void set_visited(bool b) { m_visited = b; }
};

template <class Traits>
class Dcel :
    public CGAL::Arr_dcel_base<
      CGAL::Arr_vertex_base<typename Traits::Point_2>,
      CGAL::Arr_halfedge_base<typename Traits::X_monotone_curve_2>,
      Face_with_color>
{
public:
   /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  struct rebind
  {
    typedef Dcel<T> other;
  };
};

// use same rational type across arrangements
// less specializations, and makes PointSnapper untemplated
typedef typename RationalTypes::Rational                Rational;
typedef typename RationalTypes::Rat_kernel              Rat_kernel;
typedef typename RationalTypes::Rat_point_2             Rat_point_2;

// Segments:
typedef CGAL::Arr_segment_traits_2<Rat_kernel>          Seg_traits;
typedef Dcel<Seg_traits>                                Seg_dcel;
typedef CGAL::Arrangement_with_history_2<Seg_traits, Seg_dcel>
                                                        Seg_arr;

// Polyline
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Pol_traits;
typedef Dcel<Pol_traits>                                Pol_dcel;
typedef CGAL::Arrangement_with_history_2<Pol_traits,
                                         Pol_dcel>      Pol_arr;

// Linear:
typedef CGAL::Arr_linear_traits_2<Rat_kernel>           Lin_traits;
typedef Dcel<Lin_traits>                                Lin_dcel;
typedef CGAL::Arrangement_with_history_2<Lin_traits, Lin_dcel>
                                                        Lin_arr;

#ifdef CGAL_USE_CORE
  typedef CGAL::CORE_algebraic_number_traits            Core_nt_traits;
  typedef Core_nt_traits::Integer                       Core_integer;
  typedef Core_nt_traits::Rational                      Core_rational;
  typedef Core_nt_traits::Algebraic                     Core_algebraic;
  typedef CGAL::Cartesian<Core_rational>                Core_rat_kernel;
  typedef CGAL::Cartesian<Core_algebraic>               Core_alg_kernel;

  // Conics
  typedef CGAL::Arr_conic_traits_2<
    Core_rat_kernel, Core_alg_kernel, Core_nt_traits>
                                                        Conic_traits;

  typedef Dcel<Conic_traits>                            Conic_dcel;
  typedef CGAL::Arrangement_with_history_2<Conic_traits, Conic_dcel>
                                                        Conic_arr;
  // Algebraic
  typedef CGAL::Arr_algebraic_segment_traits_2<Core_integer>
                                                        Alg_seg_traits;
  typedef Dcel<Alg_seg_traits>                          Alg_seg_dcel;
  typedef CGAL::Arrangement_with_history_2<Alg_seg_traits, Alg_seg_dcel>
                                                        Alg_seg_arr;
  // Bezier
  typedef CGAL::Arr_Bezier_curve_traits_2<
    Core_rat_kernel, Core_alg_kernel, Core_nt_traits>
                                                        Bezier_traits;
  typedef Dcel<Bezier_traits>                           Bezier_dcel;
  typedef CGAL::Arrangement_with_history_2<Bezier_traits, Bezier_dcel>
                                                        Bezier_arr;

  // Rational functions
  typedef CGAL::Algebraic_kernel_d_1<Core_integer>      AK1;
  typedef CGAL::Arr_rational_function_traits_2<AK1>     Rational_traits;
  typedef Dcel<Rational_traits>                         Rational_dcel;
  typedef CGAL::Arrangement_with_history_2<Rational_traits, Rational_dcel>
                                                        Rational_arr;
#endif
};
}

#ifdef CGAL_USE_CORE
  #define ARRANGEMENT_DEMO_SPECIALIZE_ARR_CORE(class_name)                     \
    template class class_name<demo_types::DemoTypes::Conic_arr>;               \
    template class class_name<demo_types::DemoTypes::Alg_seg_arr>;             \
    template class class_name<demo_types::DemoTypes::Bezier_arr>;              \
    template class class_name<demo_types::DemoTypes::Rational_arr>;
  #define ARRANGEMENT_DEMO_SPECIALIZE_TRAITS_CORE(class_name)                  \
    template class class_name<demo_types::DemoTypes::Conic_traits>;            \
    template class class_name<demo_types::DemoTypes::Alg_seg_traits>;          \
    template class class_name<demo_types::DemoTypes::Bezier_traits>;           \
    template class class_name<demo_types::DemoTypes::Rational_traits>;
#else
  #define ARRANGEMENT_DEMO_SPECIALIZE_ARR_CORE(class_name)
  #define ARRANGEMENT_DEMO_SPECIALIZE_TRAITS_CORE(class_name)
#endif

#define ARRANGEMENT_DEMO_SPECIALIZE_ARR(class_name)                            \
  template class class_name<demo_types::DemoTypes::Seg_arr>;                   \
  template class class_name<demo_types::DemoTypes::Pol_arr>;                   \
  template class class_name<demo_types::DemoTypes::Lin_arr>;                   \
  ARRANGEMENT_DEMO_SPECIALIZE_ARR_CORE(class_name)

#define ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(class_name)                         \
  template class class_name<demo_types::DemoTypes::Seg_traits>;                \
  template class class_name<demo_types::DemoTypes::Pol_traits>;                \
  template class class_name<demo_types::DemoTypes::Lin_traits>;                \
  ARRANGEMENT_DEMO_SPECIALIZE_TRAITS_CORE(class_name)

#endif // ARRANGEMENT_DEMO_TYPES_H
