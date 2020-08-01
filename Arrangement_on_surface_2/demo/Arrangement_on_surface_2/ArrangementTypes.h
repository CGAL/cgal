// Copyright (c) 2005,2012 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Alex Tsui <alextsui05@gmail.com>
//                 Saurabh Singh <ssingh@cs.iitr.ac.in>

#ifndef ARRANGEMENT_DEMO_TYPES_H
#define ARRANGEMENT_DEMO_TYPES_H

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Bbox_2.h>


#ifdef CGAL_USE_CORE
  #include <CGAL/CORE/BigRat.h>
  typedef CORE::BigRat                                  NT;
#elif CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
  typedef CGAL::Gmpq                                    NT;
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
  typedef CGAL::Quotient<CGAL::MP_Float>                NT;
#endif

// instead of
//typedef CGAL::Cartesian<NT>                           Kernel;
// workaround for VC++
struct Kernel : public CGAL::Cartesian<NT> {};

#include <QColor>

class Face_with_color : public CGAL::Arr_face_base
{
  QColor    m_color;
  bool      m_visited;

public:
  Face_with_color() :
      CGAL::Arr_face_base(), m_color(QColorConstants::White), m_visited(false)
  {
  }

  QColor color() const { return m_color; }
  void set_color(const QColor& c) { m_color = c; }
  bool visited() const{ return m_visited; }
  void set_visited(bool b) { m_visited = b; }
};

template <class Traits>
class Dcel :
  public CGAL::Arr_dcel_base<CGAL::Arr_vertex_base<typename Traits::Point_2>,
                             CGAL::Arr_halfedge_base<
                               typename Traits::X_monotone_curve_2>,
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

  // CREATION
  Dcel() {}
};

// Segments:
typedef CGAL::Arr_segment_traits_2<Kernel>              Seg_traits;
typedef Seg_traits::Curve_2                             Arr_seg_2;
typedef Seg_traits::X_monotone_curve_2                  Arr_xseg_2;
typedef Seg_traits::Point_2                             Arr_seg_point_2;
typedef Dcel<Seg_traits>                                Seg_dcel;
typedef CGAL::Arrangement_with_history_2<Seg_traits, Seg_dcel>
                                                        Seg_arr;
typedef Seg_arr::Halfedge                               Seg_halfedge;
typedef Seg_arr::Halfedge_handle                        Seg_halfedge_handle;
typedef Seg_arr::Face_handle                            Seg_face_handle;
typedef Seg_arr::Ccb_halfedge_circulator
  Seg_ccb_halfedge_circulator;
typedef Seg_arr::Hole_iterator                          Seg_holes_iterator;
typedef Seg_arr::Face_iterator                          Seg_face_iterator;
typedef std::list<Arr_seg_2*>                           Arr_seg_list;
typedef Arr_seg_list::const_iterator                    Arr_seg_const_iter;
typedef Arr_seg_list::iterator                          Arr_seg_iter;


// Polyline
typedef CGAL::Arr_polyline_traits_2<Seg_traits>         Pol_traits;

typedef Pol_traits::Curve_2                             Arr_pol_2;
typedef Pol_traits::X_monotone_curve_2                  Arr_xpol_2;

typedef Pol_traits::Point_2                             Arr_pol_point_2;
typedef Dcel<Pol_traits>                                Pol_dcel;
typedef CGAL::Arrangement_with_history_2<Pol_traits,
                                         Pol_dcel>      Pol_arr;
typedef Pol_arr::Halfedge_handle                        Pol_halfedge_handle;
typedef Pol_arr::Face_handle                            Pol_face_handle;
typedef Pol_arr::Ccb_halfedge_circulator
  Pol_ccb_halfedge_circulator;
typedef Pol_arr::Hole_iterator                          Pol_holes_iterator;
typedef Pol_arr::Halfedge                               Pol_halfedge;
typedef Pol_arr::Face_iterator                          Pol_face_iterator;

typedef std::list<Arr_pol_2*>                            Arr_pol_list;
typedef Arr_pol_list::const_iterator                     Arr_pol_const_iter;
typedef Arr_pol_list::iterator                           Arr_pol_iter;

// Bezier Curves typedefs

#ifdef CGAL_USE_CORE

#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Bezier_nt_traits;
typedef Bezier_nt_traits::Rational                    Bezier_rational;
typedef Bezier_nt_traits::Algebraic                   Bezier_algebraic;
typedef CGAL::Cartesian<Bezier_rational>              Bezier_rat_kernel;
typedef CGAL::Cartesian<Bezier_algebraic>             Bezier_alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<
  Bezier_rat_kernel, Bezier_alg_kernel, Bezier_nt_traits>
                                                      Bezier_traits;
typedef Bezier_rat_kernel::Point_2                    Bezier_rat_point;
typedef Bezier_traits::Curve_2                        Bezier_curve;
typedef Bezier_traits::X_monotone_curve_2             Bezier_X_monotone_curve;
typedef Bezier_traits::Point_2                        Bezier_point;
typedef Dcel<Bezier_traits>                           Bezier_dcel;
typedef CGAL::Arrangement_with_history_2<Bezier_traits, Bezier_dcel>
                                                      Bezier_arr;
#endif

// Conics

#ifdef CGAL_USE_CORE

#include <CGAL/CORE_algebraic_number_traits.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;

// instead of
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Conic_traits;
// workaround for VC++
/*
struct Conic_traits: public CGAL::Arr_conic_traits_2<Rat_kernel,
                                                     Alg_kernel,
                                                     Nt_traits>   {};
*/

typedef  Conic_traits::Curve_2                          Arr_conic_2;
typedef  Conic_traits::Rat_point_2                      Rat_point_2;
typedef  Conic_traits::Rat_segment_2                    Rat_segment_2;
typedef  Conic_traits::Rat_circle_2                     Rat_circle_2;
typedef  Conic_traits::Rat_line_2                       Rat_line_2;

typedef Conic_traits::X_monotone_curve_2                Arr_xconic_2;
typedef Conic_traits::Point_2                           Arr_conic_point_2;
typedef Dcel<Conic_traits>                              Conic_dcel;
typedef CGAL::Arrangement_with_history_2<Conic_traits, Conic_dcel>
                                                        Conic_arr;
typedef Conic_arr::Halfedge_handle                      Conic_halfedge_handle;
typedef Conic_arr::Face_handle                          Conic_face_handle;
typedef Conic_arr::Ccb_halfedge_circulator
  Conic_ccb_halfedge_circulator;
typedef Conic_arr::Hole_iterator                        Conic_holes_iterator;
//typedef CGAL::Arr_file_scanner<Conic_arr>                Arr_scanner;
typedef Conic_arr::Halfedge                             Conic_halfedge;
typedef Conic_arr::Face_iterator                        Conic_face_iterator;

typedef std::list<Arr_xconic_2*>                        Arr_xconic_list;
typedef Arr_xconic_list::const_iterator                 Arr_xconic_const_iter;
typedef Arr_xconic_list::iterator                       Arr_xconic_iter;

typedef Nt_traits::Integer                                Coefficient;
typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient> Alg_seg_traits;
typedef Dcel< Alg_seg_traits >                            Alg_seg_dcel;
typedef CGAL::Arrangement_with_history_2< Alg_seg_traits, Alg_seg_dcel >
                                                          Alg_seg_arr;
typedef Alg_seg_traits::Point_2                           Alg_seg_point_2;


#endif

// Linear:
typedef CGAL::Arr_linear_traits_2<Kernel>               Lin_traits;
typedef Lin_traits::Curve_2                             Arr_lin_2;
typedef Lin_traits::X_monotone_curve_2                  Arr_xlin_2;
typedef Lin_traits::Point_2                             Arr_lin_point_2;
typedef Dcel<Lin_traits>                                Lin_dcel;
typedef CGAL::Arrangement_with_history_2<Lin_traits, Lin_dcel>
                                                        Lin_arr;
typedef Lin_arr::Halfedge                               Lin_halfedge;
typedef Lin_arr::Halfedge_handle                        Lin_halfedge_handle;
typedef Lin_arr::Face_handle                            Lin_face_handle;
typedef Lin_arr::Ccb_halfedge_circulator
  Lin_ccb_halfedge_circulator;
typedef Lin_arr::Hole_iterator                          Lin_holes_iterator;
typedef Lin_arr::Face_iterator                          Lin_face_iterator;
typedef std::list<Arr_lin_2*>                           Arr_lin_list;
typedef Arr_lin_list::const_iterator                    Arr_lin_const_iter;
typedef Arr_lin_list::iterator                          Arr_lin_iter;

#endif // ARRANGEMENT_DEMO_TYPES_H
