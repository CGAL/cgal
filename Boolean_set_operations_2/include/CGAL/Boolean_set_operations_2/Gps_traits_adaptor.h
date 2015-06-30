// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_TRAITS_ADAPTOR_H
#define CGAL_GPS_TRAITS_ADAPTOR_H

#include <CGAL/enum.h>
#include <CGAL/assertions.h>

namespace CGAL {

template <typename Traits_>
class Gps_traits_adaptor : public Traits_ {
  typedef Traits_                     Base;
  typedef Gps_traits_adaptor<Base>    Self;

public:
  typedef typename Base::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Base::Point_2                 Point_2;
  typedef typename Base::Compare_xy_2            Compare_xy_2;
  typedef typename Base::Compare_y_at_x_right_2  Compare_y_at_x_right_2;
  typedef typename Base::Compare_endpoints_xy_2  Compare_endpoints_xy_2;
  typedef typename Base::Construct_min_vertex_2  Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2  Construct_max_vertex_2;

  /*! Default constructor. */
  Gps_traits_adaptor() : Base() {}

  /*! Constructor from a base-traits class. */
  Gps_traits_adaptor(const Base& traits) : Base (traits) {}

  class Construct_vertex_2 {
  protected:
    typedef Gps_traits_adaptor<Base>        Traits_adaptor;

    /*! The traits (in case it has state) */
    const Traits_adaptor& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Construct_vertex_2(const Traits_adaptor& traits) : m_traits(traits) {}

    friend class Gps_traits_adaptor<Base>;

  public:
    Point_2 operator()(const X_monotone_curve_2& cv, int i) const
    {
      const Base& tr = m_traits;
      Compare_endpoints_xy_2 cmp_endpoints =
        tr.compare_endpoints_xy_2_object();
      Construct_min_vertex_2 ctr_min_v = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2 ctr_max_v = tr.construct_max_vertex_2_object();
      i %= 2;
      if (i == 0) {
        // return the source
        if (cmp_endpoints(cv) == SMALLER) return (ctr_min_v(cv));

        return (ctr_max_v(cv));
      }

      // else i==1 (return the target)
      if (cmp_endpoints(cv) == SMALLER) return (ctr_max_v(cv));

      return (ctr_min_v(cv));
    }
  };

  Construct_vertex_2 construct_vertex_2_object() const
  { return Construct_vertex_2(*this); }

  class Orientation_2 {
  protected:
    typedef Gps_traits_adaptor<Base>        Traits_adaptor;

    /*! The traits (in case it has state) */
    const Traits_adaptor& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Orientation_2(const Traits_adaptor& traits) : m_traits(traits) {}

    friend class Gps_traits_adaptor<Base>;

  public:
    template <class CurveInputIteraor>
    Orientation operator()(CurveInputIteraor begin,
                           CurveInputIteraor end) const
    {
      Compare_xy_2 cmp_xy = m_traits.compare_xy_2_object();
      Compare_y_at_x_right_2 cmp_y_at_x_right =
        m_traits.compare_y_at_x_right_2_object();
      Construct_vertex_2 ctr_v = m_traits.construct_vertex_2_object();

      CurveInputIteraor from_left_most = begin;
      CurveInputIteraor into_left_most = end;

      Point_2 left_most_v = ctr_v(*from_left_most, 0);

      --into_left_most;

      CurveInputIteraor ci = from_left_most;

      for (++ci ; ci != end; ++ci) {
        Comparison_result res_xy = cmp_xy( ctr_v(*ci, 0), left_most_v);
        if (res_xy == LARGER) continue;
        if (res_xy == SMALLER) {
          left_most_v =  ctr_v(*ci, 0);
          from_left_most = into_left_most = ci;
          --into_left_most;
        }
        else {
          // res_xy == EQUAL
          CurveInputIteraor tmp_from_left_most = ci;
          CurveInputIteraor tmp_into_left_most = ci;
          --tmp_into_left_most;

          Comparison_result res_from = cmp_y_at_x_right(*from_left_most,
                                                        *tmp_from_left_most,
                                                        left_most_v);

          Comparison_result res_to = cmp_y_at_x_right(*into_left_most,
                                                      *tmp_into_left_most,
                                                      left_most_v);

          CGAL_assertion(res_from != EQUAL && res_to != EQUAL);
          if(res_from == LARGER && res_to == SMALLER)
          {
            if(cmp_y_at_x_right(*tmp_from_left_most,
                                *into_left_most,
                                left_most_v) == LARGER)
            {
              from_left_most = tmp_from_left_most;
              into_left_most = tmp_into_left_most;
            }
          }
          else
            if (res_from == SMALLER && res_to == LARGER) {
              if (cmp_y_at_x_right(*tmp_into_left_most,
                                   *from_left_most,
                                   left_most_v) == LARGER)
              {
                from_left_most = tmp_from_left_most;
                into_left_most = tmp_into_left_most;
              }
            }
        }
      }// end for
      Comparison_result res = cmp_y_at_x_right(*into_left_most,
                                              *from_left_most,
                                              left_most_v);
      CGAL_assertion(res != EQUAL);
      if(res == SMALLER)
        return (CLOCKWISE);
      return (COUNTERCLOCKWISE);
    }
  };

  Orientation_2 orientation_2_object() const { return Orientation_2(*this); }
};

} //namespace CGAL

#endif
