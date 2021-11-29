// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Efi Fogel       <efifogel@gmail.com>

#ifndef CGAL_GPS_SIMPLIFIER_TRAITS_H
#define CGAL_GPS_SIMPLIFIER_TRAITS_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_traits_decorator.h>

namespace CGAL {

class Gps_simplifier_curve_data {
protected:
  unsigned int m_bc;
  unsigned int m_twin_bc;
  unsigned int m_index;

public:
  Gps_simplifier_curve_data() {}

  Gps_simplifier_curve_data(unsigned int bc, unsigned int twin_bc,
                            unsigned int index):
    m_bc(bc),
    m_twin_bc(twin_bc),
    m_index(index)
  {}

  unsigned int bc() const { return m_bc; }

  unsigned int twin_bc() const { return m_twin_bc; }

  unsigned int index() const { return m_index; }

  unsigned int& index() { return m_index; }

  unsigned int& twin_bc() { return m_twin_bc; }

  void set_bc(unsigned int bc) { m_bc = bc; }

  void set_twin_bc(unsigned int twin_bc) { m_twin_bc = twin_bc; }

  void set_index(unsigned int index) { m_index = index; }
};

struct Gps_simplifier_point_data {
protected:
  unsigned int m_index;

public:
  Gps_simplifier_point_data() {}

  Gps_simplifier_point_data(unsigned int index) : m_index(index) {}

  unsigned int index() const { return m_index; }

  void set_index(unsigned int index) { m_index = index; }
};

template <typename Traits_>
class Gps_simplifier_traits :
  public Gps_traits_decorator<Traits_,
                              Gps_simplifier_curve_data,
                              Gps_simplifier_point_data>
{
public:
  typedef Traits_                                           Traits;
  typedef Gps_traits_decorator<Traits_,
                               Gps_simplifier_curve_data,
                               Gps_simplifier_point_data>   Base;
  typedef Gps_simplifier_traits<Traits>                     Self;
  typedef typename Traits::X_monotone_curve_2     Base_x_monotone_curve_2;
  typedef typename Traits::Point_2                Base_point_2;
  typedef typename Traits::Construct_min_vertex_2 Base_Construct_min_vertex_2;
  typedef typename Traits::Construct_max_vertex_2 Base_Construct_max_vertex_2;
  typedef typename Traits::Compare_endpoints_xy_2 Base_Compare_endpoints_xy_2;
  typedef typename Traits::Compare_xy_2           Base_Compare_xy_2;
  typedef typename Traits::Compare_y_at_x_right_2 Base_Compare_y_at_x_right_2;
  typedef typename Traits::Compare_y_at_x_2       Base_Compare_y_at_x_2;
  typedef typename Traits::Intersect_2            Base_Intersect_2;
  typedef typename Traits::Split_2                Base_Split_2;

protected:
  mutable unsigned int m_pgn_size;

public:
  typedef typename Base::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename Base::Point_2                  Point_2;
  typedef typename Base::Multiplicity             Multiplicity;

  typedef typename Base::Curve_data               Curve_data;
  typedef typename Base::Point_data               Point_data;

  Gps_simplifier_traits() {}

  Gps_simplifier_traits(const Traits& tr) : Base(tr) {}

  unsigned int polygon_size() const { return m_pgn_size; }

  void set_polygon_size(unsigned int pgn_size) const { m_pgn_size = pgn_size; }

  bool is_valid_index(unsigned int index) const
  { return (index < m_pgn_size); }

  unsigned int invalid_index() const { return (m_pgn_size); }

  class Intersect_2 {
  private:
    /*! The traits (in case it has state) */
    const Self& m_traits;

    /*! Constructor. */
    Intersect_2(const Self& tr) : m_traits(tr) {}

    friend Self;

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef const std::pair<Base_point_2, Multiplicity>
        Intersection_base_point;
      typedef boost::variant<Intersection_base_point, Base_x_monotone_curve_2>
                                                        Intersection_base_result;
      typedef const std::pair<Point_2, Multiplicity>    Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      const auto* base_traits = m_traits.m_base_traits;
      auto base_cmp_xy = base_traits->compare_xy_2_object();
      auto base_cmp_endpoints = base_traits->compare_endpoints_xy_2_object();
      auto base_ctr_min_vertex = base_traits->construct_min_vertex_2_object();
      auto base_intersect = base_traits->intersect_2_object();

      //// if the two curves are incident, do not intersect them
      //if (m_traits.is_valid_index(cv1.data().index()) &&
      //    m_traits.is_valid_index(cv2.data().index()))
      //{
      //  unsigned int index_diff =
      //    (cv1.data().index() > cv2.data().index()) ?
      //    (cv1.data().index() - cv2.data().index()):
      //    (cv2.data().index() - cv1.data().index());

      //  if(index_diff == 1 ||index_diff == m_traits.polygon_size() -1)
      //  {
      //    return (oi);
      //  }
      //}
      std::vector<Intersection_base_result> xections;
      if (base_cmp_xy(base_ctr_min_vertex(cv1.base()),
                      base_ctr_min_vertex(cv2.base())) == LARGER)
        base_intersect(cv1.base(), cv2.base(), back_inserter(xections));
      else
        base_intersect(cv2.base(), cv1.base(), back_inserter(xections));

      // convert objects that are associated with Base_x_monotone_curve_2 to
      // the extenede X_monotone_curve_2
      for (const auto& xection : xections) {
        const Intersection_base_point* base_pt =
          boost::get<Intersection_base_point>(&xection);
        if (base_pt != nullptr) {
          Point_data pt_data(m_traits.invalid_index());
          Point_2 point_plus(base_pt->first, pt_data); // the extended point
          *oi++ =
            Intersection_result(std::make_pair(point_plus, base_pt->second));
          continue;
        }

        const Base_x_monotone_curve_2* overlap_cv =
          boost::get<Base_x_monotone_curve_2>(&xection);

        CGAL_assertion(overlap_cv != nullptr);
        unsigned int ov_bc;
        unsigned int ov_twin_bc;
        if (base_cmp_endpoints(cv1) == base_cmp_endpoints(cv2)) {
          // cv1 and cv2 have the same directions
          ov_bc = cv1.data().bc() + cv2.data().bc();
          ov_twin_bc = cv1.data().twin_bc() + cv2.data().twin_bc();
        }
        else {
          // cv1 and cv2 have opposite directions
          ov_bc = cv1.data().bc() + cv2.data().twin_bc();
          ov_twin_bc = cv1.data().twin_bc() + cv2.data().bc();
        }

        if (base_cmp_endpoints(*overlap_cv) != base_cmp_endpoints(cv1)) {
          // overlap_cv, cv1 have opposite directions
          std::swap(ov_bc, ov_twin_bc);
        }

        Curve_data cv_data(ov_bc, ov_twin_bc, m_traits.invalid_index());
        *oi++ = Intersection_result(X_monotone_curve_2(*overlap_cv, cv_data));
      }

      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return Intersect_2(*this); }

  class Split_2 {
  private:
    const Self& m_traits;

    /*! Constructor. */
    Split_2(const Self& tr) : m_traits(tr) {}

    friend Self;

  public:
    void operator()(const X_monotone_curve_2& cv, const Point_2 & p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      const auto* base_traits = m_traits.m_base_traits;
      auto base_split = base_traits->split_2_object();
      base_split(cv.base(), p.base(), c1.base(), c2.base());
      const Curve_data& cv_data = cv.data();
      c1.set_data(Curve_data(cv_data.bc(), cv_data.twin_bc(),
                             m_traits.invalid_index()));

      c2.set_data(Curve_data(cv_data.bc(), cv_data.twin_bc(),
                             m_traits.invalid_index()));
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const { return Split_2(*this); }

  class Construct_min_vertex_2 {
  private:
    const Self& m_traits;

    Construct_min_vertex_2(const Self& tr) : m_traits(tr) {}

    friend Self;

  public:
    /*! Obtain the left endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      const auto* base_traits = m_traits.m_base_traits;
      auto base_ctr_min_vertex = base_traits->construct_min_vertex_2_object();

      if (! m_traits.is_valid_index(cv.data().index()))
        return Point_2(base_ctr_min_vertex(cv.base()), m_traits.invalid_index());

      auto base_cmp_endpoints = base_traits->compare_endpoints_xy_2_object();
      Comparison_result res = base_cmp_endpoints(cv);
      Point_data pt_data;
      if (res == SMALLER) {
        // min vertex is the source
        pt_data.set_index(cv.data().index());
      }
      else {
        // min vertex is the target
        pt_data.set_index((cv.data().index() + 1) % m_traits.polygon_size());
      }
      return Point_2(base_ctr_min_vertex(cv.base()), pt_data);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  { return Construct_min_vertex_2(*this); }

  class Construct_max_vertex_2 {
  private:
    const Self& m_traits;

    Construct_max_vertex_2(const Self& tr) : m_traits(tr) {}

    friend Self;

  public:
    /*! Obtain the right endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      const auto* base_traits = m_traits.m_base_traits;
      auto base_ctr_max_vertex = base_traits->construct_max_vertex_2_object();
      if (! m_traits.is_valid_index(cv.data().index()))
        return Point_2(base_ctr_max_vertex(cv.base()), m_traits.invalid_index());

      auto base_cmp_endpoints = base_traits->compare_endpoints_xy_2_object();
      Comparison_result res = base_cmp_endpoints(cv);
      Point_data pt_data;
      if (res == SMALLER) {
        // min vertex is the target
        pt_data.set_index((cv.data().index() + 1) % m_traits.polygon_size());
      }
      else {
        // min vertex is the source
        pt_data.set_index(cv.data().index());
      }
      return Point_2(base_ctr_max_vertex(cv.base()), pt_data);
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  { return Construct_max_vertex_2(*this); }

  class Compare_xy_2 {
  private:
    const Self& m_traits;

    Compare_xy_2(const Self& tr) : m_traits(tr) {}

    friend Self;

  public:
    /*! Obtain the left endpoint of the x-monotone curve (segment).
      * \param cv The curve.
      * \return The left endpoint.
      */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const auto* base_traits = m_traits.m_base_traits;
      auto base_cmp_xy = base_traits->compare_xy_2_object();

      //if one of the indexes is invalid, compare p1 and p2
      if (! m_traits.is_valid_index(p1.data().index()) ||
          ! m_traits.is_valid_index(p2.data().index()))
        return (base_cmp_xy(p1.base(), p2.base()));

      // if the two point has the same index, return EQUAL
      if (p1.data().index() == p2.data().index()) return EQUAL;

      return (base_cmp_xy(p1.base(), p2.base()));
    }
  };


  /*! Get a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const { return Compare_xy_2(*this); }
};

} //namespace CGAL

#endif
