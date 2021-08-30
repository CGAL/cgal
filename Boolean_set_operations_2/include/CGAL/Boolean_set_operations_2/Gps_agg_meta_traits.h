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

#ifndef CGAL_BSO_2_GPS_AGG_META_TRAITS_H
#define CGAL_BSO_2_GPS_AGG_META_TRAITS_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <vector>
#include <boost/mpl/assert.hpp>
#include <CGAL/Boolean_set_operations_2/Gps_traits_decorator.h>
#include <CGAL/Boolean_set_operations_2/Curve_with_halfedge.h>
#include <CGAL/Boolean_set_operations_2/Point_with_vertex.h>

namespace CGAL {

template <typename Arrangement_>
class Gps_agg_curve_data : public Curve_with_halfedge<Arrangement_>
{
protected:
  typedef Arrangement_                             Arrangement;
  typedef typename Arrangement::Halfedge_handle    Halfedge_handle;
  typedef Curve_with_halfedge<Arrangement_>        Base;

  const Arrangement* m_arr; // pointer to the arrangement containing the edge.
  unsigned int m_bc;        // the boudary counter of the halfedge with the same
                            // direction as the curve

  unsigned int m_twin_bc;   // the boudary counter of the halfedge with the same
                            // direction as the curve

public:
  Gps_agg_curve_data() :
    Base(),
    m_arr(nullptr),
    m_bc(0),
    m_twin_bc(0)
  {}

  Gps_agg_curve_data(const Arrangement* arr, Halfedge_handle he,
                     unsigned int bc, unsigned int twin_bc) :
    Base(he),
    m_arr(arr),
    m_bc(bc),
    m_twin_bc(twin_bc)
  {}

  unsigned int bc() const { return m_bc; }

  unsigned int twin_bc() const { return m_twin_bc; }

  unsigned int& bc() { return m_bc; }

  unsigned int& twin_bc() { return m_twin_bc; }

  void set_bc(unsigned int bc) { m_bc = bc; }

  void set_twin_bc(unsigned int twin_bc) { m_twin_bc = twin_bc; }

  const Arrangement* arr() const { return m_arr; }
};

template <typename Arrangement_>
class Gps_agg_meta_traits :
  public Gps_traits_decorator<typename Arrangement_::Traits_adaptor_2,
                              Gps_agg_curve_data<Arrangement_>,
                              Point_with_vertex<Arrangement_> >
{
  typedef Arrangement_                          Arrangement;
  typedef Arrangement                           Arr;

  typedef typename Arr::Traits_adaptor_2        Traits;
  typedef Traits                                Gt2;

  typedef typename Gt2::X_monotone_curve_2      Base_x_monotone_curve_2;
  typedef typename Gt2::Point_2                 Base_point_2;
  typedef typename Gt2::Construct_min_vertex_2  Base_Construct_min_vertex_2;
  typedef typename Gt2::Construct_max_vertex_2  Base_Construct_max_vertex_2;
  typedef typename Gt2::Compare_endpoints_xy_2  Base_Compare_endpoints_xy_2;
  typedef typename Gt2::Compare_xy_2            Base_Compare_xy_2;
  typedef typename Gt2::Compare_y_at_x_right_2  Base_Compare_y_at_x_right_2;
  typedef typename Gt2::Compare_y_at_x_2        Base_Compare_y_at_x_2;
  typedef typename Gt2::Intersect_2             Base_Intersect_2;
  typedef typename Gt2::Split_2                 Base_Split_2;

  typedef typename Gt2::Parameter_space_in_x_2  Base_Parameter_space_in_x_2;
  typedef typename Gt2::Compare_y_near_boundary_2
                                                Base_Compare_y_near_boundary_2;

  typedef typename Gt2::Parameter_space_in_y_2  Base_Parameter_space_in_y_2;
  typedef typename Gt2::Compare_x_near_boundary_2
                                                Base_Compare_x_near_boundary_2;

public:
  typedef typename Gt2::Multiplicity            Multiplicity;
  typedef Gps_agg_curve_data<Arr>               Curve_data;
  typedef Point_with_vertex<Arr>                Point_data;

private:
  typedef Gps_agg_meta_traits<Arrangement>                      Self;
  typedef Gps_traits_decorator<Gt2, Curve_data, Point_data>     Base;

public:
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Point_2                Point_2;
  typedef typename Gt2::Has_left_category       Has_left_category;
  typedef typename Gt2::Has_merge_category      Has_merge_category;
  typedef typename Gt2::Has_do_intersect_category
    Has_do_intersect_category;

  typedef typename Arr::Left_side_category      Left_side_category;
  typedef typename Arr::Bottom_side_category    Bottom_side_category;
  typedef typename Arr::Top_side_category       Top_side_category;
  typedef typename Arr::Right_side_category     Right_side_category;

  // a side is either oblivious or open (unbounded)
  BOOST_MPL_ASSERT((boost::mpl::or_<
                    boost::is_same<Left_side_category, Arr_oblivious_side_tag>,
                    boost::is_same<Left_side_category, Arr_open_side_tag> >));
  BOOST_MPL_ASSERT((boost::mpl::or_<
                    boost::is_same<Bottom_side_category, Arr_oblivious_side_tag>,
                    boost::is_same<Bottom_side_category, Arr_open_side_tag> >));
  BOOST_MPL_ASSERT((boost::mpl::or_<
                    boost::is_same<Top_side_category, Arr_oblivious_side_tag>,
                    boost::is_same<Top_side_category, Arr_open_side_tag> >));
  BOOST_MPL_ASSERT((boost::mpl::or_<
                    boost::is_same<Right_side_category, Arr_oblivious_side_tag>,
                    boost::is_same<Right_side_category, Arr_open_side_tag> >));

  typedef typename Arr::Halfedge_handle         Halfedge_handle;
  typedef typename Arr::Vertex_handle           Vertex_handle;

  Gps_agg_meta_traits() {}

  Gps_agg_meta_traits(const Gt2& base_tr) : Base(base_tr) {}

  class Intersect_2 {
  private:
    const Self& m_traits;

    /*! Constructor. */
    Intersect_2(const Self& traits) : m_traits(traits) {}

    friend Self;

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      // Check whether the curves are already in the same arrangement, and thus
      // must be interior-disjoint
      if (cv1.data().arr() == cv2.data().arr()) return oi;

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
          Point_2 point_plus(base_pt->first); // the extended point
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

        Curve_data cv_data(cv1.data().arr(), Halfedge_handle(),
                           ov_bc, ov_twin_bc);
        *oi++ = Intersection_result(X_monotone_curve_2(*overlap_cv, cv_data));
      }

      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  class Split_2 {
  private:
    Base_Split_2 m_base_split;

  public:
    /*! Construct. */
    Split_2(const Base_Split_2& base) : m_base_split(base) {}

    void operator()(const X_monotone_curve_2& cv, const Point_2 & p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      m_base_split(cv.base(), p.base(), c1.base(), c2.base());
      const Curve_data& cv_data = cv.data();
      c1.set_data(Curve_data(cv_data.arr(), Halfedge_handle(), cv_data.bc(),
                             cv_data.twin_bc()));

      c2.set_data(Curve_data(cv_data.arr(), Halfedge_handle(), cv_data.bc(),
                             cv_data.twin_bc()));
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const
  { return Split_2(this->m_base_traits->split_2_object()); }

  class Construct_min_vertex_2 {
  private:
    Base_Construct_min_vertex_2 m_base;

  public:
    Construct_min_vertex_2(const Base_Construct_min_vertex_2& base) :
      m_base(base)
    {}

    /*! Obtain the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      if (cv.data().halfedge() == Halfedge_handle())
        return (Point_2(m_base(cv.base())));

      CGAL_assertion
        ((Arr_halfedge_direction)cv.data().halfedge()->direction() ==
         ARR_LEFT_TO_RIGHT);
      return Point_2(m_base(cv.base()), cv.data().halfedge()->source());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2(this->m_base_traits->
                                  construct_min_vertex_2_object());
  }


  class Construct_max_vertex_2 {
  private:
    Base_Construct_max_vertex_2 m_base;

  public:
    Construct_max_vertex_2(const Base_Construct_max_vertex_2& base) :
        m_base(base)
    {}

    /*! Obtain the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    Point_2 operator()(const X_monotone_curve_2& cv) const
    {
      if (cv.data().halfedge() == Halfedge_handle())
        return (Point_2(m_base(cv.base())));

      CGAL_assertion((Arr_halfedge_direction)cv.data().halfedge()->direction() ==
                     ARR_LEFT_TO_RIGHT);
      return Point_2(m_base(cv.base()), cv.data().halfedge()->target());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return Construct_max_vertex_2(this->m_base_traits->
                                  construct_max_vertex_2_object());
  }

  class Compare_xy_2 {
  private:
    Base_Compare_xy_2 m_base;

  public:
    Compare_xy_2(const Base_Compare_xy_2& base) : m_base(base) {}

    /*! Obtain the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Point_data& inf1 = p1.data();
      const Point_data& inf2 = p2.data();

      if (inf1.vertex() == Vertex_handle() || inf2.vertex() == Vertex_handle())
        return m_base(p1.base(), p2.base());

      if (inf1.vertex() == inf2.vertex()) return EQUAL;
      return m_base(p1.base(), p2.base());
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(this->m_base_traits->compare_xy_2_object()); }

  // left-right
  class Parameter_space_in_x_2 {
  private:
    Base_Parameter_space_in_x_2 m_base;

  public:
    Parameter_space_in_x_2(const Base_Parameter_space_in_x_2& base) :
      m_base(base)
    {}

    /*! Obtain the parameter space at the end of a curve-end along the x-axis.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & cv,
                                   const Arr_curve_end& end) const
    { return m_base(cv.base(), end); }

    /*! Obtain the parameter space for a curve along the x-axis.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & cv) const
    { return m_base(cv.base()); }

    /*! Obtain the parameter space for a point along the x-axis.
     */
    Arr_parameter_space operator()(const Point_2 & pt) const
    { return m_base(pt.base()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  {
    return Parameter_space_in_x_2(this->m_base_traits->
                                  parameter_space_in_x_2_object());
  }

  class Compare_y_near_boundary_2 {
  private:
    Base_Compare_y_near_boundary_2 m_base;

  public:
    Compare_y_near_boundary_2(const Base_Compare_y_near_boundary_2& base) :
      m_base(base)
    {}

    /*! Compare the relative y-positions of two curve ends.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    { return m_base(xcv1, xcv2, ce); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  {
    return Compare_y_near_boundary_2(this->m_base_traits->
                                     compare_y_near_boundary_2_object()
    );
  }

  // TODO Compare_y_on_boundary_2
  // TODO Is_on_x_identification_2

  // bottom-top

  class Parameter_space_in_y_2 {
  private:
    Base_Parameter_space_in_y_2 m_base;

  public:
    Parameter_space_in_y_2(const Base_Parameter_space_in_y_2& base) :
      m_base(base)
    {}

    /*! Obtain the parameter space at the end of a curve-end along the y-axis.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & cv,
                                   const Arr_curve_end& end) const
    { return m_base(cv.base(), end); }

    /*! Obtain the parameter space for a curve along the x-axis.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & cv) const
    { return m_base(cv.base()); }

    /*! Obtain the parameter space for a point along the x-axis.
     */
    Arr_parameter_space operator()(const Point_2 & pt) const
    { return m_base(pt.base()); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  {
    return Parameter_space_in_y_2(this->m_base_traits->
                                  parameter_space_in_y_2_object());
  }

  class Compare_x_near_boundary_2 {
  private:
    Base_Compare_x_near_boundary_2 m_base;

  public:
    Compare_x_near_boundary_2(const Base_Compare_x_near_boundary_2& base) :
      m_base(base)
    {}

    /*! Compare the relative x-positions of a vertical curve and another given
     * curve end.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { return m_base(p, xcv, ce); }

    /*! Compare the relative x-positions of two curve ends.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    { return m_base(xcv1, ce1, xcv2, ce2); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  {
    return Compare_x_near_boundary_2(this->m_base_traits->
                                     compare_x_near_boundary_2_object());
  }

  // TODO Compare_x_on_boundary_2
  // TODO Is_on_y_identification_2
};

} // namespace CGAL

#endif
