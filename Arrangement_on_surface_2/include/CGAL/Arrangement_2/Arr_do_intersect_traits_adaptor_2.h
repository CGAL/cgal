// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@gmail.com>
//            Waqar Khan        <wkhan@mpi-inf.mpg.de>
//            Simon Giraudot    <simon.giraudot@geometryfactory.com>

#ifndef CGAL_ARR_DO_INTERSECT_TRAITS_ADAPTOR_2_H
#define CGAL_ARR_DO_INTERSECT_TRAITS_ADAPTOR_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <boost/iterator/function_output_iterator.hpp>

namespace CGAL
{

template <typename TraitsBase>
class Arr_do_intersect_traits_adaptor_2 : public TraitsBase
{
public:
  typedef TraitsBase Base;
  typedef typename Base::Point_2                      Point_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;

  class Exception { };

  class Intersect_2 {
  protected:
    typedef Arr_do_intersect_traits_adaptor_2<TraitsBase>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Construct
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend Traits;
  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may intersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      return intersect_impl (cv1, cv2, oi, typename Base::Has_do_intersect_category());
    }

  protected:
    template <typename OutputIterator>
    OutputIterator intersect_impl (const X_monotone_curve_2& cv1,
                                   const X_monotone_curve_2& cv2,
                                   OutputIterator oi,
                                   const Tag_true&) const
    {
      typename Base::Do_intersect_2
        do_intersect_2 = dynamic_cast<const Base*>(&m_traits)->do_intersect_2_object();
      if (do_intersect_2(cv1, cv2))
        throw Exception();
      return oi;
    }

    template <typename OutputIterator>
    OutputIterator intersect_impl (const X_monotone_curve_2& cv1,
                                   const X_monotone_curve_2& cv2,
                                   OutputIterator oi,
                                   const Tag_false&) const
    {
      typename Base::Intersect_2
        intersect_2 = dynamic_cast<const Base*>(&m_traits)->intersect_2_object();
      intersect_2(cv1, cv2, boost::make_function_output_iterator([](const auto&) { throw Exception(); }));
      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }
};

} // namespace CGAL

#endif // CGAL_ARR_DO_INTERSECT_TRAITS_ADAPTOR_2_H
