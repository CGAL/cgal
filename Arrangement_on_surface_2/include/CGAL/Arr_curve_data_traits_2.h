// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CURVE_DATA_TRAITS_2_H
#define CGAL_ARR_CURVE_DATA_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_curve_data_traits_2<> class template.
 */

#include<CGAL/Object.h>
#include<CGAL/tags.h>
#include<CGAL/Arr_tags.h>
#include<CGAL/Arr_geometry_traits/Curve_data_aux.h>
#include<list>

namespace CGAL {

/*! \class
 * A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with an ordinary traits
 * class, which is also used as a based traits class to inherit from.
 * It can attach data objects to Curve_2 and to X_monotone_curve_2 objects
 * (possibly of two different types).
 * The data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split.
 * When two x-monotone curves overlap, the data field to be associated with
 * the overlapping subcurve is obtained from the merge functor.
 * All other functors are inherited from the base ordinary traits class.
 */
template <class Traits_, class XMonotoneCurveData_, 
          class Merge_ = _Default_merge_func<XMonotoneCurveData_>,
          class CurveData_ = XMonotoneCurveData_,
          class Convert_ = _Default_convert_func<CurveData_,
                                                 XMonotoneCurveData_> >
class Arr_curve_data_traits_2 : public Traits_ 
{
public:

  typedef Traits_                                     Base_traits_2;
  typedef XMonotoneCurveData_                         X_monotone_curve_data;
  typedef Merge_                                      Merge;
  typedef CurveData_                                  Curve_data;
  typedef Convert_                                    Convert;
  
  typedef typename Base_traits_2::Curve_2             Base_curve_2;
  typedef typename Base_traits_2::X_monotone_curve_2  Base_x_monotone_curve_2;
  typedef typename Base_traits_2::Point_2             Point_2;

  typedef typename Base_traits_2::Has_left_category   Has_left_category;
  typedef typename Base_traits_2::Has_merge_category  Base_has_merge_category;
  typedef Tag_true                                    Has_merge_category;
  typedef typename Base_traits_2::Has_do_intersect_category
                                                      Has_do_intersect_category;

  typedef typename internal::Arr_complete_left_side_category< Base_traits_2 >::Category
                                                      Left_side_category;
  typedef typename internal::Arr_complete_bottom_side_category< Base_traits_2 >::Category
                                                      Bottom_side_category;
  typedef typename internal::Arr_complete_top_side_category< Base_traits_2 >::Category
                                                      Top_side_category;
  typedef typename internal::Arr_complete_right_side_category< Base_traits_2 >::Category
                                                      Right_side_category;

  // Representation of a curve with an addtional data field:
  typedef _Curve_data_ex<Base_curve_2, Curve_data>    Curve_2;
  
  // Representation of an x-monotone curve with an addtional data field:
  typedef _Curve_data_ex<Base_x_monotone_curve_2,
                         X_monotone_curve_data>       X_monotone_curve_2;

  typedef typename Base_traits_2::Multiplicity        Multiplicity;
  
public:
  
  /// \name Construction.
  //@{

  /*! Default constructor. */
  Arr_curve_data_traits_2 ()
  {}
  
  /*! Constructor from a base-traits class. */
  Arr_curve_data_traits_2 (const Base_traits_2 & traits) :
    Base_traits_2 (traits)
  {}
  //@}

  /// \name Overriden functors.
  //@{

  class Make_x_monotone_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Make_x_monotone_2 (const Base_traits_2 * _base) :
      base (_base)
    {}
    
    /*!
     * Cut the given curve into x-monotone subcurves and insert them to the
     * given output iterator. As segments are always x_monotone, only one
     * x-monotone curve will be contained in the iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is X_monotone_curve_2.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    {
      // Make the original curve x-monotone.
      std::list<CGAL::Object>       base_objects;
    
      base->make_x_monotone_2_object() (cv,
                                        std::back_inserter (base_objects));

      // Attach the data to each of the resulting x-monotone curves.
      typename std::list<CGAL::Object>::const_iterator  it;
      const Base_x_monotone_curve_2  *base_x_curve;
      X_monotone_curve_data           xdata = Convert()(cv.data());

      for (it = base_objects.begin(); it != base_objects.end(); ++it)
      {
        base_x_curve = object_cast<Base_x_monotone_curve_2> (&(*it));
        if (base_x_curve != NULL)
        {
          // Current object is an x-monotone curve: Attach data to it.
          *oi = make_object (X_monotone_curve_2 (*base_x_curve,
                                                 xdata));
        }
        else
        {
          // Current object is an isolated point: Leave it as is.
          CGAL_assertion (object_cast<Point_2> (&(*it)) != NULL);
          *oi = *it;
        }
        ++oi;
      }

      return (oi);
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2 (this);
  }

  class Split_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Split_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      // Split the original curve.
      base->split_2_object() (cv, p, c1, c2);

      // Attach data to the split curves.
      c1.set_data (cv.data());
      c2.set_data (cv.data());

      return;
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2 (this);
  }

  class Intersect_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Intersect_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      // Use the base functor to obtain all intersection objects.
      std::list<CGAL::Object>                   base_objects;

      base->intersect_2_object() (cv1, cv2,
                                  std::back_inserter (base_objects));

      // Stop if the list is empty:
      if (base_objects.empty())
        return (oi);

      // Go over all intersection objects and prepare the output.
      typename std::list<CGAL::Object>::const_iterator  it;
      const Base_x_monotone_curve_2                    *base_cv;

      for (it = base_objects.begin(); it != base_objects.end(); ++it)
      {
        if ((base_cv = object_cast<Base_x_monotone_curve_2> (&(*it))) != NULL)
        {
          // The current intersection object is an overlapping x-monotone
          // curve: Merge the data fields of both intersecting curves and
          // associate the result with the overlapping curve.
          X_monotone_curve_2  cv (*base_cv,
                                  Merge() (cv1.data(), cv2.data()));

          *oi = make_object (cv);
        }
        else
        {
          // The current intersection object is an intersection point:
          // Copy it as is.
          *oi = *it;
        }
        ++oi;
      }

      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2 (this);
  }

  class Are_mergeable_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Are_mergeable_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      return (_are_mergeable_base_imp (cv1, cv2, Base_has_merge_category()));
    }

  private:

    /*!
     * Implementation of the base predicate in case the HasMerge tag is true.
     */
    bool _are_mergeable_base_imp (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  Tag_true) const
    {
      // In case the two base curves are not mergeable, the extended curves
      // are not mergeable as well.
      if (! (base->are_mergeable_2_object() (cv1, cv2)))
        return (false);

      // In case the two base curves are mergeable, check that they have the
      // same data fields.
      return (cv1.data() == cv2.data());
    }

    /*!
     * Implementation of the base predicate in case the HasMerge tag is false.
     */
    bool _are_mergeable_base_imp (const X_monotone_curve_2& ,
                                  const X_monotone_curve_2& ,
                                  Tag_false) const
    {
      // Curve merging is not supported:
      return (false);
    }
  };
  
  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2 (this);
  }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Merge_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      // The function is implemented based on the base Has_merge category.
      _merge_imp (cv1, cv2, c, Base_has_merge_category());
    }

  private:

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    void _merge_imp (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c,
                     Tag_true) const
    {      
      // Merge the two base curve.
      Base_x_monotone_curve_2  base_cv;

      base->merge_2_object() (cv1, cv2, base_cv);

      // Attach data from one of the curves.
      CGAL_precondition (cv1.data() == cv2.data());

      c = X_monotone_curve_2 (base_cv, cv1.data());
      return;
    }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    void _merge_imp (const X_monotone_curve_2& ,
                     const X_monotone_curve_2& ,
                     X_monotone_curve_2& ,
                     Tag_false) const
    {
      // This function should never be called!
      CGAL_error_msg("Merging curves is not supported.");
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2 (this);
  }

  class Construct_x_monotone_curve_2
  {
  private:
    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Construct_x_monotone_curve_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return An x-monotone curve connecting p and q.
     */
    X_monotone_curve_2 operator() (const Point_2& p, const Point_2& q) const
    {
      Base_x_monotone_curve_2  base_cv =
        base->construct_x_monotone_curve_2_object() (p, q);

      return (X_monotone_curve_2 (base_cv, X_monotone_curve_data()));
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2 (this);
  }
  //@}

};

} //namespace CGAL

#endif
