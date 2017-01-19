// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_ARR_LABELED_TRAITS_2_H
#define CGAL_ARR_LABELED_TRAITS_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Labels.h>
#include <list>

namespace CGAL {

/*! \class
 * A meta-traits class that adds lables to points and to x-monotone curves,
 * such that the comparison of two points, as well as the computation of the
 * intersections between two segments can be easily filtered.
 */
template <class Traits_>
class Arr_labeled_traits_2 : public Traits_
{
private:

  typedef Traits_                                      Base_traits_2;
  typedef typename Base_traits_2::Point_2              Base_point_2;
  typedef typename Base_traits_2::X_monotone_curve_2   Base_x_monotone_curve_2;

public:

  /*! \class
   * A point extended by a label.
   */
  class Point_2 : public Base_point_2
  {
  private:

    Point_label         _label;

  public:

    /*! Default constructor. */
    Point_2 ()
    {}

    /*! Constructor from a base point. */
    Point_2 (const Base_point_2& p) :
      Base_point_2 (p),
      _label()
    {}

    /*! Constructor from a point an a label. */
    Point_2 (const Base_point_2& p, const Point_label& label) :
      Base_point_2 (p),
      _label (label)
    {}

    /*! Get the label. */
    const Point_label& label () const
    {
      return (_label);
    }
  };

  /*! \class
   * An x-monotone curve extended by a label.
   */
  class X_monotone_curve_2 : public Base_x_monotone_curve_2
  {
  private:

    X_curve_label         _label;

  public:

    /*! Default constructor. */
    X_monotone_curve_2 ()
    {}

    /*! Constructor from a base x-monotone curve. */
    X_monotone_curve_2 (const Base_x_monotone_curve_2& p) :
      Base_x_monotone_curve_2 (p),
      _label()
    {}

    /*! Constructor from an x-monotone curve an a label. */
    X_monotone_curve_2 (const Base_x_monotone_curve_2& p,
                        const X_curve_label& label) :
      Base_x_monotone_curve_2 (p),
      _label (label)
    {}

    /*! Get the label (const version). */
    const X_curve_label& label () const
    {
      return (_label);
    }

    /*! Get the label (non-const version). */
    X_curve_label& label ()
    {
      return (_label);
    }

    /*! Set the label. */
    void set_label (const X_curve_label& label)
    {
      _label = label;
      return;
    }
  };

  typedef typename Base_traits_2::Has_left_category      Has_left_category;
  typedef Tag_false                                      Has_merge_category;

  /*! Default constructor. */
  Arr_labeled_traits_2 ()
  {}

  // Inherited functors:
  typedef typename Base_traits_2::Is_vertical_2         Is_vertical_2;
  typedef typename Base_traits_2::Compare_y_at_x_2      Compare_y_at_x_2;
  typedef typename Base_traits_2::Compare_y_at_x_right_2
                                                        Compare_y_at_x_right_2;
  typedef typename Base_traits_2::Equal_2               Equal_2;

  /// \name Overriden functors.
  //@{
  class Compare_x_2
  {
  private:

    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Compare_x_2 (const Base_traits_2 * _base) :
      base (_base)
    {}

    /*!
     * Compare the x-coordinates of two points.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      // If two points have the same label, they are equal.
      if (p1.label() == p2.label())
        return (EQUAL);

      return (base->compare_x_2_object()(p1, p2));
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (this));
  }


  class Compare_xy_2
  {
  private:

    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Compare_xy_2 (const Base_traits_2 *_base) :
      base (_base)
    {}

    /*!
     * Compare two points lexigoraphically: by x, then by y.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      // If two points have the same label, they are equal.
      if (p1.label() == p2.label())
        return (EQUAL);

      return (base->compare_xy_2_object()(p1, p2));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return (Compare_xy_2 (this));
  }


  class Construct_min_vertex_2
  {
  private:

    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Construct_min_vertex_2 (const Base_traits_2 *_base) :
      base (_base)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve.
     */
    Point_2 operator() (const X_monotone_curve_2& cv) const
    {
      const Base_point_2&  pt = base->construct_min_vertex_2_object() (cv);

      if (cv.label().right_count() == 1 && cv.label().left_count() == 0)
      {
        // A curve directed from left to right:
        Point_label   label (cv.label().component(), cv.label().index());

        return (Point_2 (pt, label));
      }
      else if (cv.label().right_count() == 0 && cv.label().left_count() == 1)
      {
        // A curve directed from right to left:
        Point_label   label (cv.label().component(),
                             cv.label().is_last() ? 0 : cv.label().index()+1);

        return (Point_2 (pt, label));
      }

      // Assign an invalid label to the point.
      return (Point_2 (pt));
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return (Construct_min_vertex_2 (this));
  }


  class Construct_max_vertex_2
  {
  private:

    const Base_traits_2 * base;

  public:

    /*! Constructor. */
    Construct_max_vertex_2 (const Base_traits_2 *_base) :
      base (_base)
    {}

    /*!
     * Get the right endpoint of the x-monotone curve.
     */
    Point_2 operator() (const X_monotone_curve_2& cv) const
    {
      const Base_point_2&  pt = base->construct_max_vertex_2_object() (cv);

      if (cv.label().right_count() == 1 && cv.label().left_count() == 0)
      {
        // A curve directed from left to right:
        Point_label   label (cv.label().component(),
                             cv.label().is_last() ? 0 : cv.label().index()+1);

        return (Point_2 (pt, label));
      }
      else if (cv.label().right_count() == 0 && cv.label().left_count() == 1)
      {
        // A curve directed from right to left:
        Point_label   label (cv.label().component(), cv.label().index());

        return (Point_2 (pt, label));
      }

      // Assign an invalid label to the point.
      return (Point_2 (pt));
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return (Construct_max_vertex_2 (this));
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
     */
    void operator() (const X_monotone_curve_2& cv, const Point_2& p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      // Split the base curve into two.
      base->split_2_object() (cv, p, c1, c2);

      // Duplicate the label to both subcurves.
      c1.set_label (cv.label());
      c2.set_label (cv.label());

      return;
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return (Split_2 (this));
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
     * given output iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      // In case the curves are adjacent in their curve sequence, we do
      // not have to compute their intersection (we already know that they
      // have just one common endpoint).
      if (cv1.label().is_adjacent (cv2.label()))
        return (oi);

      // Compute the intersection.
      std::list<CGAL::Object>            base_objs;

      base->intersect_2_object() (cv1, cv2, std::back_inserter (base_objs));

      if (base_objs.empty())
        return (oi);

      // Attach labels to the intersection objects.
      std::list<CGAL::Object>::iterator             obj_it;
      const std::pair<Base_point_2, unsigned int>  *base_pt;
      const Base_x_monotone_curve_2                *base_xcv;

      for (obj_it = base_objs.begin(); obj_it != base_objs.end(); ++obj_it)
      {
        base_pt =
          object_cast<std::pair<Base_point_2, unsigned int> > (&(*obj_it));

        if (base_pt != NULL)
        {
          // Attach an invalid label to an itersection point.
          *oi = CGAL::make_object
            (std::make_pair (Point_2 (base_pt->first), base_pt->second));
          ++oi;
        }
        else
        {
          base_xcv = object_cast<Base_x_monotone_curve_2> (&(*obj_it));
          CGAL_assertion (base_xcv != NULL);

          // Attach a merged label to the overlapping curve.
          *oi = CGAL::make_object
            (X_monotone_curve_2 (*base_xcv,
                                 X_curve_label (cv1.label(), cv2.label())));
          ++oi;
        }
      }

      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (this));
  }
  //@}
};

} //namespace CGAL

#endif
