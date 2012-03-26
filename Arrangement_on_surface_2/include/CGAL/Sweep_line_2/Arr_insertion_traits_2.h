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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_INSERTION_TRAITS_2_H
#define CGAL_ARR_INSERTION_TRAITS_2_H

/*!
 * Defintion of the Arr_insertion_traits_2<Traits,Arrangement> class.
 */

#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>

namespace CGAL {

/*! \class
 * A meta-traits class that stores a halfedge handle with every x-monotone
 * curve, and a vertex handle with each point. This information is used to
 * speed up the aggregated insertion process.
 */
template <typename Traits_, typename Arrangement_>
class Arr_insertion_traits_2 : 
  public Arr_basic_insertion_traits_2<Traits_, Arrangement_>
{
public:

  typedef Traits_                                               Traits_2;
  typedef Arr_basic_insertion_traits_2<Traits_, Arrangement_>   Base;

  typedef typename Traits_2::Intersect_2              Base_intersect_2;
  typedef typename Traits_2::Split_2                  Base_split_2;
  typedef typename Base::Base_x_monotone_curve_2      Base_x_monotone_curve_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Base::Halfedge_handle              Halfedge_handle;
  typedef typename Base::Base_point_2                 Base_point_2;
  typedef typename Base::Point_2                      Point_2;
  typedef typename Base::Multiplicity                 Multiplicity;

  typedef typename Base::Has_left_category            Has_left_category;
  typedef typename Base::Has_do_intersect_category    Has_do_intersect_category;

  // should be ok, as basic_insertion (=Base) completes incomplete tags
  typedef typename Base::Left_side_category       Left_side_category;
  typedef typename Base::Bottom_side_category     Bottom_side_category;
  typedef typename Base::Top_side_category        Top_side_category;
  typedef typename Base::Right_side_category      Right_side_category;

  /* Insertion is implemented as sweep-line visitor. The sweep-line algorithm
   * never performs merging of curves. Therefore, AreMergeable_2 and
   * Merge_2 are not needed either.
   */
  typedef Tag_false                                   Has_merge_category;
  
public:

  /*! Constructor with a traits class. */
  Arr_insertion_traits_2 (const Traits_2& tr) :
    Base (tr)
  {}

  /*! A functor that compares compares the y-coordinates of two x-monotone
   * curves immediately to the right of their intersection point.
   */
  class Intersect_2 {
  protected:
    //! The base operators.
    Base_intersect_2      m_base_intersect;
    Halfedge_handle       invalid_he;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Intersect_2 (const Base_intersect_2& base) :
      m_base_intersect (base),
      invalid_he()
    {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_insertion_traits_2<Traits_2, Arrangement_>;
    
  public:
    template<typename OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi)
    {
      if(cv1.halfedge_handle() != invalid_he &&
         cv2.halfedge_handle() != invalid_he)
      {
        // The curves are interior-disjoint as both of them are already in
        //  the arrangement.
        return (oi);
      }

      OutputIterator           oi_end = m_base_intersect(cv1.base(),
                                                         cv2.base(), oi);
      const Base_x_monotone_curve_2                *base_overlap_cv;
      const std::pair<Base_point_2, unsigned int>  *intersect_p;

      // convert objects that are associated with Base_x_monotone_curve_2 to
      // X_monotone_curve_2 
      for(; oi != oi_end; ++oi)
      {
        base_overlap_cv = object_cast<Base_x_monotone_curve_2> (&(*oi));
        if (base_overlap_cv != NULL)
        {
          // Add halfedge handles to the resulting curve.
          Halfedge_handle  he;

          if (cv1.halfedge_handle() != invalid_he)
            he = cv1.halfedge_handle();
          else if (cv2.halfedge_handle() != invalid_he)
            he = cv2.halfedge_handle();

          X_monotone_curve_2    overlap_cv (*base_overlap_cv, he);

          overlap_cv.set_overlapping();
          *oi = make_object (overlap_cv);
        }
        else
        {
          intersect_p = 
            object_cast<std::pair<Base_point_2, unsigned int> > (&(*oi));

          CGAL_assertion (intersect_p != NULL);

          *oi = make_object (std::make_pair (Point_2(intersect_p->first),
                                             intersect_p->second));
        }
      }

      // Return a past-the-end iterator.
      return oi_end;
    }
  };

  /*! Obtain a Intersect_2 function object */
  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (this->m_base_traits->intersect_2_object())); 
  }

  /*! A functor that splits an arc at a point. */
  class Split_2 {
  protected:
    //! The base operator.
    Base_split_2    m_base_split;

    /*! Constructor.
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Split_2(const Base_split_2& base) : m_base_split (base) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_insertion_traits_2<Traits_2, Arrangement_>;
    
  public:
    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      m_base_split(cv.base(), p.base(), c1.base(), c2.base());
      c1.set_halfedge_handle(cv.halfedge_handle());
      c2.set_halfedge_handle(cv.halfedge_handle());
    }
  };

  /*! Obtain a plit_2 function object */
  Split_2 split_2_object () const
  {
    return (Split_2 (this->m_base_traits->split_2_object()));
  }
};

} //namespace CGAL

#endif
