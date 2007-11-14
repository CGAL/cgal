// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_ARR_INSERTION_TRAITS_2_H
#define CGAL_ARR_INSERTION_TRAITS_2_H

/*!
 * Defintion of the Arr_insertion_traits_2<Traits,Arrangement> class.
 */

#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A meta-traits class that stores a halfedge handle with every x-monotone
 * curve, and a vertex handle with each point. This information is used to
 * speed up the aggregated insertion process.
 */
template <class Traits_, class Arrangement_>
class Arr_insertion_traits_2 : 
  public Arr_basic_insertion_traits_2<Traits_, Arrangement_>
{
public:

  typedef Traits_                                     Traits_2;
  typedef Arr_basic_insertion_traits_2<Traits_,
                                       Arrangement_>  Base;

  typedef typename Traits_2::Intersect_2              Base_intersect_2;
  typedef typename Traits_2::Split_2                  Base_split_2;
  typedef typename Base::Base_x_monotone_curve_2      Base_x_monotone_curve_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Base::Halfedge_handle              Halfedge_handle;
  typedef typename Base::Base_point_2                 Base_point_2;
  typedef typename Base::Point_2                      Point_2;

  typedef typename Base::Has_boundary_category        Has_boundary_category;
  typedef typename Base::Boundary_category            Boundary_category;
  typedef typename Base::Has_left_category            Has_left_category;
  typedef Tag_false                                   Has_merge_category;

public:

  /*! Constructor with a traits class. */
  Arr_insertion_traits_2 (Traits_2& tr) :
    Base (tr)
  {}

  /*! \class
   * The Intersect_2 functor.
   */
  class Intersect_2
  {
  private:

    Base_intersect_2      m_base_intersect;
    Halfedge_handle       invalid_he;

  public:
   
    /*! Constructor. */
    Intersect_2 (const Base_intersect_2& base) :
      m_base_intersect (base),
      invalid_he()
    {}

    template<class OutputIterator>
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

  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (this->m_base_traits->intersect_2_object())); 
  }


   /*! \class
   * The Split_2 functor.
   */
  class Split_2
  {
  private:
    Base_split_2    m_base_split;

  public:

    /*! Constructor. */
    Split_2 (const Base_split_2& base) :
        m_base_split (base)
    {}

    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2)
    {
      m_base_split(cv.base(),
                   p.base(),
                   c1.base(),
                   c2.base());
      c1.set_halfedge_handle(cv.halfedge_handle());
      c2.set_halfedge_handle(cv.halfedge_handle());
    }
  };

  Split_2 split_2_object () const
  {
    return (Split_2 (this->m_base_traits->split_2_object()));
  }
};

CGAL_END_NAMESPACE

#endif
