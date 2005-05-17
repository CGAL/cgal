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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef ARR_NON_X_TRAITS_H
#define ARR_NON_X_TRAITS_H


CGAL_BEGIN_NAMESPACE

// a meta-traits class that get a traits class as a template parameter ,
// and creates a new traits class that avoid intersection between curves.

template <class Traits>
class Arr_non_x_traits : public Traits
{

public:

  typedef typename Traits::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename Traits::Point_2                  Point_2;


  class Intersect_2
  {
  public:
   
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2();
  }
};


CGAL_END_NAMESPACE

#endif


