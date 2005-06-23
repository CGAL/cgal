// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef SWEEP_LINE_2_UTILS
#define SWEEP_LINE_2_UTILS

CGAL_BEGIN_NAMESPACE

#include <CGAL/Object.h>
#include <CGAL/assertions.h>
#include <vector>
#include <algorithm>

template <class Traits,
          class CurveInputIter,
          class XCurveOutIter,
          class PointOutIter>
void  make_x_monotone(CurveInputIter begin,
                      CurveInputIter end,
                      XCurveOutIter xcurves_out,
                      PointOutIter points_out,
                      const Traits* tr)
{
  unsigned int num_of_curves = std::distance(begin, end);
  std::vector<Object> object_vec;
  object_vec.reserve(num_of_curves);
  for(CurveInputIter iter = begin; iter != end; ++iter)
  {
    tr->make_x_monotone_2_object()(*iter, std::back_inserter(object_vec));
  }

  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits::Point_2               Point_2;

  X_monotone_curve_2 xcurve;

  for(unsigned int i = 0 ; i < object_vec.size() ; ++i)
  {
    Object& obj = object_vec[i];
    if(assign(xcurve, obj))
      *xcurves_out++ = xcurve;
    else
    {
      Point_2            pt;
      CGAL_assertion(assign(pt, obj));
      assign(pt,obj);
      *points_out++ = pt;
    }
  }
}



CGAL_END_NAMESPACE

#endif
