// Copyright (c) 2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Remy Thomasse  <remy.thomasse@inria.fr>

#ifndef CGAL_RANDOM_CONVEX_HULL_TRAITS_2_H
#define CGAL_RANDOM_CONVEX_HULL_TRAITS_2_H

namespace CGAL{
	template<class R_>
class Random_convex_hull_traits_2 
{
  public:
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Segment_2		  Segment_2;
   	typedef typename R::Compare_x_2		  Compare_x_2;
   	typedef typename R::Compare_y_2		  Compare_y_2;
   	typedef typename R::Orientation_2	  Orientation_2;


    Compare_x_2
    compare_x_2_object() const
    {
    	return Compare_x_2(); 
    }

    Compare_y_2
    compare_y_2_object() const
    {
    	return Compare_y_2();
    }

    Orientation_2
    orientation_2_object() const 
    {
    	return Orientation_2();
    }

};

} //namespace CGAL
#endif //CGAL_RANDOM_CONVEX_HULL_TRAITS_2_H