// Copyright (c) 2013  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// 
//
// Author(s)     : Olivier Devillers


template<class Kernel>
typename Kernel::Circle_2
compute_circle_in_pencil(typename Kernel::Circle_2 c,typename Kernel::Circle_2 c1,typename Kernel::Circle_2 c2)
{
  typename Kernel::Point_2 origin=CGAL::ORIGIN;
  typename Kernel::FT lambda = CGAL::squared_distance(c.center(),origin) ;
  lambda -= c.squared_radius() ;
  typename Kernel::FT l1 = CGAL::squared_distance(c1.center(),origin) - c1.squared_radius() ;
  typename Kernel::FT l2 = CGAL::squared_distance(c2.center(),origin) - c2.squared_radius() ;
  l1 += -2*((c1.center()-origin)*(c.center()-origin));
  l2 += -2*((c2.center()-origin)*(c.center()-origin));
  if (l1==l2){ // degenerate case, radical axis
    return  typename Kernel::Circle_2();
  }
  lambda= - (lambda+l2)/(l1-l2);
  typename Kernel::Point_2 center=origin+lambda*(c1.center()-origin)+(1-lambda)*(c2.center()-origin);
  typename Kernel::FT sqradius=-    lambda*(CGAL::squared_distance(c1.center(),origin)-c1.squared_radius())
    -(1-lambda)*(CGAL::squared_distance(c2.center(),origin)-c2.squared_radius())
    +CGAL::squared_distance(center,origin);
  typename Kernel::Circle_2 circ(center,sqradius);
  return circ;
}


template<class Kernel>
typename Kernel::Circle_2
compute_circle_orthogonal(typename Kernel::Circle_2 c,typename Kernel::Circle_2 c1,typename Kernel::Circle_2 c2)
{
  typename Kernel::Point_2 origin=CGAL::ORIGIN;
  typename Kernel::FT z= CGAL::squared_distance(c.center(),origin)-c.squared_radius();
  typename Kernel::FT z1= CGAL::squared_distance(c1.center(),origin)-c1.squared_radius();
  typename Kernel::FT z2= CGAL::squared_distance(c2.center(),origin)-c2.squared_radius();
  typename Kernel::FT det=-(c1.center().x() * c2.center().y() - c1.center().y() * c2.center().x())
                          +(c.center().x() * c2.center().y() - c.center().y() * c2.center().x())
                          -(c.center().x() * c1.center().y() - c.center().y() * c1.center().x());
  if (det==0.0){ // degenerate casse, radical axis
    return  typename Kernel::Circle_2();
  }
  typename Kernel::FT x=(  -(z1 * c2.center().y() - c1.center().y() * z2)
	                +(z * c2.center().y() - c.center().y() * z2)
	                -(z * c1.center().y() - c.center().y() * z1))/2/det;
  typename Kernel::FT y=(  -(c1.center().x() * z2 - z1 * c2.center().x())
	                +(c.center().x() * z2 - z * c2.center().x())
	                -(c.center().x() * z1 - z * c1.center().x()))/2/det;
  typename Kernel::FT rr=-(  (c1.center().x() * c2.center().y() - c1.center().y() * c2.center().x())*z
     		         -(c.center().x() * c2.center().y() - c.center().y() * c2.center().x())*z1
		         +(c.center().x() * c1.center().y() - c.center().y() * c1.center().x())*z2)/det+x*x+y*y;
  typename Kernel::Point_2 center(x,y);
  typename Kernel::Circle_2 circ(center,rr);
  return circ;
}
