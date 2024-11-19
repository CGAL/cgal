// Copyright (c) 2009 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Sebastien Loriot
//

template <class R>
bool
_test_compare_angle_3(const R& rep)
{
  typedef typename R::FT FT;
  typedef typename R::Point_3 Point_3;
  typedef typename R::Vector_3 Vector_3;

  typename R::Compare_angle_3 compare_angle
    = rep.compare_angle_3_object();

  for(int theta1 = -170; theta1 <= 180; theta1+= 10)
  {
    const double angle1 = CGAL_PI*theta1/180.;
    Point_3 a(1, 0, 0);
    Point_3 b(0, 0, 0);
    Point_3 c((int)(std::cos(angle1)*1000), (int)(std::sin(angle1)*1000), 0);

    for(int theta2 = -170; theta2 <= 180; theta2+= 10) {
      if (theta1!=0 && theta1!=180 && abs(theta1)==abs(theta2)) continue;
      const double angle2 = CGAL_PI*theta2/180.;
      if ( CGAL::compare(abs(theta1), abs(theta2)) != CGAL::compare_angle(a, b, c, FT(std::cos(angle2))) )
        return false;
      if ( CGAL::compare(abs(theta1), abs(theta2)) != compare_angle(a, b, c, FT(std::cos(angle2))) )
        return false;

     Point_3 d((int)(std::cos(angle2)*1000), (int)(std::sin(angle2)*1000), 0);
      if ( CGAL::compare(abs(theta1), abs(theta2)) != CGAL::compare_angle(a, b, c, a, b, d ) )
        return false;
      if ( CGAL::compare(abs(theta1), abs(theta2)) != compare_angle(a, b, c, a, b, d ) )
        return false;

      Vector_3 u1(b, a), v1(b, c), v2(b, d);
      if ( CGAL::compare(abs(theta1), abs(theta2)) != CGAL::compare_angle(u1, v1, u1, v2) )
        return false;
      if ( CGAL::compare(abs(theta1), abs(theta2)) != compare_angle(u1, v1, u1, v2) )
        return false;
    } // end loop on theta2
  } // end loop and theta1
  return true;
}
