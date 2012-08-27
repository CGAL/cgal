// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion


#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < typename Point_3 >
class Orientation_3_benchmark

{

public:


  CGAL::Orientation
  operator()(const Point_3 &p, const Point_3 &q,
	     const Point_3 &r, const Point_3 &s) const
  {
      Get_approx<Point_3> get_approx; // Identity functor for all points
                                      // but lazy points.

       double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz;
      
      if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
          fit_in_double(get_approx(p).z(), pz) &&
          fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
          fit_in_double(get_approx(q).z(), qz) &&
          fit_in_double(get_approx(r).x(), rx) && fit_in_double(get_approx(r).y(), ry) &&
          fit_in_double(get_approx(r).z(), rz) &&
          fit_in_double(get_approx(s).x(), sx) && fit_in_double(get_approx(s).y(), sy) &&
          fit_in_double(get_approx(s).z(), sz))
      
      {
          double pqx = qx - px;
          double pqy = qy - py;
          double pqz = qz - pz;
          double prx = rx - px;
          double pry = ry - py;
          double prz = rz - pz;
          double psx = sx - px;
          double psy = sy - py;
          double psz = sz - pz;

          // CGAL::abs uses fabs on platforms where it is faster than (a<0)?-a:a
          // Then semi-static filter.
#if 1
          double maxx = CGAL::abs(pqx);
          double maxy = CGAL::abs(pqy);
          double maxz = CGAL::abs(pqz);

          double aprx = CGAL::abs(prx);
          double apsx = CGAL::abs(psx);

          double apry = CGAL::abs(pry);
          double apsy = CGAL::abs(psy);

          double aprz = CGAL::abs(prz);
          double apsz = CGAL::abs(psz);

          if (maxx < aprx) maxx = aprx;
          if (maxx < apsx) maxx = apsx;
          if (maxy < apry) maxy = apry;
          if (maxy < apsy) maxy = apsy;
          if (maxz < aprz) maxz = aprz;
          if (maxz < apsz) maxz = apsz;
#endif

#if 0
          // Don't call abs at all
          double maxx = pqx;
          double maxy = pqy;
          double maxz = pqz;

          double aprx = prx;
          double apsx = psx;

          double apry = pry;
          double apsy = psy;

          double aprz = prz;
          double apsz = psz;

          if (maxx < aprx) maxx = aprx;
          if (maxx < apsx) maxx = apsx;
          if (maxy < apry) maxy = apry;
          if (maxy < apsy) maxy = apsy;
          if (maxz < aprz) maxz = aprz;
          if (maxz < apsz) maxz = apsz;
#endif

#if 0
          // 3 less local variables, but no gain
          double maxx = CGAL::abs(pqx);
          double aprx = CGAL::abs(prx);
          if (maxx < aprx) maxx = aprx;
          aprx = CGAL::abs(psx);
          if (maxx < aprx) maxx = aprx;

          double maxy = CGAL::abs(pqy);
          double apry = CGAL::abs(pry);
          if (maxy < apry) maxy = apry;
          apry = CGAL::abs(psy);
          if (maxy < apry) maxy = apry;

          double maxz = CGAL::abs(pqz);
          double aprz = CGAL::abs(prz);
          if (maxz < aprz) maxz = aprz;
          aprz = CGAL::abs(psz);
          if (maxz < aprz) maxz = aprz;
#endif

#if 0
          // The slowest variant
          double maxx = (std::max)( (std::max)(CGAL::abs(pqx),CGAL::abs(prx)), CGAL::abs(psx));
          double maxy = (std::max)( (std::max)(CGAL::abs(pqy),CGAL::abs(pry)), CGAL::abs(psy));
          double maxz = (std::max)( (std::max)(CGAL::abs(pqz),CGAL::abs(prz)), CGAL::abs(psz));
#endif
          double det = 0; /* CGAL::determinant(pqx, pqy, pqz,
                                         prx, pry, prz,
                                         psx, psy, psz);*/
          //return ZERO; //  0.9 sec  because it becomes an almost empty function

          // "use" the max values
          det  = det + maxx + maxy + maxz;

          if(det >0) return POSITIVE;
          if(det < 0) return NEGATIVE;
          return ZERO; 
          // 4.9 sec without "use the .."  -  This is then just the determinant computation 
          // 10.2 sec with "use the ..."   -  This is determinant + abs + computing the max values
          // 7.3 sec with "use the ..", but with CGAL::abs commented

          // Sort maxx < maxy < maxz.
          if (maxx > maxz)
              std::swap(maxx, maxz);
          if (maxy > maxz)
              std::swap(maxy, maxz);
          else if (maxy < maxx)
              std::swap(maxx, maxy);

          // Protect against underflow in the computation of eps.
          if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
            if (maxx == 0)
              return ZERO;
          }
          // Protect against overflow in the computation of det.
          else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {
            double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;
            if (det > eps)  return POSITIVE;
            if (det < -eps) return NEGATIVE;
          }

	  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
      }
      std::cerr << "failure" << std::endl;
      return ZERO;  
      // 11 sec

      //return Base::operator()(p, q, r, s);
     
  }


};

} } } // namespace CGAL::internal::Static_filters_predicates

