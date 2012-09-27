// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Filtered_kernel/include/CGAL/Static_filters/Side_of_oriented_circle_2.h $
// $Id: Side_of_oriented_circle_2.h 47568 2008-12-21 15:56:20Z spion $
// 
//
// Author(s)     : 

#ifndef CGAL_STATIC_FILTERS_IN_CONE_3_H
#define CGAL_STATIC_FILTERS_IN_CONE_3_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/config.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>

namespace CGAL {

#define CGAL_EVALUATE_SIGN_AND_BREAK_IF_WE_HAVE_A_FAILURE(s, e, eps)\
\
xx = sign(e);\
if(xx == 0) break;\
if((xx < 0) && (xx >= -eps)) break;\
if((xx > 0) && (xx <= eps)) break;\
s = xx; \
\

/*long long hexify(double d)
{
  union {
    long long i;
    double d;
	} x;
	x.d = d;
	return x.i;
}*/

template < typename K_base >
class SF_In_cone_3
  : public K_base::In_cone_3
{
  typedef typename K_base::Point_2                      Point_2;
  typedef typename K_base::In_cone_3                    Base;

public:

	static double eps;
	
	using K_base::In_cone_3::operator();
	
  Oriented_side operator()(const Point_2 &p, const Point_2 &q,
	                   const Point_2 &s, const Point_2 &t) const
  {
      double px, py, pz, qx, qy, qz, sx, sy, sz, tx, ty, tz;

      if(internal::fit_in_double(p.x(), px) && internal::fit_in_double(p.y(), py) && internal::fit_in_double(p.z(), pz) &&
            internal::fit_in_double(q.x(), qx) && internal::fit_in_double(q.y(), qy) && internal::fit_in_double(q.z(), qz) &&
            internal::fit_in_double(s.x(), sx) && internal::fit_in_double(s.y(), sy) && internal::fit_in_double(s.z(), sz) &&
            internal::fit_in_double(t.x(), tx) && internal::fit_in_double(t.y(), ty) && internal::fit_in_double(t.z(), tz))
      {
	
	      using std::fabs;
	
	      const double A1 = determinant(qx, qy, qz, sx, sy, sz, tx, ty, tz);
	      const double A3 = determinant(px, py, pz, qx, qy, qz, tx, ty, tz);
	      const double A2 = -determinant(px, py, pz, sx, sy, sz, tx, ty, tz);
	      const double A4 = -determinant(px, py, pz, qx, qy, qz, sx, sy, sz);
	
	      const double G1 = CGAL_NTS square(px) + CGAL_NTS square(py) + CGAL_NTS square(pz);
        const double G2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz);
        const double G3 = CGAL_NTS square(sx) + CGAL_NTS square(sy) + CGAL_NTS square(sz);
	      const double G4 = CGAL_NTS square(tx) + CGAL_NTS square(ty) + CGAL_NTS square(tz);
	
	      // Comme l'erreur relative de G1 < 1, on peut utiliser du sqrt
	
        const double res = (CGAL_NTS sqrt(G1)) * A1 +
                     (CGAL_NTS sqrt(G2)) * A2 +
                     (CGAL_NTS sqrt(G3)) * A3 +
										 (CGAL_NTS sqrt(G4)) * A4;

        // We compute the semi-static bound.
				const double abspx = fabs(px);
        const double absqx = fabs(qx);
        const double abssx = fabs(sx);
        const double abstx = fabs(tx);
        const double abspy = fabs(py);
        const double absqy = fabs(qy);
        const double abssy = fabs(sy);
        const double absty = fabs(ty);
        const double abspz = fabs(pz);
        const double absqz = fabs(qz);
        const double abssz = fabs(sz);
        const double abstz = fabs(tz);

        double maxx = abspx;
        if (maxx < absqx) maxx = absqx;
        if (maxx < abssx) maxx = abssx;
        if (maxx < abstx) maxx = abstx;

        double maxy = abspy;
        if (maxy < absqy) maxy = absqy;
        if (maxy < abssy) maxy = abssy;
        if (maxy < absty) maxy = absty;

        double maxz = abspz;
        if (maxz < absqz) maxz = absqz;
        if (maxz < abssz) maxz = abssz;
        if (maxz < abstz) maxz = abstz;

         // Sort maxx < maxy < maxz.
         if (maxx > maxz)
             std::swap(maxx, maxz);
         if (maxy > maxz)
             std::swap(maxy, maxz);
         else if (maxy < maxx)
             std::swap(maxx, maxy); 

       // Protect against underflow in the computation of eps.
        if (maxx < 1e-73) {
          if(maxx == 0) return DEGENERATE;
        }	else if (maxx < 1e76) {
				  double epsb = eps * maxx * maxy * maxz * maxz;
          if (res > epsb)  {return COUNTERCLOCKWISE;}
					if (res < -epsb) {return CLOCKWISE;}
				}
      }
			/*std::cout << "example" << std::endl;
			std::cout << "0x" << std::hex << hexify(px) << " "<< "0x" << std::hex << hexify(py) << " " << "0x" << std::hex << hexify(pz) << std::endl;
			std::cout << "0x" << std::hex << hexify(qx) << " "<< "0x" << std::hex << hexify(qy) << " " << "0x" << std::hex << hexify(qz) << std::endl;
			std::cout << "0x" << std::hex << hexify(sx) << " "<< "0x" << std::hex << hexify(sy) << " " << "0x" << std::hex << hexify(sz) << std::endl;
			std::cout << "0x" << std::hex << hexify(tx) << " "<< "0x" << std::hex << hexify(ty) << " " << "0x" << std::hex << hexify(tz) << std::endl;*/
      return Base::operator()(p, q, s, t);
  }

  // Computes several epsilons for In_cone_3.
  // on supposera que nos coordonnees resteront plus petites que 2.
  static void compute_epsilons()
  {
    typedef CGAL::internal::Static_filter_error F;

		double err;

    F t1 = F(1, F::ulp()/2);         // First translation
    
    F A1 = determinant(t1, t1, t1, t1, t1, t1, t1, t1, t1);
    F G1 = t1*t1 + t1*t1 + t1*t1;
		F sqrtG1 = CGAL_NTS sqrt(G1);
		F res = A1 * sqrtG1 + A1 * sqrtG1 + A1 * sqrtG1 + A1 * sqrtG1;

		// the degree is 4
		
		// eps
		err = res.error();    
		err += err * 4 * F::ulp(); // Correction due to "eps * maxx * ...".
    std::cerr << "*** eps for In_cone_3 = " << std::setprecision(16) << err << std::endl;
	
  }
};

template < typename K_base >
double SF_In_cone_3<K_base>::eps = 9.46507951849212e-07;

} //namespace CGAL

#endif // CGAL_STATIC_FILTERS_IN_CONE_3_H