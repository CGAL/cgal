// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
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
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_H
#define CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>
#include <cmath>
#include <iostream>

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base >
class Do_intersect_3
  : public K_base::Do_intersect_3
{
  typedef typename K_base::Point_3   Point_3;
  typedef typename K_base::Segment_3 Segment_3;

public:


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else 

  // TODO: add them all!!

#endif


 typedef typename Base::result_type  result_type;


  result_type 
  operator()(const Segment_3 &s, Bbox_3& b) const
  {
    CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Do_intersect_3", tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
    // but lazy points.
    const Point_3& p = s.source(); 
    const Point_3& q = s.target(); 

    double px, py, pz, qx, qy, qz;
    double bxmin = b.xmin(), bymin = b.ymin(), bzmin = b.zmin(), 
      bxmax = b.xmax(), bymax = b.ymax(), bzmax = b.zmax();

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

        
      }

    // AF: I copy pasted the code to call the max

    // -----------------------------------
    // treat x coord
    // -----------------------------------
    double dmin, tmin, tmax;
    if ( qx >= px )  // this is input and needs no epsilon
      {
        tmin = bxmin - px;
        tmax = bxmax - px;
        dmin = qx - px;
      }
    else
      {
        tmin = px - bxmax;
        tmax = px - bxmin;
        dmin = px - qx;
      }
   
    
    double m = CGAL::abs(tmin), m2;
    m2 = CGAL::abs(tmax); if(m2 > m) { m = m2; }
    m2 = CGAL::abs(dmin); if(m2 > m) { m = m2; }

    double error =  ERROR_FOR_ONE * m;
    // tmax is a difference so we should compare with an appropriate eps
    // We also should compare tmin-dmin with another eps

    if ( tmax < - error || tmin - dmin > error )
      return false;

    double dmax = dmin;
    if ( tmin < 0 )
      {
        tmin = 0;
        dmin = 1;
      }

    if ( tmax > dmax )
      {
        tmax = 1;
        dmax = 1;
      }

    // -----------------------------------
    // treat y coord
    // -----------------------------------
    double d_, tmin_, tmax_;
    if ( qy >= py )   // this is input and needs no epsilon
      {
        tmin_ = bymin - py;
        tmax_ = bymax - py;
        d_ = qy - py;
      }
    else
      {
        tmin_ = py - bymax;
        tmax_ = py - bymin;
        d_ = py - qy;
      }
    
    m2 = CGAL::abs(tmin_); if(m2 > m) { m = m2; }
    m2 = CGAL::abs(tmax_); if(m2 > m) { m = m2; }
    m2 = CGAL::abs(d_); if(m2 > m) { m = m2; }
    double error =  ERROR_FOR_ONE * m * m;
    // epsilons needed
    if (  (((d_*tmin) - (dmin*tmax_)) > error) || (((dmax*tmin_) - (d_*tmax)) > error) )
      return false;
    
    // epsilons needed
    if( ((dmin*tmin_) - (d_*tmin)) > error )
      {
        tmin = tmin_;
        dmin = d_;
      }
    
    // epsilons needed
    if(  ((d_*tmax) - (dmax*tmax_) > error ))
      {
        tmax = tmax_;
        dmax = d_;
      }
    
    // -----------------------------------
    // treat z coord
    // -----------------------------------
    if ( qz >= pz )   // this is input and needs no epsilon
      {
        tmin_ = bzmin - pz;
        tmax_ = bzmax - pz;
        d_ = qz - pz;
      }
    else
      {
        tmin_ = pz - bzmax;
        tmax_ = pz - bzmin;
        d_ = pz - qz;
      }
    
    m2 = CGAL::abs(tmin_); if(m2 > m) { m = m2; }
    m2 = CGAL::abs(tmax_); if(m2 > m) { m = m2; }
    m2 = CGAL::abs(d_); if(m2 > m) { m = m2; }
    double error =  ERROR_FOR_ONE * m * m;
    
    // epsilons needed
    if( (((dmin*tmax_) - (d_*tmin))> error) && (((d_*tmax) - (dmax*tmin_)) > error)  )
      {
        return true;
      }    
    return Base::operator()(s,b);
  }



    // Computes the epsilon for Bbox_3_Segment_3_do_intersect.
    static double compute_epsilon_bbox_segment_3()
    {
      typedef Static_filter_error F;
      F t1 = F(1);

      // TODO: write the most complex arithmetic expression that happens
      //       in the operator above 

      F f = ((t1 - t1) * (t1 - t1)) - ((t1 - t1) * (t1 - t1));
      
      double err = f.error();
      err += err * 2 *  F::ulp(); // Correction due to "eps * m * m".  Do we need 2 ?
      std::cerr << "*** epsilon for Do_intersect_3 = " << err << std::endl;
      return err;
    }

}; // class Do_intersect_3

}  // namespace Static_filters_predicates

} // namespace internal


} //namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_H
