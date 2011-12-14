// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
// Copyright (c) 2011 GeometryFactory Sarl (France)
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
// Author(s)     : Andreas Fabri, Laurent Rineau


#ifndef CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_3_H

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


template < typename K_base, typename SFK >
class Do_intersect_3
  : public K_base::Do_intersect_3
{
  typedef typename K_base::Point_3   Point_3;
  typedef typename K_base::Ray_3     Ray_3;
  typedef typename K_base::Segment_3 Segment_3;
  typedef typename K_base::Triangle_3 Triangle_3;
  typedef typename K_base::Do_intersect_3 Base;

public:

  typedef typename Base::result_type  result_type;


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif // CGAL_CFG_MATCHING_BUG_6


  Sign sign_with_error(const double x, const double error) const {
    if(x > error) return POSITIVE;
    else if( x < - error) return NEGATIVE;
    else return ZERO;
  }


  // The internal::do_intersect(..) function 
  // only performs orientation tests on the vertices
  // of the triangle and the segment
  // By calling the do_intersect function with
  // the  statically filtered kernel we avoid
  // that doubles are put into Inteval_nt
  // to get taken out again with fit_in_double
  result_type 
  operator()(const Segment_3 &s, const Triangle_3& t) const
  {
    return internal::do_intersect(t,s, SFK());
  }

  result_type 
  operator()(const Triangle_3& t, const Segment_3 &s) const
  {
    return internal::do_intersect(t,s, SFK());
  }

  result_type 
  operator()(const Bbox_3& b, const Segment_3 &s) const
  {
    return this->operator()(s, b);
  }

  result_type 
  operator()(const Segment_3 &s, const Bbox_3& b) const
  {
    CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
                           std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
                                    // but lazy points
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

      // -----------------------------------
      // treat x coord
      // -----------------------------------
      double dmin, dmax, tmin, tmax;
      if ( qx >= px )  // this is input and needs no epsilon
      {
        if(px > bxmax) return false; // segment on the right of bbox
        if(qx < bxmin) return false; // segment on the left of bbox

        tmax = bxmax - px;
        
        dmax = qx - px;
        if ( bxmin < px ) // tmin < 0 means px is in the x-range of bbox
        {
          tmin = 0;
          dmin = 1;
        } else {
          tmin = bxmin - px;
          dmin = dmax;
        }
      }
      else
      {
        if(qx > bxmax) return false; // segment on the right of bbox
        if(px < bxmin) return false; // segment on the left of bbox

        tmax = px - bxmin;

        dmax = px - qx;
        if ( px < bxmax ) // tmin < 0 means px is in the x-range of bbox
        {
          tmin = 0;
          dmin = 1;
        } else {
          tmin = px - bxmax;
          dmin = dmax;
        }
      }   

      double m = CGAL::abs(tmin), m2;
      m2 = CGAL::abs(tmax); if(m2 > m) { m = m2; }
      m2 = CGAL::abs(dmin); if(m2 > m) { m = m2; }

      if(m < 7e-294) {
        // underflow in the computation of 'error'
        return Base::operator()(s,b);
      }

      const double EPS_1 = 3.55618e-15;

      double error =  EPS_1 * m;

      switch(sign_with_error( tmax - dmax, error)) {
      case POSITIVE:
        tmax = 1;
        dmax = 1;
        break;
      case NEGATIVE:
        break;
      default:
        // ambiguity of the comparison tmax > dmin
        // let's call the exact predicate
        return Base::operator()(s,b);
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

      if(m < 3e-147) {
        // underflow in the computation of 'error'
        return Base::operator()(s,b);
      }
      
      error =  EPS_1 * m * m;

      if(m > 1e153) { /* sqrt(max_double [hadamard]/2) */ 
        // potential overflow on the computation of 'sign1' and 'sign2'
        return Base::operator()(s,b);
      }
      Sign sign1 = sign_with_error( (d_*tmin) - (dmin*tmax_) , error);
      Sign sign2 = sign_with_error( (dmax*tmin_) - (d_*tmax) , error);

      if(sign1 == POSITIVE || sign2 == POSITIVE) 
        return false; // We are *sure* the segment is outside the box, on one
                      // side or the other.
      if(sign1 == ZERO || sign2 == ZERO) {
        return Base::operator()(s,b); // We are *unsure*: one *may be*
                                      // positive.
      }

      // Here we are sure the two signs are negative. We can continue with
      // the rest of the function...

      // epsilons needed
      switch(sign_with_error((dmin*tmin_) - (d_*tmin) , error)) {
      case POSITIVE:
        tmin = tmin_;
        dmin = d_;
        break;
      case NEGATIVE:
        break;
      default: // uncertainty
        return Base::operator()(s,b);
      }
    
      // epsilons needed
      switch(sign_with_error((d_*tmax) - (dmax*tmax_) , error)) {
      case POSITIVE:
        tmax = tmax_;
        dmax = d_;
      case NEGATIVE:
        break;
      default: // uncertainty
        return Base::operator()(s,b);
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

      // m may have changed
      error =  EPS_1 * m * m;
    
      if(m > 1e153) { /* sqrt(max_double [hadamard]/2) */ 
        // potential overflow on the computation of 'sign1' and 'sign2'
        return Base::operator()(s,b);
      }
      sign1 = sign_with_error( (dmin*tmax_) - (d_*tmin) , error);
      sign2 = sign_with_error( (d_*tmax) - (dmax*tmin_) , error);
      if(sign1 == NEGATIVE || sign2 == NEGATIVE) {
        return false; // We are *sure* the segment is outside the box, on one
                      // side or the other.
      }
      if(sign1 == ZERO || sign2 == ZERO) {
        return Base::operator()(s,b); // We are *unsure*: one *may be*
                                      // negative.
      }
      return true; // We are *sure* the two signs are positive.
    }
    return Base::operator()(s,b);
  }



  result_type 
  operator()(const Bbox_3& b, const Ray_3 &r) const
  {
    return this->operator()(r, b);
  }

  result_type 
  operator()(const Ray_3 &r, const Bbox_3& b) const
  {
    CGAL_BRANCH_PROFILER_3(std::string("semi-static failures/attempts/calls to   : ") +
                           std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Point_3> get_approx; // Identity functor for all points
    // but lazy points.
    const Point_3& p = r.source(); 
    const Point_3& q = r.second_point(); 

    double px, py, pz, qx, qy, qz;
    double bxmin = b.xmin(), bymin = b.ymin(), bzmin = b.zmin(), 
      bxmax = b.xmax(), bymax = b.ymax(), bzmax = b.zmax();

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) )   
    {
      CGAL_BRANCH_PROFILER_BRANCH_1(tmp);

      // -----------------------------------
      // treat x coord
      // -----------------------------------
      double dmin, dmax, tmin, tmax;
      if ( qx >= px )  // this is input and needs no epsilon
      {
        if(px > bxmax) return false; // if tmax < 0, ray on the right of bbox
        tmax = bxmax - px;
        
        dmax = qx - px;
        if ( bxmin < px ) // tmin < 0 means px is in the x-range of bbox
        {
          tmin = 0;
          dmin = 1;
        } else {
          if( px == qx ) return false; // if dmin == 0
          tmin = bxmin - px;
          dmin = dmax;
        }
      }
      else
      {
        if(px < bxmin) return false; // if tmax < 0, ray on the left of bbox
        tmax = px - bxmin;

        dmax = px - qx;
        if ( px < bxmax ) // tmin < 0 means px is in the x-range of bbox
        {
          tmin = 0;
          dmin = 1;
        } else {
          if( px == qx ) return false; // if dmin == 0
          tmin = px - bxmax;
          dmin = dmax;
        }
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
    
      double m = CGAL::abs(tmin), m2;
      m2 = CGAL::abs(tmax); if(m2 > m) { m = m2; }
      m2 = CGAL::abs(dmin); if(m2 > m) { m = m2; }
      m2 = CGAL::abs(tmin_); if(m2 > m) { m = m2; }
      m2 = CGAL::abs(tmax_); if(m2 > m) { m = m2; }
      m2 = CGAL::abs(d_); if(m2 > m) { m = m2; }

      if(m < 3e-147) {
        // underflow in the computation of 'error'
        return Base::operator()(r,b);
      }
      
      const double EPS_1 = 3.55618e-15;

      double error =  EPS_1 * m * m;

      if(m > 1e153) { /* sqrt(max_double [hadamard]/2) */ 
        // potential overflow on the computation of 'sign1' and 'sign2'
        return Base::operator()(r,b);
      }
      Sign sign1 = sign_with_error( (d_*tmin) - (dmin*tmax_) , error);
      Sign sign2 = sign_with_error( (dmax*tmin_) - (d_*tmax) , error);

      if(sign1 == POSITIVE || sign2 == POSITIVE) 
        return false; // We are *sure* the ray is outside the box, on one
                      // side or the other.
      if(sign1 == ZERO || sign2 == ZERO) {
        return Base::operator()(r,b); // We are *unsure*: one *may be*
                                      // positive.
      }

      // Here we are sure the two signs are negative. We can continue with
      // the rest of the function...

      // epsilons needed
      switch(sign_with_error((dmin*tmin_) - (d_*tmin) , error)) {
      case POSITIVE:
        tmin = tmin_;
        dmin = d_;
        break;
      case NEGATIVE:
        break;
      default: // uncertainty
        return Base::operator()(r,b);
      }
    
      // epsilons needed
      switch(sign_with_error((d_*tmax) - (dmax*tmax_) , error)) {
      case POSITIVE:
        tmax = tmax_;
        dmax = d_;
      case NEGATIVE:
        break;
      default: // uncertainty
        return Base::operator()(r,b);
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

      // m may have changed
      error =  EPS_1 * m * m;
    
      if(m > 1e153) { /* sqrt(max_double [hadamard]/2) */ 
        // potential overflow on the computation of 'sign1' and 'sign2'
        return Base::operator()(r,b);
      }
      sign1 = sign_with_error( (dmin*tmax_) - (d_*tmin) , error);
      sign2 = sign_with_error( (d_*tmax) - (dmax*tmin_) , error);
      if(sign1 == NEGATIVE || sign2 == NEGATIVE) {
        return false; // We are *sure* the ray is outside the box, on one
                      // side or the other.
      }
      if(sign1 == ZERO || sign2 == ZERO) {
        return Base::operator()(r,b); // We are *unsure*: one *may be*
                                      // negative.
      }
      return true; // We are *sure* the two signs are positive.
    }
    return Base::operator()(r,b);
  }



  // Computes the epsilon for Bbox_3_Segment_3_do_intersect.
  static double compute_epsilon_bbox_segment_3()
  {
    typedef Static_filter_error F;
    F t1 = F(1);
    F f = ((t1 - t1) * (t1 - t1)) - ((t1 - t1) * (t1 - t1));
    F f1 = (t1 - t1);
    F f1bis = (t1 - t1) - (t1 - t1);
    F f2 = f1*f1;
    F f3 = f2 - f2;
    std::cerr << "epsilons:\n"
              << "  degre " << f1.degree() << ": " <<  f1.error() << "\n"
              << "  degre " << f1bis.degree() << ": " <<  f1bis.error() << "\n"
              << "  degre " << f2.degree() << ": " <<  f2.error() << "\n"
              << "  degre " << f3.degree() << ": " <<  f3.error() << "\n";
      
    double err = f.error();
    err += err * 2 *  F::ulp(); // Correction due to "eps * m * m".  Do we need 2 ?
    std::cerr << "*** epsilon for Do_intersect_3(Bbox_3, Segment_3) = "
              << err << std::endl;
    std::cerr << "\n"
              << "Now for underflow/overflows...\n"
              << "        min_double/eps = " 
              << std::numeric_limits<double>::min() / err << std::endl
              << "  sqrt(min_double/eps) = "
              << CGAL::sqrt(std::numeric_limits<double>::min() / err) << std::endl;
    return err;
  }

}; // class Do_intersect_3

} // end namespace Static_filters_predicates

} // end namespace internal


} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_3_H
