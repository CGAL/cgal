// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Sebastien Loriot
// This predicates was generated using the fpg tool written by Andreas Meyer.
//

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COMPARE_WEIGHTED_SQUARED_RADIUS_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COMPARE_WEIGHTED_SQUARED_RADIUS_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <cmath>


namespace CGAL { namespace internal { namespace Static_filters_predicates {
  
  
  template <typename K_base>
  class Compare_weighted_squared_radius_3:
    public K_base::Compare_weighted_squared_radius_3
  {
    typedef typename K_base::Weighted_point_3 Weighted_point_3;
    typedef typename K_base::FT FT;
    typedef typename K_base::Compare_weighted_squared_radius_3 Base;
  public:
    typedef typename Base::result_type result_type;
  
    #ifndef CGAL_CFG_MATCHING_BUG_6
    using Base::operator();
    #else 
    result_type
    operator()(const Weighted_point_3& p, const FT& w) const
    { 
      return Base::operator()(p,w);
    }
    #endif
    
    
    result_type operator() (
        const Weighted_point_3& p, 
        const Weighted_point_3& q, 
        const Weighted_point_3& r, 
        const Weighted_point_3& s,
        const FT& w
    ) const {

      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_weighted_squared_radius_3 with 4 wpoints", tmp);
      
      double px, py, pz, pw, qx, qy, qz, qw, rx, ry, rz, rw, sx, sy, sz, sw, alpha;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pw) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qw) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)      &&
          fit_in_double(r.z(), rz) && fit_in_double(r.weight(), rw) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy)      &&
          fit_in_double(s.z(), sz) && fit_in_double(s.weight(), sw) && 
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        
        double qpx = (qx - px);
        double qpy = (qy - py);
        double qpz = (qz - pz);
        double pw_qw = (pw - qw);
        double qp2 = (((CGAL::square( qpx ) + CGAL::square( qpy )) + CGAL::square( qpz )) + pw_qw);
        double rpx = (rx - px);
        double rpy = (ry - py);
        double rpz = (rz - pz);
        double pw_rw = (pw - rw);
        double rp2 = (((CGAL::square( rpx ) + CGAL::square( rpy )) + CGAL::square( rpz )) + pw_rw);
        double spx = (sx - px);
        double spy = (sy - py);
        double spz = (sz - pz);
        double pw_sw = (pw - sw);
        double sp2 = (((CGAL::square( spx ) + CGAL::square( spy )) + CGAL::square( spz )) + pw_sw);
        double num_x = CGAL::determinant( qpy, qpz, qp2, rpy, rpz, rp2, spy, spz, sp2 );
        double num_y = CGAL::determinant( qpx, qpz, qp2, rpx, rpz, rp2, spx, spz, sp2 );
        double num_z = CGAL::determinant( qpx, qpy, qp2, rpx, rpy, rp2, spx, spy, sp2 );
        double den = CGAL::determinant( qpx, qpy, qpz, rpx, rpy, rpz, spx, spy, spz );
        double alpha_pw = (alpha + pw);
        CGAL::Sign int_tmp_result;
        double eps;
        double double_tmp_result = (((alpha_pw * 4.00000000000000000000e+00) * CGAL::square( den )) - ((CGAL::square( num_x ) + CGAL::square( num_y )) + CGAL::square( num_z )));
        double max4 = CGAL::abs(qpy);

        double aqpx = CGAL::abs(qpx);

        double arpx = CGAL::abs(rpx);
        double arpy = CGAL::abs(rpy);
        double arpz = CGAL::abs(rpz);

        double aspx = CGAL::abs(spx);
        double aspy = CGAL::abs(spy);
        double aspz = CGAL::abs(spz);

        double apw_rw = CGAL::abs(pw_rw);
        double apw_sw = CGAL::abs(pw_sw);

        double aalpha_pw = CGAL::abs(alpha_pw);

        if( max4 < arpy ) max4 = arpy;
        if( max4 < aspy ) max4 = aspy; 
        double max2 = max4;
        if( max2 < aqpx ) max2 = aqpx; 
        if( max2 < arpx ) max2 = arpx; 
        if( max2 < aspx ) max2 = aspx; 
        double max1 = max2;
        double max3 = max4;
        double max5 = CGAL::abs(qpz);
        if( max5 < arpz ) max5 = arpz;  
        if( max5 < aspz ) max5 = aspz;  
        if( max3 < max5 )      max3 = max5;       
        if( max1 < max3 )      max1 = max3;       
        if( max1 < max4 )      max1 = max4;       
        if( max1 < max5 )      max1 = max5;      
        double max6 = CGAL::abs(pw_qw);
        if( max6 < apw_rw ) max6 = apw_rw;
        if( max6 < apw_sw ) max6 = apw_sw;
        double max7 = max6;
        if( max7 < aalpha_pw ) max7 = aalpha_pw; 
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( max2 < lower_bound_1 ) lower_bound_1 = max2;
        if( max3 < lower_bound_1 ) lower_bound_1 = max3;
        if( max4 < lower_bound_1 ) lower_bound_1 = max4;
        if( max5 < lower_bound_1 ) lower_bound_1 = max5;
        double lower_bound_2;
        double upper_bound_2;
        lower_bound_2 = max6;
        upper_bound_2 = max6;
        if( max7 < lower_bound_2 ) lower_bound_2 = max7;
        else 
        {
            if( max7 > upper_bound_2 ) upper_bound_2 = max7;
        } 
        if( ((lower_bound_1 < 8.99995159231796093217e-38) || (lower_bound_2 < 8.09991286640666077573e-75)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,s,w);
        } 
        else 
        {
            if( ((upper_bound_1 > 2.59614842926741294957e+33) || (upper_bound_2 > 6.73998666678765545893e+66)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,s,w);
            } 
            eps = (5.16919651938288206243e-12 * (((((max2 * max3) * max1) * max1) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max7, (max1 * max1) )) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max6, (max4 * max5) ), (max1 * max1) )));
            if( (double_tmp_result > eps) )
            {
                int_tmp_result = CGAL::NEGATIVE;
            } 
            else 
            {
                if( (double_tmp_result < -eps) )
                {
                    int_tmp_result = CGAL::POSITIVE;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,s,w);
                } 
            } 
        } 
        return int_tmp_result;
      }
      else
        return Base::operator()(p,q,r,s,w);
    }

    result_type
    operator() (
      const Weighted_point_3& p, 
      const Weighted_point_3& q , 
      const Weighted_point_3& r,
      const FT& w
    ) const {
      
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_weighted_squared_radius_3 with 3 wpoints", tmp);
      
      double px, py, pz, pw, qx, qy, qz, qw, rx, ry, rz, rw, alpha;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pw) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qw) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)      &&
          fit_in_double(r.z(), rz) && fit_in_double(r.weight(), rw) &&
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        
        double qpx = (qx - px);
        double qpy = (qy - py);
        double qpz = (qz - pz);
        double pw_qw = (pw - qw);
        double qp2 = (((CGAL::square( qpx ) + CGAL::square( qpy )) + CGAL::square( qpz )) + pw_qw);
        double rpx = (rx - px);
        double rpy = (ry - py);
        double rpz = (rz - pz);
        double pw_rw = (pw - rw);
        double rp2 = (((CGAL::square( rpx ) + CGAL::square( rpy )) + CGAL::square( rpz )) + pw_rw);
        double sx = ((qpy * rpz) - (qpz * rpy));
        double sy = ((qpz * rpx) - (qpx * rpz));
        double sz = ((qpx * rpy) - (qpy * rpx));
        double num_x = ((qp2 * CGAL::determinant( rpy, rpz, sy, sz )) - (rp2 * CGAL::determinant( qpy, qpz, sy, sz )));
        double num_y = ((qp2 * CGAL::determinant( rpx, rpz, sx, sz )) - (rp2 * CGAL::determinant( qpx, qpz, sx, sz )));
        double num_z = ((qp2 * CGAL::determinant( rpx, rpy, sx, sy )) - (rp2 * CGAL::determinant( qpx, qpy, sx, sy )));
        double den = CGAL::determinant( qpx, qpy, qpz, rpx, rpy, rpz, sx, sy, sz );
        double alpha_pw = (alpha + pw);
        CGAL::Sign int_tmp_result;
        double eps;
        double double_tmp_result = ((alpha_pw * CGAL::square( den )) - (2.50000000000000000000e-01 * ((CGAL::square( num_x ) + CGAL::square( num_y )) + CGAL::square( num_z ))));
        double max2 = CGAL::abs(qpx);

        double aqpy = CGAL::abs(qpy);
        double aqpz = CGAL::abs(qpz);

        double arpx = CGAL::abs(rpx);
        double arpy = CGAL::abs(rpy);
        double arpz = CGAL::abs(rpz);

        double apw_rw = CGAL::abs(pw_rw);

        double aalpha_pw = CGAL::abs(alpha_pw);

        if( max2 < aqpy )  max2 = aqpy;
        if( max2 < arpx )  max2 = arpx;
        if( max2 < arpy )  max2 = arpy;
        double max1 = max2;
        if( max1 < aqpz )  max1 = aqpz;
        if( max1 < arpz )  max1 = arpz;
        double max3 = CGAL::abs(pw_qw);
        if( max3 < apw_rw ) max3 = apw_rw;
        double max4 = max3;
        if( max4 < aalpha_pw ) max4 = aalpha_pw;
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( max2 < lower_bound_1 ) lower_bound_1 = max2;
        double lower_bound_2;
        double upper_bound_2;
        lower_bound_2 = max3;
        upper_bound_2 = max3;
        if( max4 < lower_bound_2 )  lower_bound_2 = max4;
        else 
        {
          if( max4 > upper_bound_2 ) upper_bound_2 = max4;
        } 
        if( ((lower_bound_1 < 2.13354839219409151500e-30) || (lower_bound_2 < 4.55202874183399304187e-60)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,w);
        } 
        else 
        {
            if( ((upper_bound_1 > 1.23794003928538000002e+27) || (upper_bound_2 > 1.53249554086588817779e+54)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,w);
            } 
            eps = (1.13846439714120896721e-11 * (((((((max1 * max2) * max1) * max1) * max1) * max1) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max4, (max1 * max1) )) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max3, (max1 * max1) )));
            if( (double_tmp_result > eps) )
            {
                int_tmp_result = CGAL::NEGATIVE;
            } 
            else 
            {
                if( (double_tmp_result < -eps) )
                {
                    int_tmp_result = CGAL::POSITIVE;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,w);
                } 
            } 
        } 
        return int_tmp_result;
      }
      else
        return Base::operator()(p,q,r,w);
    }

    result_type
    operator() (
      const Weighted_point_3& p, 
      const Weighted_point_3& q,
      const FT& w
    ) const {
      
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_weighted_squared_radius_3 with 2 wpoints", tmp);
      
      double px, py, pz, pw, qx, qy, qz, qw, alpha;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pw) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qw) &&
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        
        double qpx = (qx - px);
        double qpy = (qy - py);
        double qpz = (qz - pz);
        double qp2 = ((CGAL::square( qpx ) + CGAL::square( qpy )) + CGAL::square( qpz ));
        double alpha_pw = (alpha + pw);
        double pw_qw = (pw - qw);
        CGAL::Sign int_tmp_result;
        double eps;
        double double_tmp_result = (((4.00000000000000000000e+00 * qp2) * alpha_pw) - CGAL::square( (qp2 + pw_qw) ));
        double max1 = CGAL::abs(qpx);

        double aqpy = CGAL::abs(qpy);
        double aqpz = CGAL::abs(qpz);
        double aalpha_pw = CGAL::abs(alpha_pw);

        if( max1 < aqpy ) max1 = aqpy;
        if( max1 < aqpz ) max1 = aqpz;
        double max2;
        double max3 = CGAL::abs(pw_qw);
        max2 = max3;
        if( max2 < aalpha_pw ) max2 = aalpha_pw;
        double lower_bound_2;
        double upper_bound_2;
        lower_bound_2 = max2;
        upper_bound_2 = max2;
        if( (max3 < lower_bound_2) ) lower_bound_2 = max3;
        if( ((max1 < 3.12497063152273477754e-74) || (lower_bound_2 < 9.76544144787960039738e-148)) ){
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,w);
        }
        else 
        {
            if( ((max1 > 6.42775217703595896130e+60) || (upper_bound_2 > 4.13159980493905099125e+121)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,w);
            } 
            eps = (2.33324675561911025753e-14 * (CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max2, (max1 * max1) ) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max3, (max1 * max1) )));
            if( (double_tmp_result > eps) )
            {
                int_tmp_result = CGAL::NEGATIVE;
            } 
            else 
            {
                if( (double_tmp_result < -eps) )
                {
                    int_tmp_result = CGAL::POSITIVE;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,w);
                } 
            } 
        } 
        return int_tmp_result;        
      }
      else
        return Base::operator()(p,q,w);
    }
    
  };
  
  
    } } }//namespace CGAL::internal::Static_filters_predicates


#endif //CGAL_INTERNAL_STATIC_FILTERS_COMPARE_WEIGHTED_SQUARED_RADIUS_3_H
