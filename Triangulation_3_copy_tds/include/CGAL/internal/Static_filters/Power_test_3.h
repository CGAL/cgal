// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Sebastien Loriot
// This predicates was generated using the fpg tool written by Andreas Meyer.
//

#ifndef CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <cmath>


namespace CGAL { namespace internal { namespace Static_filters_predicates {
  
  template <typename K_base>
  class Power_test_3:
    public K_base::Power_test_3
  {
    typedef typename K_base::Weighted_point_3 Weighted_point_3;
    typedef typename K_base::FT FT;
    typedef typename K_base::Power_test_3 Base;
  public:
    typedef typename Base::result_type result_type;
  
    #ifndef CGAL_CFG_MATCHING_BUG_6
    using Base::operator();
    #else 
    result_type
    operator()(const Weighted_point_3& p, const Weighted_point_3& q) const
    { 
      return Base::operator()(p,q);
    }
    #endif
    
    
    result_type operator() ( const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & r,
                             const Weighted_point_3 & s,
                             const Weighted_point_3 & t) const
    {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Power_test_3 with 4+1 wpoints", tmp);
      using std::fabs;
      
      double px, py, pz, pwt, qx, qy, qz, qwt, rx, ry, rz, rwt, sx, sy, sz, swt, tx, ty, tz, twt;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)       &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pwt) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)       &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qwt) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)       &&
          fit_in_double(r.z(), rz) && fit_in_double(r.weight(), rwt) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy)       &&
          fit_in_double(s.z(), sz) && fit_in_double(s.weight(), swt) && 
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty)       &&
          fit_in_double(t.z(), tz) && fit_in_double(t.weight(), twt) 
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        using std::fabs;
        
        double dpx = (px - tx);
        double dpy = (py - ty);
        double dpz = (pz - tz);
        double twt_pwt = (twt - pwt);
        double dpt = (((square( dpx ) + square( dpy )) + square( dpz )) + twt_pwt);
        double dqx = (qx - tx);
        double dqy = (qy - ty);
        double dqz = (qz - tz);
        double twt_qwt = (twt - qwt);
        double dqt = (((square( dqx ) + square( dqy )) + square( dqz )) + twt_qwt);
        double drx = (rx - tx);
        double dry = (ry - ty);
        double drz = (rz - tz);
        double twt_rwt = (twt - rwt);
        double drt = (((square( drx ) + square( dry )) + square( drz )) + twt_rwt);
        double dsx = (sx - tx);
        double dsy = (sy - ty);
        double dsz = (sz - tz);
        double twt_swt = (twt - swt);
        double dst = (((square( dsx ) + square( dsy )) + square( dsz )) + twt_swt);
        result_type int_tmp_result;
        double eps;
        double RT_tmp_result = CGAL::determinant( dpx, dpy, dpz, dpt, dqx, dqy, dqz, dqt, drx, dry, drz, drt, dsx, dsy, dsz, dst );
        double max2 = fabs(dpx);
        if( (max2 < fabs(dqx)) ) max2 = fabs(dqx);
        if( (max2 < fabs(drx)) ) max2 = fabs(drx);
        if( (max2 < fabs(dsx)) ) max2 = fabs(dsx);
        double max1 = max2; 
        double max3 = fabs(dpy);
        if( (max3 < fabs(dqy)) ) max3 = fabs(dqy);
        if( (max3 < fabs(dry)) ) max3 = fabs(dry);
        if( (max3 < fabs(dsy)) ) max3 = fabs(dsy);
        if( (max1 < max3) )      max1 = max3;
        double max4 = fabs(dpz);
        if( (max4 < fabs(dqz)) ) max4 = fabs(dqz);
        if( (max4 < fabs(drz)) ) max4 = fabs(drz);
        if( (max4 < fabs(dsz)) ) max4 = fabs(dsz);
        if( (max1 < max4) )      max1 = max4;
        double max5 = fabs(twt_pwt);
        if( (max5 < fabs(twt_qwt)) ) max5 = fabs(twt_qwt);
        if( (max5 < fabs(twt_rwt)) ) max5 = fabs(twt_rwt);
        if( (max5 < fabs(twt_swt)) ) max5 = fabs(twt_swt);
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max2 < lower_bound_1) ) lower_bound_1 = max2;
        if( (max3 < lower_bound_1) ) lower_bound_1 = max3;
        if( (max4 < lower_bound_1) ) lower_bound_1 = max4;
        //handwritten workaround to handle case where all weights are equal
        //if( ((lower_bound_1 < 1.05893684636332247212e-59) || (max5 < 1.12134724458589890404e-118)) )
        if( ((lower_bound_1 < 1.05893684636332247212e-59) || (max5 < 1.12134724458589890404e-118 && max5!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,s,t);
        } 
        else 
        {
            if( ((upper_bound_1 > 3.21387608851797948065e+60) || (max5 > 1.03289995123476274781e+121)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,s,t);
            } 
            eps = (1.67106803095990471147e-13 * (((max2 * max3) * max4) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max5, (max1 * max1) )));
            if( (RT_tmp_result > eps) )
            {
                int_tmp_result = NEGATIVE;
            } 
            else 
            {
                if( (RT_tmp_result < -eps) )
                {
                    int_tmp_result = POSITIVE;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,s,t);
                } 
            } 
        } 
        return int_tmp_result;
      }
      else
        return Base::operator()(p,q,r,s,t);
    }

    result_type operator() ( const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & r,
                             const Weighted_point_3 & t) const
    {

      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Power_test_3 with 3+1 wpoints", tmp);
      using std::fabs;

      double px, py, pz, pwt, qx, qy, qz, qwt, rx, ry, rz, rwt, tx, ty, tz, twt;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)       &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pwt) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)       &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qwt) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)       &&
          fit_in_double(r.z(), rz) && fit_in_double(r.weight(), rwt) &&
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty)       &&
          fit_in_double(t.z(), tz) && fit_in_double(t.weight(), twt) 
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        using std::fabs;
        
        double dpx = (px - tx);
        double dpy = (py - ty);
        double dpz = (pz - tz);
        double twt_pwt = (twt - pwt);
        double dpt = (((square( dpx ) + square( dpy )) + square( dpz )) + twt_pwt);
        double dqx = (qx - tx);
        double dqy = (qy - ty);
        double dqz = (qz - tz);
        double twt_qwt = (twt - qwt);
        double dqt = (((square( dqx ) + square( dqy )) + square( dqz )) + twt_qwt);
        double drx = (rx - tx);
        double dry = (ry - ty);
        double drz = (rz - tz);
        double twt_rwt = (twt - rwt);
        double drt = (((square( drx ) + square( dry )) + square( drz )) + twt_rwt);
        int cmp;
        int int_tmp_result;
        double eps;
        double RT_tmp_result = CGAL::determinant( dpx, dpy, dpt, dqx, dqy, dqt, drx, dry, drt );
        double max7 = fabs(dpz);
        if( (max7 < fabs(dqz)) ) max7 = fabs(dqz);
        if( (max7 < fabs(drz)) ) max7 = fabs(drz);
        double max1 = max7;
        double max2 = fabs(dpx);
        if( (max2 < fabs(dqx)) ) max2 = fabs(dqx);
        if( (max2 < fabs(drx)) ) max2 = fabs(drx);
        if( (max1 < max2) )      max1 = max2;
        double max3 = fabs(dpy);
        if( (max3 < fabs(dqy)) ) max3 = fabs(dqy);
        if( (max3 < fabs(dry)) ) max3 = fabs(dry);
        if( (max1 < max3) )      max1 = max3;
        double max4 = fabs(twt_pwt);
        if( (max4 < fabs(twt_qwt)) ) max4 = fabs(twt_qwt);
        if( (max4 < fabs(twt_rwt)) ) max4 = fabs(twt_rwt);
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max2 < lower_bound_1) ) lower_bound_1 = max2;
        if( (max3 < lower_bound_1) ) lower_bound_1 = max3;
        //handwritten workaround to handle case where all weights are equal        
        //if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148)) )
        if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148 && max4!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,t);
        } 
        else 
        {
            if( ((upper_bound_1 > 7.23700557733225980357e+75) || (max4 > 5.23742497263382350320e+151)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            eps = (3.04426660386257731823e-14 * ((max2 * max3) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max4, (max1 * max1) )));
            if( (RT_tmp_result > eps) )
            {
                int_tmp_result = 1;
            } 
            else 
            {
                if( (RT_tmp_result < -eps) )
                {
                    int_tmp_result = -1;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
            } 
        } 
        cmp = int_tmp_result;
        double px_rx = (px - rx);
        double qy_ry = (qy - ry);
        double qx_rx = (qx - rx);
        double py_ry = (py - ry);
        double max5 = fabs(px_rx);
        if( (max5 < fabs(qx_rx)) ) max5 = fabs(qx_rx);
        double max6 = fabs(qy_ry);
        if( (max6 < fabs(py_ry)) ) max6 = fabs(py_ry);
        if( (cmp != 0) )
        {
            int int_tmp_result_FFWKCAA;
            double double_tmp_result;
            double_tmp_result = ((px_rx * qy_ry) - (qx_rx * py_ry));
            lower_bound_1 = max5;
            upper_bound_1 = max5;
            if( (max6 < lower_bound_1) ) lower_bound_1 = max6;
            else 
            {
                if( (max6 > upper_bound_1) ) upper_bound_1 = max6;
            } 
            if( (lower_bound_1 < 5.00368081960964690982e-147) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            else 
            {
                if( (upper_bound_1 > 7.23700557733225980357e+75) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
                eps = (8.88720573725927976811e-16 * (max5 * max6));
                if( (double_tmp_result > eps) ) int_tmp_result_FFWKCAA = 1;
                else 
                {
                    if( (double_tmp_result < -eps) )
                    {
                        int_tmp_result_FFWKCAA = -1;
                    } 
                    else 
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return Base::operator()(p,q,r,t);
                    } 
                } 
            } 
            return static_cast<result_type>(cmp * int_tmp_result_FFWKCAA);
        } 
        int int_tmp_result_k60Ocge;
        double RT_tmp_result_3SPBwDj = CGAL::determinant( dpx, dpz, dpt, dqx, dqz, dqt, drx, drz, drt );
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max7 < lower_bound_1) ) lower_bound_1 = max7;
        if( (max2 < lower_bound_1) ) lower_bound_1 = max2;
        //handwritten workaround to handle case where all weights are equal
        //if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148)) )
        if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148 && max4!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,t);
        } 
        else 
        {
            if( ((upper_bound_1 > 7.23700557733225980357e+75) || (max4 > 5.23742497263382350320e+151)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            eps = (3.04426660386257731823e-14 * ((max2 * max7) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max4, (max1 * max1) )));
            if( (RT_tmp_result_3SPBwDj > eps) )
            {
                int_tmp_result_k60Ocge = 1;
            } 
            else 
            {
                if( (RT_tmp_result_3SPBwDj < -eps) )
                {
                    int_tmp_result_k60Ocge = -1;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
            } 
        } 
        cmp = int_tmp_result_k60Ocge;
        double qz_rz = (qz - rz);
        double pz_rz = (pz - rz);
        double max8 = fabs(qz_rz);
        if( (max8 < fabs(pz_rz)) ) max8 = fabs(pz_rz);
        if( (cmp != 0) )
        {
            int int_tmp_result_k3Lzf6g;
            double double_tmp_result_Gx4H;
            double_tmp_result_Gx4H = ((px_rx * qz_rz) - (qx_rx * pz_rz));
            lower_bound_1 = max5;
            upper_bound_1 = max5;
            if( (max8 < lower_bound_1) ) lower_bound_1 = max8;
            else 
            {
                if( (max8 > upper_bound_1) ) upper_bound_1 = max8;
            } 
            if( (lower_bound_1 < 5.00368081960964690982e-147) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            else 
            {
                if( (upper_bound_1 > 7.23700557733225980357e+75) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
                eps = (8.88720573725927976811e-16 * (max5 * max8));
                if( (double_tmp_result_Gx4H > eps) )
                {
                    int_tmp_result_k3Lzf6g = 1;
                } 
                else 
                {
                    if( (double_tmp_result_Gx4H < -eps) )
                    {
                        int_tmp_result_k3Lzf6g = -1;
                    } 
                    else 
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return Base::operator()(p,q,r,t);
                    } 
                } 
            } 
            return static_cast<result_type>(cmp * int_tmp_result_k3Lzf6g);
        } 
        int int_tmp_result_AvrrXBP;
        double RT_tmp_result_feLwnHn = CGAL::determinant( dpy, dpz, dpt, dqy, dqz, dqt, dry, drz, drt );
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max7 < lower_bound_1) ) lower_bound_1 = max7;
        if( (max3 < lower_bound_1) ) lower_bound_1 = max3;
        //handwritten workaround to handle case where all weights are equal
        //if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148)) )
        if( ((lower_bound_1 < 2.92391967062015793913e-74) || (max4 < 8.54930624023949352313e-148 && max4!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,t);
        } 
        else 
        {
            if( ((upper_bound_1 > 7.23700557733225980357e+75) || (max4 > 5.23742497263382350320e+151)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            eps = (3.04426660386257731823e-14 * ((max3 * max7) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max4, (max1 * max1) )));
            if( (RT_tmp_result_feLwnHn > eps) )
            {
                int_tmp_result_AvrrXBP = 1;
            } 
            else 
            {
                if( (RT_tmp_result_feLwnHn < -eps) )
                {
                    int_tmp_result_AvrrXBP = -1;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
            } 
        } 
        cmp = int_tmp_result_AvrrXBP;
        int int_tmp_result_agX3WsT;
        double double_tmp_result_Dw20Kqh = ((py_ry * qz_rz) - (qy_ry * pz_rz));
        lower_bound_1 = max8;
        upper_bound_1 = max8;
        if( (max6 < lower_bound_1) )
        {
            lower_bound_1 = max6;
        } 
        else 
        {
            if( (max6 > upper_bound_1) ) upper_bound_1 = max6;
        } 
        if( (lower_bound_1 < 5.00368081960964690982e-147) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,t);
        } 
        else 
        {
            if( (upper_bound_1 > 7.23700557733225980357e+75) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,t);
            } 
            eps = (8.88720573725927976811e-16 * (max6 * max8));
            if( (double_tmp_result_Dw20Kqh > eps) )
            {
                int_tmp_result_agX3WsT = 1;
            } 
            else 
            {
                if( (double_tmp_result_Dw20Kqh < -eps) )
                {
                    int_tmp_result_agX3WsT = -1;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,r,t);
                } 
            } 
        } 
        return static_cast<result_type>(cmp * int_tmp_result_agX3WsT);        
      }
      else
        return Base::operator()(p,q,r,t);
    }

    result_type operator() ( const Weighted_point_3 & p,
                             const Weighted_point_3 & q,
                             const Weighted_point_3 & t) const
    {

      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Power_test_3 with 2+1 wpoints", tmp);
      using std::fabs;
      
      double px, py, pz, pwt, qx, qy, qz, qwt, tx, ty, tz, twt;
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)       &&
          fit_in_double(p.z(), pz) && fit_in_double(p.weight(), pwt) && 
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)       &&
          fit_in_double(q.z(), qz) && fit_in_double(q.weight(), qwt) &&
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty)       &&
          fit_in_double(t.z(), tz) && fit_in_double(t.weight(), twt) 
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
          
        double dpx = (px - tx);
        double dpy = (py - ty);
        double dpz = (pz - tz);
        double twt_pwt = (twt - pwt);
        double dpt = (((square( dpx ) + square( dpy )) + square( dpz )) + twt_pwt);
        double dqx = (qx - tx);
        double dqy = (qy - ty);
        double dqz = (qz - tz);
        double twt_qwt = (twt - qwt);
        double dqt = (((square( dqx ) + square( dqy )) + square( dqz )) + twt_qwt);
        int cmp = ((px > qx) ? 1 : ((px < qx) ? -1 : 0));
        double eps;
        double max1;
        double max4 = fabs(dpy);
        if( (max4 < fabs(dqy)) ) max4 = fabs(dqy);
        max1 = max4;
        double max5 = fabs(dpz); 
        if( (max5 < fabs(dqz)) ) max5 = fabs(dqz);
        if( (max1 < max5) ) max1 = max5;
        double max2 = fabs(dpx);
        if( (max2 < fabs(dqx)) ) max2 = fabs(dqx);
        if( (max1 < max2) ) max1 = max2;
        double max3 = fabs(twt_pwt);
        if( (max3 < fabs(twt_qwt)) ) max3 = fabs(twt_qwt);
        double lower_bound_1;
        double upper_bound_1;
        if( (cmp != 0) )
        {
            int int_tmp_result;
            double double_tmp_result;
            double_tmp_result = ((dpx * dqt) - (dqx * dpt));
            lower_bound_1 = max1;
            upper_bound_1 = max1;
            if( (max2 < lower_bound_1) ) lower_bound_1 = max2;
            //handwritten workaround to handle case where all weights are equal
            //if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195)) )
            if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195 && max3!=0)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,t);
            } 
            else 
            {
                if( ((upper_bound_1 > 5.59936185544450928309e+101) || (max3 > 3.13528531882069776730e+203)) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,t);
                } 
                eps = (6.88858782307641768480e-15 * (max2 * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max3, (max1 * max1) )));
                if( (double_tmp_result > eps) )
                {
                    int_tmp_result = 1;
                } 
                else 
                {
                    if( (double_tmp_result < -eps) )
                    {
                        int_tmp_result = -1;
                    } 
                    else 
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return Base::operator()(p,q,t);
                    } 
                } 
            } 
            return static_cast<result_type>(cmp * int_tmp_result);
        } 
        cmp = ((py > qy) ? 1 : ((py < qy) ? -1 : 0));
        if( (cmp != 0) )
        {
          int int_tmp_result_FFWKCAA;
          double double_tmp_result_k60Ocge;
          double_tmp_result_k60Ocge = ((dpy * dqt) - (dqy * dpt));
          lower_bound_1 = max1;
          upper_bound_1 = max1;
          if( (max4 < lower_bound_1) ) lower_bound_1 = max4;
          //handwritten workaround to handle case where all weights are equal
          //if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195)) )
          if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195 && max3!=0)) )
          {
            CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
            return Base::operator()(p,q,t);
          } 
          else 
          {
              if( ((upper_bound_1 > 5.59936185544450928309e+101) || (max3 > 3.13528531882069776730e+203)) )
              {
                CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                return Base::operator()(p,q,t);
              } 
              eps = (6.88858782307641768480e-15 * (max4 * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max3, (max1 * max1) )));
              if( (double_tmp_result_k60Ocge > eps) )
              {
                  int_tmp_result_FFWKCAA = 1;
              } 
              else 
              {
                  if( (double_tmp_result_k60Ocge < -eps) )
                  {
                      int_tmp_result_FFWKCAA = -1;
                  } 
                  else 
                  {
                    CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                    return Base::operator()(p,q,t);
                  } 
              } 
          } 
          return static_cast<result_type>(cmp * int_tmp_result_FFWKCAA);
        } 
        cmp = ((pz > qz) ? 1 : ((pz < qz) ? -1 : 0));
        int int_tmp_result_3SPBwDj;
        double double_tmp_result_k3Lzf6g = ((dpz * dqt) - (dqz * dpt));
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max5 < lower_bound_1) )
        {
            lower_bound_1 = max5;
        } 
        //handwritten workaround to handle case where all weights are equal
        //if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195)) )
        if( ((lower_bound_1 < 4.89808663633813414271e-98) || (max3 < 2.39912526970742181620e-195 && max3!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,t);
        } 
        else 
        {
            if( ((upper_bound_1 > 5.59936185544450928309e+101) || (max3 > 3.13528531882069776730e+203)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,t);
            } 
            eps = (6.88858782307641768480e-15 * (max5 * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION( max3, (max1 * max1) )));
            if( (double_tmp_result_k3Lzf6g > eps) )
            {
                int_tmp_result_3SPBwDj = 1;
            } 
            else 
            {
                if( (double_tmp_result_k3Lzf6g < -eps) )
                {
                    int_tmp_result_3SPBwDj = -1;
                } 
                else 
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return Base::operator()(p,q,t);
                } 
            } 
        } 
        return static_cast<result_type>(cmp * int_tmp_result_3SPBwDj);           
      }
      else
        return Base::operator()(p,q,t);
    }
  
  };
  
} } }//namespace CGAL::internal::Static_filters_predicates

#endif //CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_3_H
