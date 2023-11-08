// Copyright (c) 2023  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot
// This predicates was generated using the fpg tool written by Andreas Meyer.
//

#ifndef CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/Static_filter_error.h>
#include <CGAL/Filtered_kernel/internal/Static_filters/tools.h>
#include <cmath>


namespace CGAL { namespace internal { namespace Static_filters_predicates {

  template <typename K_base>
  class Power_side_of_oriented_power_circle_2:
    public K_base::Power_side_of_oriented_power_circle_2
  {
    typedef typename K_base::Weighted_point_2 Weighted_point_2;
    typedef typename K_base::FT FT;
    typedef typename K_base::Power_side_of_oriented_power_circle_2 Base;
  public:
    typedef typename Base::result_type result_type;

    using Base::operator();

    result_type operator() ( const Weighted_point_2 & p,
                             const Weighted_point_2 & q,
                             const Weighted_point_2 & r,
                             const Weighted_point_2 & t) const
    {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Power_side_of_power_circle_2 with 3+1 wpoints", tmp);

      double px, py, pwt, qx, qy, qwt, rx, ry, rwt, tx, ty, twt;
      init_double(px, py, pwt, qx, qy, qwt, rx, ry, rwt, (FT*)(0));
      init_double(tx, ty, twt, (FT*)(0));
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)       &&
          fit_in_double(p.weight(), pwt) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)       &&
          fit_in_double(q.weight(), qwt) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)       &&
          fit_in_double(r.weight(), rwt) &&
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty)       &&
          fit_in_double(t.weight(), twt)
          )
        {
          CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
          double dpx = (px - tx);
          double dpy = (py - ty);
          double dpz = (((square( dpx ) + square( dpy )) - pwt) + twt);
          double dqx = (qx - tx);
          double dqy = (qy - ty);
          double dqz = (((square( dqx ) + square( dqy )) - qwt) + twt);
          double drx = (rx - tx);
          double dry = (ry - ty);
          double drz = (((square( drx ) + square( dry )) - rwt) + twt);
          result_type int_tmp_result;
          double RT_tmp_result;
          double eps;
          RT_tmp_result = CGAL::determinant( dpx, dpy, dpz, dqx, dqy, dqz, drx, dry, drz );
          double max1;
          double max3 = CGAL::abs(dpy);
          if( (max3 < CGAL::abs(dqy)) )
            {
              max3 = CGAL::abs(dqy);
            }
          if( (max3 < CGAL::abs(dry)) )
            {
              max3 = CGAL::abs(dry);
            }
          max1 = max3;
          double max2 = CGAL::abs(dpx);
          if( (max2 < CGAL::abs(dqx)) )
            {
              max2 = CGAL::abs(dqx);
            }
          if( (max2 < CGAL::abs(drx)) )
            {
              max2 = CGAL::abs(drx);
            }
          if( (max1 < max2) )
            {
              max1 = max2;
            }
          double max4 = CGAL::abs(pwt);
          if( (max4 < CGAL::abs(qwt)) )
            {
              max4 = CGAL::abs(qwt);
            }
          if( (max4 < CGAL::abs(rwt)) )
            {
              max4 = CGAL::abs(rwt);
            }
          if( (max4 < CGAL::abs(twt)) )
            {
              max4 = CGAL::abs(twt);
            }
          double lower_bound_1;
          double upper_bound_1;
          lower_bound_1 = max3;
          upper_bound_1 = max3;
          if( (max1 < lower_bound_1) )
            {
              lower_bound_1 = max1;
            }
          else
            {
              if( (max1 > upper_bound_1) )
                {
                  upper_bound_1 = max1;
                }
            }
          if( (max2 < lower_bound_1) )
            {
              lower_bound_1 = max2;
            }
          if( ((lower_bound_1 < 2.99168207048872973507e-74) || (max4 < 8.95016161088373414772e-148)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return  Base::operator()(p,q,r,t);
            }
          else
            {
              if( ((upper_bound_1 > 1.44740111546645180002e+76) || (max4 > 2.09496998905352916869e+152)) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return  Base::operator()(p,q,r,t);
                }
              eps = (2.77768297369183927919e-14 * ((max2 * max3) * (CGAL::max)( max4, (max1 * max1) )));
              if( (RT_tmp_result > eps) )
                {
                  int_tmp_result = POSITIVE;
                }
              else
                {
                  if( (RT_tmp_result < -eps) )
                    {
                      int_tmp_result = NEGATIVE;
                    }
                  else
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return Base::operator()(p,q,r,t);
                    }
                }
            }
          return int_tmp_result;

        }
      else
        return Base::operator()(p,q,r,t);
    }


    result_type operator() ( const Weighted_point_2 & p,
                             const Weighted_point_2 & q,
                             const Weighted_point_2 & t) const
    {

      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Power_side_of_oriented_power_circle_2 with 2+1 wpoints", tmp);

      double px, py, pwt, qx, qy, qwt, tx, ty, twt;
      init_double(px, py, pwt, qx, qy, qwt, (FT*)(0));
      init_double(  tx, ty, twt, (FT*)(0));
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)       &&
          fit_in_double(p.weight(), pwt) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)       &&
          fit_in_double(q.weight(), qwt) &&
          fit_in_double(t.x(), tx) && fit_in_double(t.y(), ty)       &&
          fit_in_double(t.weight(), twt)
          )
        {
          CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
          double dpx = (px - tx);
          double dpy = (py - ty);
          double dpz = (((square( dpx ) + square( dpy )) - pwt) + twt);
          double dqx = (qx - tx);
          double dqy = (qy - ty);
          double dqz = (((square( dqx ) + square( dqy )) - qwt) + twt);
          int cmpx;
          cmpx = ((px > qx) ? 1 : ((px < qx) ? -1 : 0));
          double eps;
          double max1;
          double max4 = CGAL::abs(dpy);
          if( (max4 < CGAL::abs(dqy)) )
            {
              max4 = CGAL::abs(dqy);
            }
          max1 = max4;
          double max2 = CGAL::abs(dpx);
          if( (max2 < CGAL::abs(dqx)) )
            {
              max2 = CGAL::abs(dqx);
            }
          if( (max1 < max2) )
            {
              max1 = max2;
            }
          double max3 = CGAL::abs(pwt);
          if( (max3 < CGAL::abs(qwt)) )
            {
              max3 = CGAL::abs(qwt);
            }
          if( (max3 < CGAL::abs(twt)) )
            {
              max3 = CGAL::abs(twt);
            }
          double lower_bound_1;
          double upper_bound_1;
          if( (cmpx != 0) )
            {
              result_type int_tmp_result;
              double RT_tmp_result;
              RT_tmp_result = CGAL::determinant( dpx, dpz, dqx, dqz );
              lower_bound_1 = max2;
              upper_bound_1 = max2;
              if( (max1 < lower_bound_1) )
                {
                  lower_bound_1 = max1;
                }
              else
                {
                  if( (max1 > upper_bound_1) )
                    {
                      upper_bound_1 = max1;
                    }
                }
              if( ((lower_bound_1 < 1.54785988882306922244e-98) || (max3 < 2.39587023542736329316e-196)) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return  Base::operator()(p,q,t);
                }
              else
                {
                  if( ((upper_bound_1 > 5.59936185544450866143e+101) || (max3 > 3.13528531882069741921e+203)) )
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return  Base::operator()(p,q,t);
                    }
                  eps = (5.99997572250729588410e-15 * (max2 * (CGAL::max)( max3, (max1 * max1) )));
                  if( (RT_tmp_result > eps) )
                    {
                      int_tmp_result = POSITIVE;
                    }
                  else
                    {
                      if( (RT_tmp_result < -eps) )
                        {
                          int_tmp_result = NEGATIVE;
                        }
                      else
                        {
                          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                          return  Base::operator()(p,q,t);
                        }
                    }
                }
              return static_cast<result_type>(cmpx * int_tmp_result);
            }
          int cmpy;
          cmpy = ((py > qy) ? 1 : ((py < qy) ? -1 : 0));
          result_type int_tmp_result_FFWKCAA;
          double RT_tmp_result_k60Ocge = CGAL::determinant( dpy, dpz, dqy, dqz );
          lower_bound_1 = max4;
          upper_bound_1 = max4;
          if( (max1 < lower_bound_1) )
            {
              lower_bound_1 = max1;
            }
          else
            {
              if( (max1 > upper_bound_1) )
                {
                  upper_bound_1 = max1;
                }
            }
          if( ((lower_bound_1 < 1.54785988882306922244e-98) || (max3 < 2.39587023542736329316e-196)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return  Base::operator()(p,q,t);
            }
          else
            {
              if( ((upper_bound_1 > 5.59936185544450866143e+101) || (max3 > 3.13528531882069741921e+203)) )
                {
                  CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                  return  Base::operator()(p,q,t);
                }
              eps = (5.99997572250729588410e-15 * (max4 * (CGAL::max)( max3, (max1 * max1) )));
              if( (RT_tmp_result_k60Ocge > eps) )
                {
                  int_tmp_result_FFWKCAA = POSITIVE;
                }
              else
                {
                  if( (RT_tmp_result_k60Ocge < -eps) )
                    {
                      int_tmp_result_FFWKCAA = NEGATIVE;
                    }
                  else
                    {
                      CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
                      return  Base::operator()(p,q,t);
                    }
                }
            }
          return static_cast<result_type>(cmpy * int_tmp_result_FFWKCAA);

        }
      else
        return Base::operator()(p,q,t);
    }


      };

    } } }//namespace CGAL::internal::Static_filters_predicates

#endif //CGAL_INTERNAL_STATIC_FILTERS_POWER_TEST_2_H
