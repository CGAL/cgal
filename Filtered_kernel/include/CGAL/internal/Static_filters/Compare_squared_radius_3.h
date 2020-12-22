// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COMPARE_SQUARED_RADIUS_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COMPARE_SQUARED_RADIUS_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <cmath>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

  template <typename K_base>
  class Compare_squared_radius_3
    : public K_base::Compare_squared_radius_3
  {
    typedef typename K_base::Point_3 Point_3;
    typedef typename K_base::FT FT;
    typedef typename K_base::Compare_squared_radius_3 Base;
  public:
    typedef typename Base::result_type result_type;

    using Base::operator();

    result_type operator() (
        const Point_3& p,
        const Point_3& q,
        const Point_3& r,
        const Point_3& s,
        const FT& w
    ) const {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_squared_radius_3 with 4 points", tmp);

      double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, alpha;
      init_double(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz, alpha, (FT*)(0));

      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) &&
          fit_in_double(r.x(), rx) && fit_in_double(r.y(), ry)      &&
          fit_in_double(r.z(), rz) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy)      &&
          fit_in_double(s.z(), sz) &&
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        double qpx;
        qpx = (qx - px);
        double qpy;
        qpy = (qy - py);
        double qpz;
        qpz = (qz - pz);
        double qp2;
        qp2 = ((square( qpx ) + square( qpy )) + square( qpz ));
        double rpx;
        rpx = (rx - px);
        double rpy;
        rpy = (ry - py);
        double rpz;
        rpz = (rz - pz);
        double rp2;
        rp2 = ((square( rpx ) + square( rpy )) + square( rpz ));
        double spx;
        spx = (sx - px);
        double spy;
        spy = (sy - py);
        double spz;
        spz = (sz - pz);
        double sp2;
        sp2 = ((square( spx ) + square( spy )) + square( spz ));
        double num_x;
        num_x = CGAL::determinant( qpy, qpz, qp2, rpy, rpz, rp2, spy, spz, sp2 );
        double num_y;
        num_y = CGAL::determinant( qpx, qpz, qp2, rpx, rpz, rp2, spx, spz, sp2 );
        double num_z;
        num_z = CGAL::determinant( qpx, qpy, qp2, rpx, rpy, rp2, spx, spy, sp2 );
        double den;
        den = CGAL::determinant( qpx, qpy, qpz, rpx, rpy, rpz, spx, spy, spz );
        int int_tmp_result;
        double double_tmp_result;
        double eps;
        double_tmp_result = (((square( num_x ) + square( num_y )) + square( num_z )) - ((alpha * 4.00000000000000000000e+00) * square( den )));
        double max1;
        double aqpx = CGAL::abs(qpx);
        double aqpz = CGAL::abs(qpz);

        double arpx = CGAL::abs(rpx);
        double arpy = CGAL::abs(rpy);
        double arpz = CGAL::abs(rpz);

        double aspx = CGAL::abs(spx);
        double aspy = CGAL::abs(spy);
        double aspz = CGAL::abs(spz);

        double max2 = CGAL::abs(qpy);
        if( (max2 < aqpz) )
        {
            max2 = aqpz;
        }
        if( (max2 < arpy) )
        {
            max2 = arpy;
        }
        if( (max2 < arpz) )
        {
            max2 = arpz;
        }
        if( (max2 < aspy) )
        {
            max2 = aspy;
        }
        if( (max2 < aspz) )
        {
            max2 = aspz;
        }
        max1 = max2;
        if( (max1 < aqpx) )
        {
            max1 = aqpx;
        }
        if( (max1 < arpx) )
        {
            max1 = arpx;
        }
        if( (max1 < aspx) )
        {
            max1 = aspx;
        }
        double max3 = CGAL::abs(alpha);
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max2 < lower_bound_1) )
        {
            lower_bound_1 = max2;
        }
        //handwritten workaround to handle case of alpha=0 (variable alone in its group)
        //if( ((lower_bound_1 < 1.00913582207214915294e-37) || (max3 < 1.01835510738923227819e-74)) )
        if( ((lower_bound_1 < 1.00913582207214915294e-37) || (max3 < 1.01835510738923227819e-74 && max3!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,r,s,w);
        }
        else
        {
            if( ((upper_bound_1 > 2.59614842926741294957e+33) || (max3 > 6.73998666678765545893e+66)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,r,s,w);
            }
            eps = (2.92119946853791784711e-12 * ((((((max1 * max1) * max1) * max1) * max1) * max2) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max3, (max1 * max1) )));
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
                  return Base::operator()(p,q,r,s,w);
                }
            }
        }
        return static_cast<Sign>(int_tmp_result);
      }
      else
        return Base::operator()(p,q,r,s,w);
    }

    result_type operator() (
        const Point_3& p,
        const Point_3& q,
        const Point_3& s,
        const FT& w
    ) const {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_squared_radius_3 with 3 points", tmp);

      double px, py, pz, qx, qy, qz, sx, sy, sz, alpha;
      init_double(px, py, pz, qx, qy, qz, sx, sy, sz, alpha, (FT*)(0));

      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) &&
          fit_in_double(s.x(), sx) && fit_in_double(s.y(), sy)      &&
          fit_in_double(s.z(), sz) &&
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        double psx;
        psx = (px - sx);
        double psy;
        psy = (py - sy);
        double psz;
        psz = (pz - sz);
        double ps2;
        ps2 = ((square( psx ) + square( psy )) + square( psz ));
        double qsx;
        qsx = (qx - sx);
        double qsy;
        qsy = (qy - sy);
        double qsz;
        qsz = (qz - sz);
        double qs2;
        qs2 = ((square( qsx ) + square( qsy )) + square( qsz ));
        double rsx;
        rsx = ((psy * qsz) - (psz * qsy));
        double rsy;
        rsy = ((psz * qsx) - (psx * qsz));
        double rsz;
        rsz = ((psx * qsy) - (psy * qsx));
        double num_x;
        num_x = ((ps2 * CGAL::determinant( qsy, qsz, rsy, rsz )) - (qs2 * CGAL::determinant( psy, psz, rsy, rsz )));
        double num_y;
        num_y = ((ps2 * CGAL::determinant( qsx, qsz, rsx, rsz )) - (qs2 * CGAL::determinant( psx, psz, rsx, rsz )));
        double num_z;
        num_z = ((ps2 * CGAL::determinant( qsx, qsy, rsx, rsy )) - (qs2 * CGAL::determinant( psx, psy, rsx, rsy )));
        double den;
        den = CGAL::determinant( psx, psy, psz, qsx, qsy, qsz, rsx, rsy, rsz );
        int int_tmp_result;
        double double_tmp_result;
        double eps;
        double_tmp_result = (((square( num_x ) + square( num_y )) + square( num_z )) - ((alpha * 4.00000000000000000000e+00) * square( den )));
        double max1;
        double max2 = CGAL::abs(psx);
        double apsy = CGAL::abs(psy);
        double apsz = CGAL::abs(psz);

        double aqsx = CGAL::abs(qsx);
        double aqsy = CGAL::abs(qsy);
        double aqsz = CGAL::abs(qsz);

        if( (max2 < apsy) )
        {
            max2 = apsy;
        }
        if( (max2 < aqsx) )
        {
            max2 = aqsx;
        }
        if( (max2 < aqsy) )
        {
            max2 = aqsy;
        }
        max1 = max2;
        if( (max1 < apsz) )
        {
            max1 = apsz;
        }
        if( (max1 < aqsz) )
        {
            max1 = aqsz;
        }
        double max3 = CGAL::abs(alpha);
        double lower_bound_1;
        double upper_bound_1;
        lower_bound_1 = max1;
        upper_bound_1 = max1;
        if( (max2 < lower_bound_1) )
        {
            lower_bound_1 = max2;
        }

        //handwritten workaround to handle case of alpha=0 (variable alone in its group)
        //if( ((lower_bound_1 < 2.26156385701827020260e-30) || (max3 < 5.11467107937135531427e-60)) )
        if( ((lower_bound_1 < 2.26156385701827020260e-30) || (max3 < 5.11467107937135531427e-60 && max3!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,s,w);
        }
        else
        {
            if( ((upper_bound_1 > 1.23794003928538000002e+27) || (max3 > 1.53249554086588817779e+54)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,s,w);
            }
            eps = (6.35705373458387935514e-12 * ((((((((max1 * max1) * max2) * max1) * max1) * max1) * max1) * max1) * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max3, (max1 * max1) )));
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
                  return Base::operator()(p,q,s,w);
                }
            }
        }
        return static_cast<result_type>(int_tmp_result);
      }
      else
        return Base::operator()(p,q,s,w);
    }


    result_type operator() (
        const Point_3& p,
        const Point_3& q,
        const FT& w
    ) const {
      CGAL_BRANCH_PROFILER_3("semi-static failures/attempts/calls to   : Compare_squared_radius_3 with 2 points", tmp);

      double px, py, pz, qx, qy, qz, alpha;
      init_double(px, py, pz, qx, qy, qz, alpha, (FT*)(0));
      if( fit_in_double(p.x(), px) && fit_in_double(p.y(), py)      &&
          fit_in_double(p.z(), pz) &&
          fit_in_double(q.x(), qx) && fit_in_double(q.y(), qy)      &&
          fit_in_double(q.z(), qz) &&
          fit_in_double(w, alpha)
        )
      {
        CGAL_BRANCH_PROFILER_BRANCH_1(tmp);
        double px_qx = (px - qx);
        double py_qy = (py - qy);
        double pz_qz = (pz - qz);
        int int_tmp_result;
        double double_tmp_result;
        double eps;
        double_tmp_result = (((square( px_qx ) + square( py_qy )) + square( pz_qz )) - (alpha * 4.00000000000000000000e+00));
        double max1 = CGAL::abs(px_qx);
        double apy_qy = CGAL::abs(py_qy);
        double apz_qz = CGAL::abs(pz_qz);
        if( (max1 < apy_qy) )
        {
            max1 = apy_qy;
        }
        if( (max1 < apz_qz) )
        {
            max1 = apz_qz;
        }
        double max2 = CGAL::abs(alpha);

        //handwritten workaround to handle case of alpha=0 (variable alone in its group)
        //if( ((max1 < 8.85464260923320109378e-147) || (max2 < 7.84046957372481590760e-293)) )
        if( ((max1 < 8.85464260923320109378e-147) || (max2 < 7.84046957372481590760e-293 && max2!=0)) )
        {
          CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
          return Base::operator()(p,q,w);
        }
        else
        {
            if( ((max1 > 2.23974474217780371323e+102) || (max2 > 5.01645651011311642768e+204)) )
            {
              CGAL_BRANCH_PROFILER_BRANCH_2(tmp);
              return Base::operator()(p,q,w);
            }
            eps = (2.11094186805729591487e-15 * CGAL::max BOOST_PREVENT_MACRO_SUBSTITUTION ( max2, (max1 * max1) ));
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
                  return Base::operator()(p,q,w);
                }
            }
        }
        return static_cast<Sign>(int_tmp_result);
      }
      else
        return Base::operator()(p,q,w);
    }

  };

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COMPARE_SQUARED_RADIUS_3_H
