// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Static_filters/Orientation_2.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_STATIC_FILTERS_ORIENTATION_2_H
#define CGAL_STATIC_FILTERS_ORIENTATION_2_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>
// #include <CGAL/Static_filter_error.h> // Only used to precompute constants

CGAL_BEGIN_NAMESPACE

template <class Point>
class SF_Orientation_2
{
  double _static_epsilon;

  // Computes the epsilon for Orientation_2.
  static double ori_2()
  {
    typedef Static_filter_error F;
    F t1 = F(1)-F(1);         // First translation
    F det = det2x2_by_formula(t1, t1,
                              t1, t1); // Full det
    double err = det.error();
    std::cerr << "*** epsilon for Orientation_2 = " << err << std::endl;
    return err;
  }

  static const double epsilon = 3.55271e-15; // ori_2();

protected:

  template < class R >
  friend class Static_filters;

  // These operations are reserved to Static_filters<>, because the context of
  // a predicate is linked to the one of the Static_filter<> it is a member of.
  SF_Orientation_2(const SF_Orientation_2 &) {}

  SF_Orientation_2 & operator=(const SF_Orientation_2 &) {}

  SF_Orientation_2()
  {
      _static_epsilon = HUGE_VAL;
  }

public:
  typedef Orientation result_type;

  void update(double dx, double dy)
  {
      _static_epsilon = dx*dy*epsilon;
  }

  Orientation operator()(const Point &p, const Point &q, const Point &r) const
  {
    return opti_orientationC2(to_double(p.x()), to_double(p.y()),
	                      to_double(q.x()), to_double(q.y()),
	                      to_double(r.x()), to_double(r.y()));
  }

private:

  Orientation
  opti_orientationC2(double px, double py,
                     double qx, double qy,
		     double rx, double ry) const
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("Orientation_2 calls"); ++calls;
#endif

    double pqx = qx-px;
    double pqy = qy-py;
    double prx = rx-px;
    double pry = ry-py;

    double det = det2x2_by_formula(pqx, pqy,
                                   prx, pry);

    // Fully static filter first.
    if (det >  _static_epsilon) return POSITIVE;
    if (det < -_static_epsilon) return NEGATIVE;

#ifdef CGAL_PROFILE
    static Profile_counter st_fail("Orientation_2 static failures"); ++st_fail;
#endif

    // Then semi-static filter.
    double maxx = fabs(px);
    if (maxx < fabs(qx)) maxx = fabs(qx);
    if (maxx < fabs(rx)) maxx = fabs(rx);
    double maxy = fabs(py);
    if (maxy < fabs(qy)) maxy = fabs(qy);
    if (maxy < fabs(ry)) maxy = fabs(ry);
    double eps = epsilon*maxx*maxy;

    if (det > eps)  return POSITIVE;
    if (det < -eps) return NEGATIVE;

#ifdef CGAL_PROFILE
    static Profile_counter fail("Orientation_2 semi-static failures"); ++fail;
#endif

    typedef Simple_cartesian<Filtered_exact<double, MP_Float> >::Point_2 P;

    Orientation oooo = orientation(P(px,py), P(qx,qy), P(rx,ry));
#ifdef CGAL_PROFILE
    if (oooo == ZERO) {
        static Profile_counter det_is_null("Orientation_2 det_is_null");
        ++det_is_null;
    }
#endif
    return oooo;
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_ORIENTATION_2_H
