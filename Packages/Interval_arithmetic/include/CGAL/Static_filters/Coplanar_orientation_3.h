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
// file          : include/CGAL/Static_filters/Coplanar_orientation_3.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H
#define CGAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Profile_counter.h>

CGAL_BEGIN_NAMESPACE

// This one is easy : it's just 3 Orientation_2.
// Note that it will mix the profiler stats.

template <class Point, class Orientation_2>
class SF_Coplanar_orientation_3
{
  Orientation_2 oxy, oyz, oxz;

protected:

  template < class R >
  friend class Static_filters;

  // These operations are reserved to Static_filters<>, because the context of
  // a predicate is linked to the one of the Static_filter<> it is a member of.
  SF_Coplanar_orientation_3(const SF_Coplanar_orientation_3 &s)
      : oxy(s.oxy), oyz(s.oyz), oxz(s.oxz) {}

  SF_Coplanar_orientation_3 & operator=(const SF_Coplanar_orientation_3 &s)
  {
      oxy = s.oxy;
      oyz = s.oyz;
      oxz = s.oxz;
      return *this;
  }

  SF_Coplanar_orientation_3() {}

public:
  typedef Orientation result_type;

  void update(double dx, double dy, double dz)
  {
      oxy.update(dx, dy);
      oyz.update(dy, dz);
      oxz.update(dx, dz);
  }

  Orientation operator()(const Point &p, const Point &q, const Point &r) const
  {
      return opti_coplanar_orientationC3(
	    to_double(p.x()), to_double(p.y()), to_double(p.z()),
	    to_double(q.x()), to_double(q.y()), to_double(q.z()),
	    to_double(r.x()), to_double(r.y()), to_double(r.z()));
  }

  Orientation operator()(const Point &p, const Point &q,
                         const Point &r, const Point &s) const
  {
      return opti_coplanar_orientationC3(
	    to_double(p.x()), to_double(p.y()), to_double(p.z()),
	    to_double(q.x()), to_double(q.y()), to_double(q.z()),
	    to_double(r.x()), to_double(r.y()), to_double(r.z()),
	    to_double(s.x()), to_double(s.y()), to_double(s.z()));
  }

private:
  Orientation
  opti_coplanar_orientationC3(double px, double py, double pz,
                              double qx, double qy, double qz,
		              double rx, double ry, double rz) const
  {
      CGAL_PROFILER(calls, "Coplanar_orientation_3 #1 calls")

      Orientation oxy_pqr = oxy.opti_orientationC2(px,py,qx,qy,rx,ry);
      if (oxy_pqr != COLLINEAR)
          return oxy_pqr;

      CGAL_PROFILER(step2, "Coplanar_orientation_3 #1 step2")

      Orientation oyz_pqr = oyz.opti_orientationC2(py,pz,qy,qz,ry,rz);
      if (oyz_pqr != COLLINEAR)
          return oyz_pqr;

      CGAL_PROFILER(step3, "Coplanar_orientation_3 #1 step3")

      return oxz.opti_orientationC2(px,pz,qx,qz,rx,rz);
  }

  Orientation
  opti_coplanar_orientationC3(double px, double py, double pz,
                              double qx, double qy, double qz,
		              double rx, double ry, double rz,
		              double sx, double sy, double sz) const
  {
      CGAL_PROFILER(calls, "Coplanar_orientation_3 #2 calls")

      Orientation oxy_pqr = oxy.opti_orientationC2(px,py,qx,qy,rx,ry);
      if (oxy_pqr != COLLINEAR)
          return Orientation( oxy_pqr *
		              oxy.opti_orientationC2(px,py,qx,qy,sx,sy));

      CGAL_PROFILER(step2, "Coplanar_orientation_3 #2 step2")

      Orientation oyz_pqr = oyz.opti_orientationC2(py,pz,qy,qz,ry,rz);
      if (oyz_pqr != COLLINEAR)
          return Orientation( oyz_pqr *
		              oyz.opti_orientationC2(py,pz,qy,qz,sy,sz));

      CGAL_PROFILER(step3, "Coplanar_orientation_3 #2 step3")

      Orientation oxz_pqr = oxz.opti_orientationC2(px,pz,qx,qz,rx,rz);
      CGAL_kernel_assertion(oxz_pqr != COLLINEAR);
      return Orientation( oxz_pqr * oxz.opti_orientationC2(px,pz,qx,qz,sx,sz));
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_COPLANAR_ORIENTATION_3_H
