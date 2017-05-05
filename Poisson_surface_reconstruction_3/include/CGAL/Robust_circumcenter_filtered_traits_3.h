// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Stéphane Tayeb
//                 Laurent Saboret

#ifndef CGAL_ROBUST_CIRCUMCENTER_FILTERED_TRAITS_3_H_INCLUDED
#define CGAL_ROBUST_CIRCUMCENTER_FILTERED_TRAITS_3_H_INCLUDED

#include <CGAL/license/Poisson_surface_reconstruction_3.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {

namespace Reconstruction {


template < typename K >
class Robust_filtered_construct_circumcenter_3
{
public:
  typedef typename K::Point_3                        Point_3;
  typedef typename K::FT                             FT;
  typedef typename K::Sphere_3                       Sphere_3;
  typedef Point_3                                    result_type;

  typedef Exact_predicates_exact_constructions_kernel   EK2;
  typedef EK2                                           EK;
  typedef Cartesian_converter<typename K::Kernel, EK2>  To_exact;
  typedef Cartesian_converter<EK2, typename K::Kernel>  Back_from_exact;


  Point_3 operator() ( const Point_3 & p,
                       const Point_3 & q,
                       const Point_3 & r,
                       const Point_3 & s ) const
  {
    typename K::Construct_circumcenter_3 circumcenter =
        K().construct_circumcenter_3_object();
    typename K::Side_of_bounded_sphere_3 side_of_bounded_sphere =
        K().side_of_bounded_sphere_3_object();

    // Compute denominator to swith to exact if it is (close to) 0.
    // TODO: replace hard coded comparison with 1E-13 by static filter.
    const FT denom = compute_denom(p,q,r,s);
    if (denom < -1E-13 || denom > 1E-13)
    {
      Point_3 point = circumcenter(p,q,r,s);

      // Fast output
      if ( side_of_bounded_sphere(p,q,r,s,point) == ON_BOUNDED_SIDE )
        return point;
    }
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_circumcenter_3 exact_circumcenter =
        EK().construct_circumcenter_3_object();

    return back_from_exact(exact_circumcenter(to_exact(p),
                                              to_exact(q),
                                              to_exact(r),
                                              to_exact(s)));
  }

  Point_3 operator() ( const Point_3 & p,
                       const Point_3 & q,
                       const Point_3 & r ) const
  {
    typename K::Construct_circumcenter_3 circumcenter =
      K().construct_circumcenter_3_object();
    typename K::Side_of_bounded_sphere_3 side_of_bounded_sphere =
      K().side_of_bounded_sphere_3_object();

    // Compute denominator to swith to exact if it is (close to) 0.
    // TODO: replace hard coded comparison with 1E-13 by static filter.
    const FT denom = compute_denom(p,q,r);
    if (denom < -1E-13 || denom > 1E-13)
    {
      Point_3 point = circumcenter(p,q,r);

      // Fast output
      if ( side_of_bounded_sphere(p,q,r,point) == ON_BOUNDED_SIDE )
        return point;
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_circumcenter_3 exact_circumcenter =
        EK().construct_circumcenter_3_object();

    return back_from_exact(exact_circumcenter(to_exact(p),
                                              to_exact(q),
                                              to_exact(r)));
  }

  Point_3 operator() ( const Point_3 & p,
                       const Point_3 & q ) const
  {
    typename K::Construct_circumcenter_3 circumcenter =
      K().construct_circumcenter_3_object();
    typename K::Side_of_bounded_sphere_3 side_of_bounded_sphere =
      K().side_of_bounded_sphere_3_object();

    // No division here
    Point_3 point = circumcenter(p,q);

    // Fast output
    if ( side_of_bounded_sphere(p,q,point) == ON_BOUNDED_SIDE )
      return point;

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_circumcenter_3 exact_circumcenter =
      EK().construct_circumcenter_3_object();

    return back_from_exact(exact_circumcenter(to_exact(p),
                                              to_exact(q)));
  }

private:

  FT compute_denom(const Point_3 & p,
                   const Point_3 & q,
                   const Point_3 & r,
                   const Point_3 & s) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z(),
                         s.x(),s.y(),s.z());
  }

  FT compute_denom(const Point_3 & p,
                   const Point_3 & q,
                   const Point_3 & r) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z());
  }

  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &rx, const FT &ry, const FT &rz,
                   const FT &sx, const FT &sy, const FT &sz) const
  {
    const FT qpx = qx-px;
    const FT qpy = qy-py;
    const FT qpz = qz-pz;
    const FT rpx = rx-px;
    const FT rpy = ry-py;
    const FT rpz = rz-pz;
    const FT spx = sx-px;
    const FT spy = sy-py;
    const FT spz = sz-pz;

    return determinant(qpx,qpy,qpz,
                       rpx,rpy,rpz,
                       spx,spy,spz);
  }

  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &rx, const FT &ry, const FT &rz) const
  {
    const FT qpx = qx-px;
    const FT qpy = qy-py;
    const FT qpz = qz-pz;
    const FT rpx = rx-px;
    const FT rpy = ry-py;
    const FT rpz = rz-pz;
    const FT sx = qpy*rpz-qpz*rpy;
    const FT sy = qpz*rpx-qpx*rpz;
    const FT sz = qpx*rpy-qpy*rpx;

    return determinant(qpx,qpy,qpz,
                       rpx,rpy,rpz,
                       sx,sy,sz);
  }

}; // end class Robust_filtered_construct_circumcenter_3


template < typename K >
class Robust_filtered_compute_squared_radius_3
{
public:
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_3;
  typedef typename K::Sphere_3    Sphere_3;
  typedef typename K::Circle_3    Circle_3;
  typedef FT                      result_type;

  typedef Exact_predicates_exact_constructions_kernel   EK2;
  typedef EK2                                           EK;
  typedef Cartesian_converter<typename K::Kernel, EK2>  To_exact;
  typedef Cartesian_converter<EK2, typename K::Kernel>  Back_from_exact;

  FT operator()( const Point_3& p,
                 const Point_3& q,
                 const Point_3& r,
                 const Point_3& s) const
  {
    typename K::Compute_squared_radius_3 squared_radius =
        K().compute_squared_radius_3_object();

    // Compute denominator to swith to exact if it is (close to) 0.
    // TODO: replace hard coded comparison with 1E-13 by static filter.
    const FT denom = compute_denom(p,q,r,s);
    if (denom < -1E-13 || denom > 1E-13)
    {
      // Fast output
      return squared_radius(p,q,r,s);
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Compute_squared_radius_3 exact_squared_radius =
        EK().compute_squared_radius_3_object();

    return back_from_exact(exact_squared_radius(to_exact(p),
                                                to_exact(q),
                                                to_exact(r),
                                                to_exact(s)));
  }

  FT operator()( const Point_3& p,
                 const Point_3& q,
                 const Point_3& r) const
  {
    typename K::Compute_squared_radius_3 squared_radius =
      K().compute_squared_radius_3_object();

    // Compute denominator to swith to exact if it is (close to) 0.
    // TODO: replace hard coded comparison with 1E-13 by static filter.
    const FT denom = compute_denom(p,q,r);
    if (denom < -1E-13 || denom > 1E-13)
    {
      // Fast output
      return squared_radius(p,q,r);
    }

    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Compute_squared_radius_3 exact_squared_radius =
        EK().compute_squared_radius_3_object();

    return back_from_exact(exact_squared_radius(to_exact(p),
                                                to_exact(q),
                                                to_exact(r)));
  }

  FT operator()( const Point_3& p,
                 const Point_3& q) const
  {
    // No division here
    return squared_radiusC3(p.x(), p.y(), p.z(),
                            q.x(), q.y(), q.z());
  }

  result_type
  operator()( const Sphere_3& s) const
  { return s.rep().squared_radius(); }

  result_type
  operator()( const Circle_3& c) const
  { return c.rep().squared_radius(); }

  result_type
  operator()( const Point_3& ) const
  { return FT(0); }

private:

  FT compute_denom(const Point_3 & p,
                   const Point_3 & q,
                   const Point_3 & r,
                   const Point_3 & s) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z(),
                         s.x(),s.y(),s.z());
  }

  FT compute_denom(const Point_3 & p,
                   const Point_3 & q,
                   const Point_3 & r) const
  {
    return compute_denom(p.x(),p.y(),p.z(),
                         q.x(),q.y(),q.z(),
                         r.x(),r.y(),r.z());
  }

  // Compute the denominator computed by squared_radiusC3()
  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &rx, const FT &ry, const FT &rz,
                   const FT &sx, const FT &sy, const FT &sz) const
  {
    // Translate p to origin to simplify the expression.
    FT qpx = qx-px;
    FT qpy = qy-py;
    FT qpz = qz-pz;
    FT rpx = rx-px;
    FT rpy = ry-py;
    FT rpz = rz-pz;
    FT spx = sx-px;
    FT spy = sy-py;
    FT spz = sz-pz;

    return     determinant( qpx,qpy,qpz,
                            rpx,rpy,rpz,
                            spx,spy,spz);
  }

  // Compute the denominator computed by squared_radiusC3()
  FT compute_denom(const FT &px, const FT &py, const FT &pz,
                   const FT &qx, const FT &qy, const FT &qz,
                   const FT &sx, const FT &sy, const FT &sz) const
  {
    // Translate s to origin to simplify the expression.
    FT psx = px-sx;
    FT psy = py-sy;
    FT psz = pz-sz;
    FT qsx = qx-sx;
    FT qsy = qy-sy;
    FT qsz = qz-sz;
    FT rsx = psy*qsz-psz*qsy;
    FT rsy = psz*qsx-psx*qsz;
    FT rsz = psx*qsy-psy*qsx;

    return     determinant( psx,psy,psz,
                            qsx,qsy,qsz,
                            rsx,rsy,rsz);
  }

}; // end class Robust_filtered_compute_squared_radius_3

} //end namespace Reconstruction


/**
 * \internal
 * Robust_circumcenter_filtered_traits_3
 * overrides construct_circumcenter_3_object() and compute_squared_radius_3_object()
 * to get robust ones when called on slivers.
 */
template<class K>
struct Robust_circumcenter_filtered_traits_3
: public K
{
  typedef CGAL::Reconstruction::Robust_filtered_construct_circumcenter_3<K>
                                            Construct_circumcenter_3;

  typedef CGAL::Reconstruction::Robust_filtered_compute_squared_radius_3<K>
                                            Compute_squared_radius_3;

  Construct_circumcenter_3
  construct_circumcenter_3_object() const
  { return Construct_circumcenter_3(); }

  Compute_squared_radius_3
  compute_squared_radius_3_object() const
  { return Compute_squared_radius_3(); }

};  // end class Robust_circumcenter_filtered_traits_3


}  // end namespace CGAL

#endif // CGAL_ROBUST_CIRCUMCENTER_FILTERED_TRAITS_3_H_INCLUDED
