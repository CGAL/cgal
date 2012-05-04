// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_DISTANCE_PREDICATESH3_H
#define CGAL_DISTANCE_PREDICATESH3_H

namespace CGAL {

template <class R>
Comparison_result
compare_distance_to_point(const PointH3<R>& ,
                          const PointH3<R>& ,
                          const PointH3<R>& );

template <class R>
bool
has_larger_distance_to_point(const PointH3<R>& ,
                             const PointH3<R>& ,
                             const PointH3<R>& );

template <class R>
bool
has_smaller_distance_to_point(const PointH3<R>& ,
                              const PointH3<R>& ,
                              const PointH3<R>& );

template <class R>
Comparison_result
compare_signed_distance_to_plane(const PlaneH3<R>& ,
                                 const PointH3<R>& ,
                                 const PointH3<R>& );

template <class R>
bool
has_larger_signed_distance_to_plane(const PlaneH3<R>& ,
                                    const PointH3<R>& ,
                                    const PointH3<R>& );

template <class R>
bool
has_smaller_signed_distance_to_plane(const PlaneH3<R>&,
                                     const PointH3<R>& ,
                                     const PointH3<R>& );

template <class R>
Comparison_result
compare_signed_distance_to_plane(const PointH3<R>& ,
                                 const PointH3<R>& ,
                                 const PointH3<R>& ,
                                 const PointH3<R>& ,
                                 const PointH3<R>& );

template <class R>
bool
has_larger_signed_distance_to_plane(const PointH3<R>& ,
                                    const PointH3<R>& ,
                                    const PointH3<R>& ,
                                    const PointH3<R>& ,
                                    const PointH3<R>& );

template <class R>
bool
has_smaller_signed_distance_to_plane(const PointH3<R>& ,
                                     const PointH3<R>& ,
                                     const PointH3<R>& ,
                                     const PointH3<R>& ,
                                     const PointH3<R>& );

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_distance_to_point(const PointH3<R>& p,
                          const PointH3<R>& q,
                          const PointH3<R>& r)
{
  typedef typename R::RT RT;

  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT rhx = r.hx();
  const RT rhy = r.hy();
  const RT rhz = r.hz();
  const RT rhw = r.hw();
  const RT RT0 = RT(0);
  const RT RT2 = RT(2);

  RT dosd =   // difference of squared distances

    rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy + qhz*qhz )
                - RT2 * qhw * ( phx*qhx + phy*qhy + phz*qhz )
              )
  - qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy + rhz*rhz )
                - RT2 * rhw * ( phx*rhx + phy*rhy + phz*rhz )
              );

  if ( RT0 < dosd )
  { return LARGER; }
  else
  { return (dosd < RT0) ? SMALLER : EQUAL; }
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_distance_to_point(const PointH3<R>& p,
                             const PointH3<R>& q,
                             const PointH3<R>& r)
{
  typedef typename R::RT RT;

  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT rhx = r.hx();
  const RT rhy = r.hy();
  const RT rhz = r.hz();
  const RT rhw = r.hw();
  const RT RT0 = RT(0);
  const RT RT2 = RT(2);

  RT dosd =   // difference of squared distances

    rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy + qhz*qhz )
                - RT2 * qhw * ( phx*qhx + phy*qhy + phz*qhz )
              )
  - qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy + rhz*rhz )
                - RT2 * rhw * ( phx*rhx + phy*rhy + phz*rhz )
              );

  return ( RT0 < dosd );
}

template < class R>
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_distance_to_point(const PointH3<R>& p,
                              const PointH3<R>& q,
                              const PointH3<R>& r)
{
  typedef typename R::RT RT;

  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT rhx = r.hx();
  const RT rhy = r.hy();
  const RT rhz = r.hz();
  const RT rhw = r.hw();
  const RT RT0 = RT(0);
  const RT RT2 = RT(2);

  RT dosd =   // difference of squared distances

    rhw*rhw * (         phw * ( qhx*qhx + qhy*qhy + qhz*qhz )
                - RT2 * qhw * ( phx*qhx + phy*qhy + phz*qhz )
              )
  - qhw*qhw * (         phw * ( rhx*rhx + rhy*rhy + rhz*rhz )
                - RT2 * rhw * ( phx*rhx + phy*rhy + phz*rhz )
              );

  return ( dosd < RT0 );
}

template < class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_signed_distance_to_plane(const PlaneH3<R>& pl,
                                 const PointH3<R>& p,
                                 const PointH3<R>& q)
{
  typedef typename R::RT RT;

  const RT pla = pl.a();
  const RT plb = pl.b();
  const RT plc = pl.c();
  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT RT0 = RT(0);

  RT scaled_dist_p_minus_scaled_dist_q =
      pla*( phx*qhw - qhx*phw )
    + plb*( phy*qhw - qhy*phw )
    + plc*( phz*qhw - qhz*phw );



  if ( scaled_dist_p_minus_scaled_dist_q < RT0 )
  { return SMALLER; }
  else
  { return (RT0 < scaled_dist_p_minus_scaled_dist_q ) ?  LARGER : EQUAL;}
}

template <class R>
bool
has_larger_signed_distance_to_plane(const PlaneH3<R>& pl,
                                    const PointH3<R>& p,
                                    const PointH3<R>& q )
{
  typedef typename R::RT RT;

  const RT pla = pl.a();
  const RT plb = pl.b();
  const RT plc = pl.c();
  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT RT0 = RT(0);

  RT scaled_dist_p_minus_scaled_dist_q =
      pla*( phx*qhw - qhx*phw )
    + plb*( phy*qhw - qhy*phw )
    + plc*( phz*qhw - qhz*phw );


  return ( RT0 < scaled_dist_p_minus_scaled_dist_q );
}

template <class R>
bool
has_smaller_signed_distance_to_plane(const PlaneH3<R>& pl,
                                     const PointH3<R>& p,
                                     const PointH3<R>& q )
{
  typedef typename R::RT RT;

  const RT pla = pl.a();
  const RT plb = pl.b();
  const RT plc = pl.c();
  const RT phx = p.hx();
  const RT phy = p.hy();
  const RT phz = p.hz();
  const RT phw = p.hw();
  const RT qhx = q.hx();
  const RT qhy = q.hy();
  const RT qhz = q.hz();
  const RT qhw = q.hw();
  const RT RT0 = RT(0);

  RT scaled_dist_p_minus_scaled_dist_q =
      pla*( phx*qhw - qhx*phw )
    + plb*( phy*qhw - qhy*phw )
    + plc*( phz*qhw - qhz*phw );


  return ( scaled_dist_p_minus_scaled_dist_q < RT0 );
}

template <class R>
inline
Comparison_result
compare_signed_distance_to_plane(const PointH3<R>& p,
                                 const PointH3<R>& q,
                                 const PointH3<R>& r,
                                 const PointH3<R>& s,
                                 const PointH3<R>& t)
{
  CGAL_kernel_precondition( !collinear(p,q,r) );
  PlaneH3<R> P(p,q,r);
  return cmp_signed_dist_to_plane( P, s, t);
}

template <class R>
inline
bool
has_larger_signed_distance_to_plane(const PointH3<R>& p,
                                    const PointH3<R>& q,
                                    const PointH3<R>& r,
                                    const PointH3<R>& s,
                                    const PointH3<R>& t)
{ return cmp_signed_dist_to_plane(p,q,r,s,t) == LARGER; }


template <class R>
inline
bool
has_smaller_signed_distance_to_plane(const PointH3<R>& p,
                                     const PointH3<R>& q,
                                     const PointH3<R>& r,
                                     const PointH3<R>& s,
                                     const PointH3<R>& t)
{ return cmp_signed_dist_to_plane(p,q,r,s,t) == SMALLER; }

} //namespace CGAL

#endif //CGAL_DISTANCE_PREDICATESH3_H
