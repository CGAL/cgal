// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 2000, October 15
// 
// source        : distance_predicatesH3.fw
// file          : distance_predicatesH3.h
// package       : H3 (2.14)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.14
// revision_date : 15 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_DISTANCE_PREDICATESH3_H
#define CGAL_DISTANCE_PREDICATESH3_H

CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
Comparison_result
cmp_dist_to_point(const PointH3<FT,RT>& ,
                  const PointH3<FT,RT>& ,
                  const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_larger_dist_to_point(const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_smaller_dist_to_point(const PointH3<FT,RT>& ,
                          const PointH3<FT,RT>& ,
                          const PointH3<FT,RT>& );

template <class FT, class RT>
Comparison_result
cmp_signed_dist_to_plane(const PlaneH3<FT,RT>& ,
                         const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_larger_signed_dist_to_plane(const PlaneH3<FT,RT>& ,
                                const PointH3<FT,RT>& ,
                                const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_smaller_signed_dist_to_plane(const PlaneH3<FT,RT>&,
                                 const PointH3<FT,RT>& ,
                                 const PointH3<FT,RT>& );

template <class FT, class RT>
Comparison_result
cmp_signed_dist_to_plane(const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& ,
                         const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_larger_signed_dist_to_plane(const PointH3<FT,RT>& ,
                                const PointH3<FT,RT>& ,
                                const PointH3<FT,RT>& ,
                                const PointH3<FT,RT>& ,
                                const PointH3<FT,RT>& );

template <class FT, class RT>
bool
has_smaller_signed_dist_to_plane(const PointH3<FT,RT>& ,
                                 const PointH3<FT,RT>& ,
                                 const PointH3<FT,RT>& ,
                                 const PointH3<FT,RT>& ,
                                 const PointH3<FT,RT>& );

template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
cmp_dist_to_point(const PointH3<FT,RT>& p,
                  const PointH3<FT,RT>& q,
                  const PointH3<FT,RT>& r)
{
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

template < class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_dist_to_point(const PointH3<FT,RT>& p,
                         const PointH3<FT,RT>& q,
                         const PointH3<FT,RT>& r)
{
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

template < class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_dist_to_point(const PointH3<FT,RT>& p,
                          const PointH3<FT,RT>& q,
                          const PointH3<FT,RT>& r)
{
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

template < class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
cmp_signed_dist_to_plane(const PlaneH3<FT,RT>& pl,
                         const PointH3<FT,RT>& p,
                         const PointH3<FT,RT>& q)
{
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

template <class FT, class RT>
bool
has_larger_signed_dist_to_plane(const PlaneH3<FT,RT>& pl,
                                const PointH3<FT,RT>& p,
                                const PointH3<FT,RT>& q )
{
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

template <class FT, class RT>
bool
has_smaller_signed_dist_to_plane(const PlaneH3<FT,RT>& pl,
                                 const PointH3<FT,RT>& p,
                                 const PointH3<FT,RT>& q )
{
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

template <class FT, class RT>
Comparison_result
cmp_signed_dist_to_plane(const PointH3<FT,RT>& p,
                         const PointH3<FT,RT>& q,
                         const PointH3<FT,RT>& r,
                         const PointH3<FT,RT>& s,
                         const PointH3<FT,RT>& t)
{
  CGAL_kernel_precondition( !collinear(p,q,r) );
  PlaneH3<FT,RT> P(p,q,r);
  return cmp_signed_dist_to_plane( P, s, t);
}

template <class FT, class RT>
bool
has_larger_signed_dist_to_plane(const PointH3<FT,RT>& p,
                                const PointH3<FT,RT>& q,
                                const PointH3<FT,RT>& r,
                                const PointH3<FT,RT>& s,
                                const PointH3<FT,RT>& t)
{ return cmp_signed_dist_to_plane(p,q,r,s,t) == LARGER; }


template <class FT, class RT>
bool
has_smaller_signed_dist_to_plane(const PointH3<FT,RT>& p,
                                 const PointH3<FT,RT>& q,
                                 const PointH3<FT,RT>& r,
                                 const PointH3<FT,RT>& s,
                                 const PointH3<FT,RT>& t)
{ return cmp_signed_dist_to_plane(p,q,r,s,t) == SMALLER; }


CGAL_END_NAMESPACE


#endif //CGAL_DISTANCE_PREDICATESH3_H
