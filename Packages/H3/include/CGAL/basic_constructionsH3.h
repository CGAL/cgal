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
// release_date  : 
// 
// file          : basic_constructionsH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_BASIC_CONSTRUCTIONSH3_H
#define CGAL_BASIC_CONSTRUCTIONSH3_H

CGAL_BEGIN_NAMESPACE

template <class R>
PointH3<R>
_projection(const PointH3<R>& p, const PlaneH3<R>& pl)
{
  typedef typename R::RT RT;
  if ( pl.has_on(p) ) return p;

  RT A = pl.a();
  RT B = pl.b();
  RT C = pl.c();
  RT D = pl.d();
  RT phx = p.hx();
  RT phy = p.hy();
  RT phz = p.hz();
  RT phw = p.hw();

  RT num = A * phx  +  B * phy  +  C * phz  +  D * phw;
  RT den = A * A    +  B * B    +  C * C;

  return PointH3<R>( num * A  -  den * phx,
                              num * B  -  den * phy,
                              num * C  -  den * phz,
                             -den );
}

template <class R>
PointH3<R>
midpoint( const PointH3<R>& p, const PointH3<R>& q )
{
  typedef typename R::RT RT;
  RT phw = p.hw();
  RT qhw = q.hw();
  return PointH3<R>( p.hx()*qhw + q.hx()*phw,
                              p.hy()*qhw + q.hy()*phw,
                              p.hz()*qhw + q.hz()*phw,
                              RT(2) * phw * qhw );
}


template <class R>
PointH3<R>
gp_linear_intersection(const PlaneH3<R> &f,
                       const PlaneH3<R> &g,
                       const PlaneH3<R> &h)
{
  return PointH3<R>(
                      det3x3_by_formula(-f.d(), f.b(), f.c(),
                                        -g.d(), g.b(), g.c(),
                                        -h.d(), h.b(), h.c()),
                      det3x3_by_formula( f.a(),-f.d(), f.c(),
                                         g.a(),-g.d(), g.c(),
                                         h.a(),-h.d(), h.c()),
                      det3x3_by_formula( f.a(), f.b(),-f.d(),
                                         g.a(), g.b(),-g.d(),
                                         h.a(), h.b(),-h.d()),
                      det3x3_by_formula( f.a(), f.b(), f.c(),
                                         g.a(), g.b(), g.c(),
                                         h.a(), h.b(), h.c()));
}


template <class R>
CGAL_KERNEL_INLINE
typename R::FT
squared_distance( PointH3<R> const& p, PointH3<R> const& q)
{ return (p-q)*(p-q); }


template <class R>
inline
PlaneH3<R>
bisector( PointH3<R> const& p,
          PointH3<R> const& q)
{ return PlaneH3<R>( midpoint(p,q), q-p); }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
PointH3<R>
centroid( PointH3<R> const& p,
          PointH3<R> const& q,
          PointH3<R> const& r,
          PointH3<R> const& s)
{
   typedef typename R::RT  RT;
   const RT phw(p.hw());
   const RT qhw(q.hw());
   const RT rhw(r.hw());
   const RT shw(s.hw());
   RT hx(p.hx()*qhw*rhw*shw + q.hx()*phw*rhw*shw + r.hx()*phw*qhw*shw
         + s.hx()*phw*qhw*rhw);
   RT hy(p.hy()*qhw*rhw*shw + q.hy()*phw*rhw*shw + r.hy()*phw*qhw*shw
         + s.hy()*phw*qhw*rhw);
   RT hz(p.hz()*qhw*rhw*shw + q.hz()*phw*rhw*shw + r.hz()*phw*qhw*shw
         + s.hz()*phw*qhw*rhw);
   RT hw( phw*qhw*rhw*shw * 4);
   return PointH3<R>(hx, hy, hz, hw);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
PointH3<R>
centroid( PointH3<R> const& p,
          PointH3<R> const& q,
          PointH3<R> const& r)
{
   typedef typename R::RT  RT;
   const RT phw(p.hw());
   const RT qhw(q.hw());
   const RT rhw(r.hw());
   RT hx(p.hx()*qhw*rhw + q.hx()*phw*rhw + r.hx()*phw*qhw);
   RT hy(p.hy()*qhw*rhw + q.hy()*phw*rhw + r.hy()*phw*qhw);
   RT hz(p.hz()*qhw*rhw + q.hz()*phw*rhw + r.hz()*phw*qhw);
   RT hw( phw*qhw*rhw * 3);
   return PointH3<R>(hx, hy, hz, hw);
}

template <class R>
PointH3<R>
circumcenter( PointH3<R> const& p,
              PointH3<R> const& q,
              PointH3<R> const& r,
              PointH3<R> const& s)
{
  typedef typename R::RT RT;

  RT phw( p.hw() );
  RT qhw( q.hw() );
  RT rhw( r.hw() );
  RT shw( s.hw() );

  RT phx( p.hx() );
  RT phy( p.hy() );
  RT phz( p.hz() );
  RT qhx( q.hx() );
  RT qhy( q.hy() );
  RT qhz( q.hz() );
  RT rhx( r.hx() );
  RT rhy( r.hy() );
  RT rhz( r.hz() );
  RT shx( s.hx() );
  RT shy( s.hy() );
  RT shz( s.hz() );

  RT pssq( phx*phx + phy*phy + phz*phz );
  RT qssq( qhx*qhx + qhy*qhy + qhz*qhz );
  RT rssq( rhx*rhx + rhy*rhy + rhz*rhz );
  RT sssq( shx*shx + shy*shy + shz*shz );

  phx *= phw;
  phy *= phw;
  phz *= phw;
  phw *= phw;
  qhx *= qhw;
  qhy *= qhw;
  qhz *= qhw;
  qhw *= qhw;
  rhx *= rhw;
  rhy *= rhw;
  rhz *= rhw;
  rhw *= rhw;
  shx *= shw;
  shy *= shw;
  shz *= shw;
  shw *= shw;

  RT chx =  det4x4_by_formula(phy, phz, pssq, phw,
                                   qhy, qhz, qssq, qhw,
                                   rhy, rhz, rssq, rhw,
                                   shy, shz, sssq, shw );
  RT chy =  det4x4_by_formula(phx, phz, pssq, phw,
                                   qhx, qhz, qssq, qhw,
                                   rhx, rhz, rssq, rhw,
                                   shx, shz, sssq, shw );
  RT chz =  det4x4_by_formula(phx, phy, pssq, phw,
                                   qhx, qhy, qssq, qhw,
                                   rhx, rhy, rssq, rhw,
                                   shx, shy, sssq, shw );
  RT chw =  det4x4_by_formula(phx, phy, phz, phw,
                                   qhx, qhy, qhz, qhw,
                                   rhx, rhy, rhz, rhw,
                                   shx, shy, shz, shw );

 return PointH3<R>( chx, -chy, chz, RT(2)*chw);
}

template <class R>
CGAL_KERNEL_INLINE
PointH3<R>
circumcenter( PointH3<R> const& p,
              PointH3<R> const& q,
              PointH3<R> const& r)
{
  return gp_linear_intersection( PlaneH3<R>(p,q,r),
                                 bisector(p,q),
                                 bisector(p,r));
}

template <class R>
CGAL_KERNEL_INLINE
typename R::FT
squared_radius( const PointH3<R>& p,
                const PointH3<R>& q,
                const PointH3<R>& r )
{
  return squared_distance(p, circumcenter(p, q, r));
}

template <class R>
CGAL_KERNEL_INLINE
typename R::FT
squared_radius( const PointH3<R>& p,
                const PointH3<R>& q,
                const PointH3<R>& r,
                const PointH3<R>& s )
{
  return squared_distance(p, circumcenter(p, q, r, s));
}

CGAL_END_NAMESPACE

#endif // CGAL_BASIC_CONSTRUCTIONSH3_H
