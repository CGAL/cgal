// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_PREDICATES_H
#define CGAL_SPHERE_PREDICATES_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 23
#include <CGAL/Nef_2/debug.h>
#include <vector>

#undef CGAL_forall_iterators
#define CGAL_forall_iterators(x,L)\
for(x = (L).begin(); x != (L).end(); ++x)


namespace CGAL {

/* |spherical_orientation| takes three points of one hemisphere and
returns the orientation of $p_3$ with respect to the halfcircle
through $p_1$ and $p_2$. We map this to the 3d orientation predicate
of the origin, $p_1$, $p_2$, and $p_3$ */

template <class R>
int spherical_orientation(const Sphere_point<R>& p1,
                          const Sphere_point<R>& p2,
                          const Sphere_point<R>& p3)
{ return CGAL::orientation(typename R::Point_3(0,0,0),
                           (typename R::Point_3)p1,
                           (typename R::Point_3)p2,
                           (typename R::Point_3)p3); }


/* |spherical_compare| codes our order of points during the sweep. The
south pole is the first point, the north pole is the last point, then
we order lexicographically. We cover the special case where both
points are part of the equator first. Otherwise we sort according to
the angle of the halfcircle through $S$, $N$, and the points with
respect to the xy-plane. If both lie on the same halfcircle then the
angle with respect to $OS$ decides. The parameter $pos=1$ does
everthing in the positive halfsphere. If $pos=-1$ then we rotate the
whole scenery around the y-axis by $\pi$. Then the x-axis points left
and the z-axis into the equatorial plane. */

template <class R>
bool is_south(const Sphere_point<R>& p, int axis) {
  if(axis==1)
    return (p.hz() >  0 &&
            p.hx() == 0 &&
            p.hy() == 0);

  return (p.hy() <  0 &&
          p.hx() == 0 &&
          p.hz() == 0);
}

template <class R>
bool is_north(const Sphere_point<R>& p, int axis) {
  if(axis==1)
    return (p.hz() <  0 &&
            p.hx() == 0 &&
            p.hy() == 0);

  return (p.hy() >  0 &&
          p.hx() == 0 &&
          p.hz() == 0);
}

template <class R>
int spherical_compare(const Sphere_point<R>& p1,
                      const Sphere_point<R>& p2,
                      int axis, int pos) {

  Sphere_point<R> pS, pN;
  CGAL_assertion(axis>=0 && axis<=2);
  switch(axis) {
  case 0:
    pS=Sphere_point<R>(0,-1,0);
    //    pN=Sphere_point<R>(0,1,0);
    break;
  case 1:
    pS=Sphere_point<R>(0,0,1);
    //    pN=Sphere_point<R>(0,0,-1);
    break;
  case 2:
    pS=Sphere_point<R>(0,-1,0);
    //    pN=Sphere_point<R>(0,1,0);
    break;
  }
  typename R::Direction_3
    d1(p1-CGAL::ORIGIN),
    d2(p2-CGAL::ORIGIN);
  if (d1 == d2) return 0;
  if(is_south(p1,axis) || is_north(p2,axis)) return -1;
  if(is_south(p2,axis) || is_north(p1,axis)) return 1;
  //  if (d1 == dS || d2 == dN) return -1;
  //  if (d1 == dN || d2 == dS) return  1;
  // now no more special cases
  if (axis==0 && (p1.hx()==static_cast<typename R::RT>(0) &&
                  p2.hx()==static_cast<typename R::RT>(0))) {
      int s1 = CGAL_NTS sign(p1.hz());
      int s2 = CGAL_NTS sign(p2.hz());
      if ( s1 != s2 ) return pos * -s1;
      // now s1 == s2
      return -s1 * CGAL::spherical_orientation(p1,Sphere_point<R>(1,0,0),p2);
  } else
  if (axis==1 && (p1.hy()==static_cast<typename R::RT>(0) &&
                  p2.hy()==static_cast<typename R::RT>(0))) {
    int s1 = CGAL_NTS sign(p1.hx());
    int s2 = CGAL_NTS sign(p2.hx());
    if ( s1 != s2 ) return pos * s1;
    // now s1 == s2
    return s1 * CGAL::spherical_orientation(p1,Sphere_point<R>(0,1,0),p2);
  } else
  if (axis==2 && (p1.hz()==static_cast<typename R::RT>(0) &&
                  p2.hz()==static_cast<typename R::RT>(0))) {
    int s1 = CGAL_NTS sign(p1.hx());
    int s2 = CGAL_NTS sign(p2.hx());
    if ( s1 != s2 ) return pos * s1;
    // now s1 == s2
    return s1 * CGAL::spherical_orientation(p1,Sphere_point<R>(0,0,1),p2);
  }
  int sor = CGAL::spherical_orientation(pS,p1,p2);
  if ( sor ) return sor;
  if(axis==0)
    return CGAL::spherical_orientation(Sphere_point<R>(0,0,pos),p2,p1);
  return CGAL::spherical_orientation(Sphere_point<R>(-pos,0,0),p2,p1);
}

/* the next two functions partition the segments in a list into a list
that carries the segment parts that are only part of a halfspace.
Halfcircles are again divided into two equally sized pieces. */

template <class R, class I>
void partition(const Sphere_circle<R>& c, I start, I beyond,
               std::list< Sphere_segment<R> >& Lpos)
{ CGAL_NEF_TRACEN("partition ");
  Sphere_segment<R> s1,s2,s3;
  while ( start != beyond ) { CGAL_NEF_TRACEN("  "<<*start);
    int i = start->intersection(c,s1,s2);
    if (i>1) Lpos.push_back(s2);
    if (i>0) Lpos.push_back(s1);
    ++start;
  }
}

template <class R, class I>
void partition_xy(I start, I beyond,
                  std::list< Sphere_segment<R> >& L, int pos)
{
  Sphere_circle<R> xy_circle(0,0,1), yz_circle(1,0,0);
  Sphere_point<R> S(0,-1,0),N(0,1,0);
  Sphere_segment<R> s1,s2;
  if (pos > 0)  partition(xy_circle,start,beyond,L);
  else partition(xy_circle.opposite(),start,beyond,L);
  typename std::list< Sphere_segment<R> >::iterator it,itl;
  CGAL_NEF_TRACEN("partition_xy ");
  CGAL_forall_iterators(it,L) {
    CGAL_NEF_TRACEN("  "<<*it);
    if ( equal_as_sets(it->sphere_circle(),xy_circle) ) {
      CGAL_NEF_TRACEN("  splitting xy seg "<<*it);
      int n1 =  it->intersection(yz_circle,s1,s2);
      if (n1 > 1 && !s2.is_degenerate()) L.insert(it,s2);
      if (n1 > 0 && !s1.is_degenerate()) L.insert(it,s1);
      int n2 =  it->intersection(yz_circle.opposite(),s1,s2);
      if (n2 > 1 && !s2.is_degenerate()) L.insert(it,s2);
      if (n2 > 0 && !s1.is_degenerate()) L.insert(it,s1);
      itl = it; --it; L.erase(itl);
      // at least one item was appended
    }
  }
  CGAL_forall_iterators(it,L) {
    if ( it->is_halfcircle() ) {
      CGAL_NEF_TRACEN("  splitting halfcircle "<<*it);
      Sphere_segment<R> s1,s2;
      it->split_halfcircle(s1,s2);
      *it = s2; L.insert(it,s1);
    }
  }
  // append 4 xy-equator segments:
  Sphere_segment<R> sp(S,N,xy_circle);
  Sphere_segment<R> sm(S,N,xy_circle.opposite());
  Sphere_segment<R> s[4];
  sp.split_halfcircle(s[0],s[1]);
  sm.split_halfcircle(s[2],s[3]);
  L.insert(L.end(),s,s+4);
}


/* |intersection| calculates the intersection of a halfspace $h$
(defined via a great circle $c$) with a sphere segment
$s$. Interesting are the boundary cases. If an endpoint is part of
$h^0$, but the part of $s$ incident to the endpoint is not part of
$h^{0+}$ then we return a trivial segment.

*/

/*
template <typename R>
int Sphere_segment<R>::
intersection(const CGAL::Sphere_circle<R>& c, std::vector<Sphere_segment<R> >& s) const {

  CGAL_NEF_TRACEN("    intersection "<<*this<<" "<<c);
  if ( is_degenerate() ) {
    CGAL_NEF_TRACEN("    degenerate");
    s.push_back(*this);
    return 1;
  }

  CGAL::Oriented_side or1 = c.oriented_side(source());
  CGAL::Oriented_side or2 = c.oriented_side(target());

  CGAL_NEF_TRACEN("    Orientation " << or1 << " " << or2);

  if ( or1 == CGAL::opposite(or2) && or1 != or2 ) {
    // it is sure that $s$ intersects $h$ in its interior. the
    //   question is which part is in the halfspace $h^+$.
    CGAL_NEF_TRACEN("    opposite");
    Sphere_point<R> i1 = CGAL::intersection(ptr()->c_,c);
    if ( !has_on(i1) )
      i1 = i1.antipode();
    s.push_back(Sphere_segment<R>(source(),i1,sphere_circle()));
    s.push_back(Sphere_segment<R>(i1,target(),sphere_circle()));
    return 2;
  }
  else if ( or1 == CGAL::ON_ORIENTED_BOUNDARY &&
            or2 == CGAL::ON_ORIENTED_BOUNDARY ) {
    // if both ends of $s$ are part of $h$ then $s$ is a halfcircle or
    //   $s$ is fully part of $h$.  In this case we have to check if the
    //   halfcircle is not part of $h^-$.  This can be formulated by an
    //   orientation test of the point $p$ at the tip of the normal of
    //   |s.sphere_circle()| with respect to the plane through the
    //   endpoints of $s$ and the tip of the normal of $h$.
    CGAL_NEF_TRACEN("    both in plane");
    s.push_back(*this);
    return 1;
  }
  else if ( or1 != CGAL::ON_NEGATIVE_SIDE &&
            or2 != CGAL::ON_NEGATIVE_SIDE ) {
    // this case covers the endpoints of $s$ as part of the closed
    //   oriented halfspace $h^{0+}$. At least one is part of
    //   $h^{+}$. If $s$ is not long then the segment must be fully part
    //   of $h^{0+}$. Otherwise if $s$ is long, then we at both ends
    //   there are subsegments as part of $h^{0+}$ (one might be
    //   degenerate).
    if ( is_long() ) {
      Sphere_point<R> i1 = CGAL::intersection(ptr()->c_,c);
      Sphere_point<R> i2 = i1.antipode();
      Sphere_segment<R> so(i1,i2,sphere_circle());
      if ( so.has_on(source()) && so.has_on(target()) )
        std::swap(i1,i2);
      // now source(),i1,i2,target() are circularly ordered
      s.push_back(Sphere_segment(source(),i1,sphere_circle()));
      s.push_back(Sphere_segment(i1,i2,sphere_circle()));
      s.push_back(Sphere_segment(i2,target(),sphere_circle()));
      //      CGAL_NEF_TRACEN("    both >= plane, long "<<s1<<s2);
      return 3;
    } // now short:
    CGAL_NEF_TRACEN("    both >= plane, short ");
    s.push_back(*this);
    return 1;
  }
  else if ( or1 != CGAL::ON_POSITIVE_SIDE &&
            or2 != CGAL::ON_POSITIVE_SIDE ) {
    // either both endpoints of $s$ are in $h^-$ or one is in $h^-$
    //   and one on $h^0$.
    if ( is_long() ) {
      Sphere_point<R> i1 = CGAL::intersection(ptr()->c_,c);
      Sphere_point<R> i2 = i1.antipode();
      Sphere_segment so(i1,i2,sphere_circle());
      //      CGAL_NEF_TRACEN("    both <= plane, long"<<so);
      if ( so.has_on(source()) && so.has_on(target()) )
      { so = so.complement(); }
      //      CGAL_NEF_TRACEN("    both <= plane, long"<<so);
      s.push_back(Sphere_segment(source(),so.source(),sphere_circle()));
      s.push_back(so);
      s.push_back(Sphere_segment(so.target(),target(),sphere_circle()));
      return 3;
    } // now short
    CGAL_NEF_TRACEN("    both <= plane, short");
    s.push_back(*this);
    return 1;
  }

  CGAL_error_msg("Oops, forgot some case.");
  return -1;
}
*/

template <typename R>
int Sphere_segment<R>::
intersection(const CGAL::Sphere_circle<R>& c,
             Sphere_segment<R>& s1, Sphere_segment<R>& s2) const
{
  CGAL_NEF_TRACEN("    intersection "<<*this<<" "<<c);
  if ( is_degenerate() ) { CGAL_NEF_TRACEN("    degenerate");
    if ( !c.has_on_negative_side(source()) )
    { s1 = *this; return 1; }
    return 0;
  }
  CGAL::Oriented_side or1 = c.oriented_side(source());
  CGAL::Oriented_side or2 = c.oriented_side(target());
  if ( or1 == CGAL::opposite(or2) && or1 != or2 ) {
// it is sure that $s$ intersects $h$ in its interior. the
//       question is which part is in the halfspace $h^+$.
    CGAL_NEF_TRACEN("    opposite");
    Sphere_point<R> i1 = CGAL::intersection(this->ptr()->c_,c);
    if ( !has_on(i1) ) i1 = i1.antipode();
    if ( or1 == CGAL::ON_POSITIVE_SIDE )
      s1 = Sphere_segment<R>(source(),i1,sphere_circle());
    else if ( or2 == CGAL::ON_POSITIVE_SIDE )
      s1 = Sphere_segment<R>(i1,target(),sphere_circle());
    else CGAL_error_msg("no intersection.");
    return 1;
  }
  else if ( or1 == CGAL::ON_ORIENTED_BOUNDARY &&
            or2 == CGAL::ON_ORIENTED_BOUNDARY ) {
// if both ends of $s$ are part of $h$ then $s$ is a halfcircle or
//    $s$ is fully part of $h$.  In this case we have to check if the
//    halfcircle is not part of $h^-$.  This can be formulated by an
//    orientation test of the point $p$ at the tip of the normal of
//    |s.sphere_circle()| with respect to the plane through the
//   endpoints of $s$ and the tip of the normal of $h$.
    CGAL_NEF_TRACEN("    both in plane");
    if ( source() != target().antipode() ) {
      s1 = *this; return 1;
    }
    // now this is a halfcircle
    bool halfcircle_notin_hminus =
      (CGAL::orientation(source(),target(),
                         CGAL::ORIGIN + c.orthogonal_vector(),
                         CGAL::ORIGIN + sphere_circle().orthogonal_vector())
       != CGAL::POSITIVE);
    CGAL_NEF_TRACE("    ");CGAL_NEF_TRACEV(halfcircle_notin_hminus);
    if ( halfcircle_notin_hminus ) { s1 = *this; return 1; }
    else {
      s1 = Sphere_segment(source(),source(),sphere_circle());
      s2 = Sphere_segment(target(),target(),sphere_circle());
      return 2;
    }
  }
  else if ( or1 != CGAL::ON_NEGATIVE_SIDE &&
            or2 != CGAL::ON_NEGATIVE_SIDE ) {
// this case covers the endpoints of $s$ as part of the closed
//   oriented halfspace $h^{0+}$. At least one is part of
//   $h^{+}$. If $s$ is not long then the segment must be fully part
//   of $h^{0+}$. Otherwise if $s$ is long, then we at both ends
//   there are subsegments as part of $h^{0+}$ (one might be
//   degenerate).
    if ( is_long() ) {
      Sphere_point<R> i1 = CGAL::intersection(this->ptr()->c_,c);
      Sphere_point<R> i2 = i1.antipode();
      Sphere_segment<R> so(i1,i2,sphere_circle());
      if ( so.has_on(source()) && so.has_on(target()) )
        std::swap(i1,i2);
      // now source(),i1,i2,target() are circularly ordered
      s1 = Sphere_segment(source(),i1,sphere_circle());
      s2 = Sphere_segment(i2,target(),sphere_circle());
      CGAL_NEF_TRACEN("    both >= plane, long "<<s1<<s2);
      return 2;
    } // now short:
    CGAL_NEF_TRACEN("    both >= plane, short ");
    s1=*this; return 1;
  }
  else if ( or1 != CGAL::ON_POSITIVE_SIDE &&
            or2 != CGAL::ON_POSITIVE_SIDE ) {
// either both endpoints of $s$ are in $h^-$ or one is in $h^-$
//     and one on $h^0$.
    if ( is_long() ) {
      Sphere_point<R> i1 = CGAL::intersection(this->ptr()->c_,c);
      Sphere_point<R> i2 = i1.antipode();
      Sphere_segment<R> so(i1,i2,sphere_circle());
      CGAL_NEF_TRACEN("    both <= plane, long"<<so);
      if ( so.has_on(source()) && so.has_on(target()) )
      { so = so.complement(); }
      CGAL_NEF_TRACEN("    both <= plane, long"<<so);
      s1 = so; return 1;
    } // now short
    CGAL_NEF_TRACEN("    both <= plane, short");
    if ( or1 == CGAL::ON_ORIENTED_BOUNDARY )
    { s1 = Sphere_segment<R>(source(),source(),sphere_circle()); return 1; }
    if ( or2 == CGAL::ON_ORIENTED_BOUNDARY )
    { s1 = Sphere_segment<R>(target(),target(),sphere_circle()); return 1; }
    return 0;
  }

  CGAL_error_msg("Oops, forgot some case.");
  return -1;
}

} //namespace CGAL

#endif //CGAL_SPHERE_PREDICATES_H
