// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_SEGMENT_H
#define CGAL_SPHERE_SEGMENT_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <vector>

namespace CGAL {

template <class R_> class Sphere_segment_rep 
{ 
  typedef typename R_::Point_3 Point_3;
  typedef typename R_::Plane_3 Plane_3;
  typedef Sphere_point<R_>  Point;
  typedef Sphere_circle<R_> Circle;
  typedef Sphere_segment_rep<R_> Rep;
  friend class Sphere_segment<R_>;
public:
Sphere_segment_rep() { ps_ = pt_ = Point(); c_ = Circle(); }

Sphere_segment_rep(const Point& p1, const Point& p2,
		   bool shorter_arc=true) :
  ps_(p1), pt_(p2), c_(Plane_3(p1,p2,Point_3(0,0,0))) 
{ // warning stays as reminder that one gets an arbitrary plane equation
  // in this degenerate case
  CGAL_warning(p1 != p2.antipode());
  CGAL_assertion(p1 != p2.antipode());
  if ( p1 == p2 ) { 
    Plane_3 h(Point_3(0,0,0),(p1-CGAL::ORIGIN));
    c_ = Sphere_circle<R_>(Plane_3(Point_3(0,0,0),h.base1()));
  }
  if (!shorter_arc) c_ = c_.opposite(); 
  CGAL_exactness_assertion(c_.has_on(p1) && c_.has_on(p2));
}

Sphere_segment_rep(const Point& p1, const Point& p2, const Circle& c) :
  ps_(p1), pt_(p2), c_(c) 
{ CGAL_assertion(c.has_on(p1)&&c.has_on(p2)); }

Sphere_segment_rep(const Circle& c1, 
                   const Circle& c2) : c_(c1) 
{ CGAL_assertion(!equal_as_sets(c1,c2)); 
  ps_ = intersection(c1,c2);
  pt_ = ps_.antipode();
  if ( R_::orientation(Point_3(0,0,0),ps_,pt_,
                   CGAL::ORIGIN + c_.orthogonal_vector()) !=
       CGAL::POSITIVE ) std::swap(ps_,pt_);
}

Sphere_segment_rep(const Rep& r) : 
  ps_(r.ps_), pt_(r.pt_), c_(r.c_) {}

Rep& operator=(const Rep& r) 
{ ps_=r.ps_; pt_=r.pt_; c_=r.c_; return *this; }

protected:
  Sphere_point<R_>  ps_,pt_;
  Sphere_circle<R_> c_;
};


/*{\Moptions print_title=yes }*/ 
/*{\Manpage{Sphere_segment}{R}{Segments on the unit sphere}{s}}*/

template <class R_> 
class Sphere_segment : 
  public Handle_for< Sphere_segment_rep<R_> > {

/*{\Mdefinition An object |\Mvar| of type |\Mname| is a segment in the
surface of a unit sphere that is part of a great circle trough the
origin. Sphere segments are represented by two sphere points $p$ and
$q$ plus an oriented plane $h$ that contains $p$ and $q$. The plane
determines the sphere segment. Let $c$ be the circle in the
intersection of $h$ and $S_2$. Then $s$ is that part of $c$ that is
swept, when we rotate $p$ into $q$ in counterclockwise rotation around
the normal vector of $h$ as seen from the positive halfspace.}*/  

public:

/*{\Mtypes 6}*/

typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R_::RT RT;
/*{\Mtypemember ring number type.}*/

typedef Sphere_segment_rep<R_> Rep;
typedef Handle_for<Rep>        Base;

typedef typename R_::Point_3  Point_3;
typedef typename R_::Vector_3 Vector_3;
typedef typename R_::Plane_3  Plane_3;
typedef typename R_::Line_3   Line_3;
typedef Sphere_segment<R_> Self;

/*{\Mcreation 4}*/

Sphere_segment() : Base() {}

Sphere_segment(const Sphere_point<R>& p1, const Sphere_point<R>& p2,
               bool shorter_arc=true) :  Base(Rep(p1,p2,shorter_arc)) {}
/*{\Mcreate creates a spherical segment spanning the shorter arc
from |p1| to |p2| if |shorter_arc == true|. Otherwise the longer
arc is created. \precond |p1 != p2| and |p1 != p2.antipode()|.}*/


Sphere_segment(const Sphere_point<R>& p1, const Sphere_point<R>& p2,
               const Sphere_circle<R>& c) : Base(Rep(p1,p2,c)) {}
/*{\Mcreate creates a spherical segment spanning the arc
from |p1| to |p2| as part of the oriented circle |c| 
(|p1 == p2| or |p1 == p2.antipode()| are possible.)
\precond |p1| and |p2| are contained in |c|.}*/

Sphere_segment(const Sphere_circle<R>& c1, 
               const Sphere_circle<R>& c2) : Base(Rep(c1,c2)) {}
/*{\Mcreate creates the spherical segment as part of |c1|
that is part of the halfsphere left of the oriented circle |c2|.
\precond |c1 != c2| as unoriented circles.}*/

Sphere_segment(const Self& s) : Base(s) {}

/*{\Moperations 4 2}*/

const Sphere_point<R>& source() const { return this->ptr()->ps_; }
/*{\Mop the source point of |\Mvar|.}*/
const Sphere_point<R>& target() const { return this->ptr()->pt_; }
/*{\Mop the target point of |\Mvar|.}*/
const Sphere_circle<R>& sphere_circle() const { return this->ptr()->c_; }
/*{\Mop the great circle supporting |\Mvar|.}*/

Sphere_segment<R> opposite() const 
/*{\Mop returns the sperical segment oriented from |target()|
  to |source()| with the same point set as |\Mvar|. }*/
{ return Sphere_segment<R>(
    target(),source(),sphere_circle().opposite()); }
  
Sphere_segment<R> complement() const 
/*{\Mop returns the sperical segment oriented from |target()|
  to |source()| with the point set completing |\Mvar| to a
  full circle. }*/
{ return Sphere_segment<R>(target(),source(),sphere_circle()); }


int intersection(const Sphere_circle<R>& c,
                 std::vector<Sphere_segment<R> >& s) const;

int intersection(const Sphere_circle<R>& c,
                 Sphere_segment<R>& s1, Sphere_segment<R>& s2) const;
/*{\Mop returns the number of non-trivial connected components
of the intersection of |\Mvar| and the closed halfsphere left of
|c|.}*/

Sphere_point<R> intersection(const Sphere_segment<R>& so) const;
/*{\Mop returns the point of intersection of |\Mvar| and
|so|. \precond |\Mvar| and |so| do intersect.}*/

void split_halfcircle(Sphere_segment<R>& s1, 
                      Sphere_segment<R>& s2) const
/*{\Mop splits a halfcircle into two equally sized segments.
\precond |\Mvar| is a halfcircle.}*/
{ CGAL_assertion( is_halfcircle() );
  Plane_3 h(Point_3(0,0,0),(target()-CGAL::ORIGIN));
  Sphere_point<R> p = 
    CGAL::intersection(sphere_circle(),Sphere_circle<R>(h));
  if ( !has_on(p) ) p = p.antipode();
  s1 = Sphere_segment<R>(this->ptr()->ps_,p,this->ptr()->c_);
  s2 = Sphere_segment<R>(p,this->ptr()->pt_,this->ptr()->c_);
}

bool is_short() const 
/*{\Mop a segment is short iff it is shorter than a halfcircle.}*/
{ 
  return R().orientation_3_object()(Point_3(0,0,0), source(), target(),
				    CGAL::ORIGIN + this->ptr()->c_.orthogonal_vector()) 
    == CGAL::POSITIVE; }

bool is_long() const 
/*{\Mop a segment is long iff it is longer than a halfcircle.}*/
{ return R().orientation_3_object()(Point_3(0,0,0), source(), target(),
				    CGAL::ORIGIN + this->ptr()->c_.orthogonal_vector()) 
    == CGAL::NEGATIVE; }

bool is_degenerate() const { return source() == target(); }
/*{\Mop return true iff |\Mvar| is degenerate.}*/

bool is_halfcircle() const { return source().antipode() == target(); }
/*{\Mop return true iff |\Mvar| is a halfcircle.}*/

bool has_on(const Sphere_point<R>& p) const;
/*{\Mop return true iff |\Mvar| contains |p|.}*/
bool has_on_after_intersection(const Sphere_point<R>& p) const;

bool has_in_relative_interior(const Sphere_point<R>& p) const;
/*{\Mop return true iff |\Mvar| contains |p| in
its relative interior.}*/

bool operator==(const Sphere_segment<R>& so) const
{ return source() == so.source() && target() == so.target() && 
         (source() == target() || 
	  sphere_circle() == so.sphere_circle()); }

bool operator!=(const Sphere_segment<R>& so) const
{ return !operator==(so); }

};

template <typename R>
bool do_intersect_internally(const Sphere_segment<R>& s1,
                             const Sphere_segment<R>& s2,
                             Sphere_point<R>& p);
/*{\Mfunc return true iff |s1| and |s2| intersect internally
(non-degenerately).  If |true| the parameter |p| returns the point of
intersection.}*/


template <typename R>
std::ostream& operator<<(std::ostream& os, 
                         const CGAL::Sphere_segment<R>& s)
{ os << s.source()<<" "<<s.target()<<" "<<
  s.sphere_circle().plane()<<" "; return os; }

template <typename R>
std::istream& operator>>(std::istream& is, 
                         CGAL::Sphere_segment<R>& s)
{ CGAL::Sphere_point<R>  p1,p2;
  CGAL::Sphere_circle<R> c;
  if ( !(is >> p1 >> p2 >> c) ) return is; 
  s = CGAL::Sphere_segment<R>(p1,p2,c); 
  return is; }


template <typename R>
std::pair< Sphere_segment<R>,Sphere_segment<R> >
Sphere_circle<R>::split_at(const Sphere_point<R>& p) const
{ CGAL_assertion(has_on(p));
  Sphere_point<R> q(p.antipode());
  return Sphere_segment_pair(
    Sphere_segment<R>(p,q,*this),
    Sphere_segment<R>(p,q,this->opposite()));
}

template <typename R>
std::pair< Sphere_segment<R>,Sphere_segment<R> >
Sphere_circle<R>::split_at_xy_plane() const
{ Self xycircle(0,0,1), yzcircle(1,0,0);
  if ( !equal_as_sets(xycircle,*this) )
    return split_at(intersection(*this,xycircle));
  else 
    return split_at(intersection(*this,yzcircle));
}


/* Contains maps to two orientation checks with the wedge
spanned by the source and the target with planes orthogonal
to the supporting plane of $p$ and $q$. The logic depends on
the segments length: long or short. */

template <typename R> 
bool Sphere_segment<R>::
has_on(const CGAL::Sphere_point<R>& p) const 
{ if ( !sphere_circle().has_on(p) ) return false;
  if ( !is_long() ) { 
    return orientation(Point_3(0,0,0),
                       CGAL::ORIGIN + sphere_circle().orthogonal_vector(),
                       source(),p) !=
           CGAL::NEGATIVE &&
           orientation(Point_3(0,0,0),target(),
                       CGAL::ORIGIN + 
                       sphere_circle().orthogonal_vector(),p) !=
           CGAL::NEGATIVE;
  } else { 
    return orientation(Point_3(0,0,0),
                       CGAL::ORIGIN + sphere_circle().orthogonal_vector(),
                       source(),p) !=
           CGAL::NEGATIVE ||
           orientation(Point_3(0,0,0),target(),
                       CGAL::ORIGIN + 
                       sphere_circle().orthogonal_vector(),p) !=
           CGAL::NEGATIVE;
  }
}

template <typename R> 
bool Sphere_segment<R>::
has_on_after_intersection(const CGAL::Sphere_point<R>& p) const {
    return orientation(Point_3(0,0,0),
                       CGAL::ORIGIN + sphere_circle().orthogonal_vector(),
                       source(),p) !=
           CGAL::NEGATIVE &&
           orientation(Point_3(0,0,0),target(),
                       CGAL::ORIGIN + 
                       sphere_circle().orthogonal_vector(),p) !=
           CGAL::NEGATIVE;
}

template <typename R> 
bool Sphere_segment<R>::
has_in_relative_interior(const CGAL::Sphere_point<R>& p) const 
{ if ( !sphere_circle().has_on(p) ) return false;
  if ( !is_long() ) { 
    return orientation(Point_3(0,0,0),
                       CGAL::ORIGIN + sphere_circle().orthogonal_vector(),
                       source(),p) == CGAL::POSITIVE &&
           orientation(Point_3(0,0,0),target(),
                       CGAL::ORIGIN + 
                       sphere_circle().orthogonal_vector(),p) ==
           CGAL::POSITIVE;
  } else { 
    return orientation(Point_3(0,0,0),
                       CGAL::ORIGIN + sphere_circle().orthogonal_vector(),
                       source(),p) == CGAL::POSITIVE ||
           orientation(Point_3(0,0,0),target(),
                       CGAL::ORIGIN + 
                       sphere_circle().orthogonal_vector(),p) ==
           CGAL::POSITIVE;
  }
}



/* Intersection of two sphere segments. It does not work if the two
   involved planes are equal as sets. */

template <typename R> 
Sphere_point<R> Sphere_segment<R>::
intersection(const Sphere_segment<R>& s) const
{
  CGAL_assertion(!equal_as_sets(sphere_circle(),s.sphere_circle()));
  Sphere_point<R> res = 
    CGAL::intersection(sphere_circle(),s.sphere_circle());
  if ( has_on(res) && s.has_on(res) ) return res;
  return res.antipode();
}

template <typename R> 
bool do_intersect_internally(const Sphere_segment<R>& s1,
                             const Sphere_segment<R>& s2,
                             Sphere_point<R>& p)
{
  if ( equal_as_sets(s1.sphere_circle(),s2.sphere_circle()) )
    return false;
  p = CGAL::intersection(s1.sphere_circle(),s2.sphere_circle());
  if ( s1.has_in_relative_interior(p) && 
       s2.has_in_relative_interior(p) ) return true;
  p = p.antipode();
  if ( s1.has_in_relative_interior(p) && 
       s2.has_in_relative_interior(p) ) return true;
  return false;
}


template <typename R> 
bool do_intersect_internally(const Sphere_circle<R>& c1,
                             const Sphere_segment<R>& s2,
                             Sphere_point<R>& p)
{
  if ( equal_as_sets(c1,s2.sphere_circle()) )
    return false;
  p = CGAL::intersection(c1,s2.sphere_circle());
  if ( s2.has_in_relative_interior(p) ) return true;
  p = p.antipode();
  if ( s2.has_in_relative_interior(p) ) return true;
  return false;
}


} //namespace CGAL
#endif //CGAL_SPHERE_SEGMENT_H
