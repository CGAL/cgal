#ifndef CGAL_INTERSECTION_OBJECTSHD_H
#define CGAL_INTERSECTION_OBJECTSHD_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

/*{\Manpage{Line_line_intersectionHd}{R}{intersecting two lines}}*/

template <class R> 
class Line_line_intersectionHd {

typedef typename R::RT RT;
typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d Point_d;

public:
enum Intersection_result { NO, POINT, LINE };

Intersection_result operator()(
  const Point_d& s1, const Point_d& t1,
  const Point_d& s2, const Point_d& t2,
  Point_d& p, FT& l1, FT& l2) 
/*{\Mfunop returns |NO| if the lines which are represented by |s1t1|
and |s2t2| don't intersect, returns |POINT| if they intersect in a
unique point, and returns LINE if they are identical. In the |POINT|
case the point of intersection is assigned to |p|.  Then |p = s1 + l1
* t1-s1| and |p = s2 + l2 * t2-s2|. \precond none of the point pairs
is degenerate.}*/
{ 
  int d = s1.dimension(); 
  CGAL_assertion_msg(d==s2.dimension(),
    "intersection: dimensions disagree!"); 
  typename LA::Matrix M(d,2),S; 
  typename LA::Vector b(d), lambda(2), c; 
  RT D; 
  
  RT s1w = s1.homogeneous(d); 
  RT t1w = t1.homogeneous(d); 
  RT s2w = s2.homogeneous(d); 
  RT t2w = t2.homogeneous(d); 
  RT g1w = s1w*t1w; 
  RT g2w = s2w*t2w; 
  RT t12w = t1w*t2w; 

  /* init $d \times 2$ - matrix |M| and $d$ - vector |b| */
  for (int i = 0; i < d; i++) { 
    M(i,0) = g2w * (t1.homogeneous(i) * s1w - s1.homogeneous(i) * t1w); 
    M(i,1) = g1w * (s2.homogeneous(i) * t2w - t2.homogeneous(i) * s2w); 
    b[i]   = t12w * (s2.homogeneous(i) * s1w - s1.homogeneous(i) * s2w); 
  }

  if (LA::linear_solver(M,b,lambda,D,S,c)) { 
    if (S.column_dimension()>0) return LINE;
    l1 = R::make_FT(lambda[0],D);
    l2 = R::make_FT(lambda[1],D);
    p = s1 + l1 * (t1 - s1); 
    return POINT; 
  }
  return NO; 
}

};

/*{\Manpage {Line_hyperplane_intersectionHd}{R} 
{intersecting a line and a hyperplane}}*/

template <class R> 
class Line_hyperplane_intersectionHd {
typedef typename R::RT RT;
typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d      Point_d;
typedef typename R::Vector_d     Vector_d;
typedef typename R::Hyperplane_d Hyperplane_d;

public:
enum Intersection_result { NO, POINT, LINE };

Intersection_result operator()(const Point_d& s, const Point_d& t,
  const Hyperplane_d& h, Point_d& p, FT& lambda) 
/*{\Mfunop returns |NO| if the line represented by |s1t1| and the
hyperplane |h| don't intersect, returns |POINT| if they intersect in a
unique point, and returns LINE if the line is part of the
hyperplane. In the |POINT| case the point of intersection is assigned
to |p|.  Then |p = s1 + lambda * t1-s1|.  \precond the point pair is
not degenerate.}*/
{ 
  CGAL_assertion_msg((h.dimension()==s.dimension() && 
                      h.dimension()==t.dimension()), 
  "Line_hyperplane_intersection_d: dimensions do not agree.");

  int d = h.dimension();
  RT S(0),T(0);
  for (int i=0; i<=d; ++i) {
    S += h[i]*s.homogeneous(i);
    T += h[i]*t.homogeneous(i);
  }
  bool s_contained = CGAL_NTS is_zero(S), 
       t_contained = CGAL_NTS is_zero(T);
  if (s_contained && t_contained) { p = s; return LINE; }
  if (s_contained) { p = s; return POINT; }
  if (t_contained) { p = t; return POINT; }
  // now the simple cases are done 

  RT D = S * t.homogeneous(d) - T * s.homogeneous(d);
  if (CGAL_NTS is_zero(D)) return NO;

  typename LA::Vector homog(d + 1);
  for (int i = 0; i < d; ++i)
    homog[i] = S * t.homogeneous(i) - T * s.homogeneous(i); 
  homog[d] = D;
  p = Point_d(d,homog.begin(),homog.end());
  lambda = R::make_FT(S * t.homogeneous(d), D);
  return POINT;
}

};


template <class R>
class Line_d_Line_d_pairH {
public:
  enum Intersection_result {NO, POINT, LINE};
  typedef typename R::Point_d Point_d;
  typedef typename R::Line_d  Line_d;
  typedef typename R::FT FT;
protected:
  Line_d _l1, _l2;
  bool _known;
  Intersection_result _result;
  Point_d _ip;
public:
  Line_d_Line_d_pairH() : _known(false) {}
  Line_d_Line_d_pairH(const Line_d& l1, const Line_d& l2) 
    : _l1(l1), _l2(l2), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Line_d& result);
};

template <class R>
Line_d_Line_d_pairH<R>::Intersection_result
Line_d_Line_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  //CGAL_assertion(!_l1.is_degenerate()&&!_l2.is_degenerate());
  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_l1.point(0),_l1.point(1),
              _l2.point(0),_l2.point(1),
              _ip,l1,l2);

  if (res == Int_obj_type::LINE)  { return _result = LINE; }
  if (res == Int_obj_type::POINT) { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool Line_d_Line_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool Line_d_Line_d_pairH<R>::intersection(Line_d& l)
{ if (!_known) intersection_type();
  if (_result != LINE) return false;
  l = _l1; return true;
}

template <class R>
class Line_d_Ray_d_pairH {
public:
  enum Intersection_result {NO, POINT, RAY};
  typedef typename R::Point_d Point_d;
  typedef typename R::Ray_d   Ray_d;
  typedef typename R::Line_d  Line_d;
  typedef typename R::FT FT;
protected:
  Line_d _l; Ray_d _r;
  bool _known;
  Intersection_result _result;
  Point_d _ip; 
public:
  Line_d_Ray_d_pairH() : _known(false) {}
  Line_d_Ray_d_pairH(const Line_d& l, const Ray_d& r) 
    : _l(l), _r(r), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Ray_d& result);
};

template <class R>
Line_d_Ray_d_pairH<R>::Intersection_result
Line_d_Ray_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  //CGAL_assertion(!_l.is_degenerate()&&!_r.is_degenerate());
  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_l.point(0),_l.point(1),
              _r.point(0),_r.point(1),
              _ip,l1,l2);

  if ( res == Int_obj_type::LINE )  { return _result = RAY; }
  if ( res == Int_obj_type::POINT &&
       l2 >= FT(0) ) 
  { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool Line_d_Ray_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool Line_d_Ray_d_pairH<R>::intersection(Ray_d& r)
{ if (!_known) intersection_type();
  if (_result != RAY) return false;
  r = _r; return true;
}

template <class R>
class Line_d_Segment_d_pairH {
public:
  enum Intersection_result {NO, POINT, SEGMENT};
  typedef typename R::Point_d Point_d;
  typedef typename R::Segment_d Segment_d;
  typedef typename R::Line_d  Line_d;
  typedef typename R::FT FT;
protected:
  Line_d _l; Segment_d _s;
  bool _known;
  Intersection_result _result;
  Point_d _ip; 
public:
  Line_d_Segment_d_pairH() : _known(false) {}
  Line_d_Segment_d_pairH(const Line_d& l, const Segment_d& s) 
    : _l(l), _s(s), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Segment_d& result);
};

template <class R>
Line_d_Segment_d_pairH<R>::Intersection_result
Line_d_Segment_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  //CGAL_assertion(!_l.is_degenerate());
  if ( _s.is_degenerate() ) {
    if ( _l.has_on(_s.point(0)) ) {
      _ip = _s.point(0);
      return _result = POINT;
    }
    return _result = NO; 
  }
  // _s not degenerate
  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_l.point(0),_l.point(1),
              _s.point(0),_s.point(1),
              _ip,l1,l2);

  if ( res == Int_obj_type::LINE )  { return _result = SEGMENT; }
  if ( res == Int_obj_type::POINT &&
       FT(0) <= l2 && l2 <= FT(1) ) 
  { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool Line_d_Segment_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool Line_d_Segment_d_pairH<R>::intersection(Segment_d& s)
{ if (!_known) intersection_type();
  if (_result != SEGMENT) return false;
  s = _s; return true;
}


template <class R>
class Ray_d_Ray_d_pairH {
public:
  enum Intersection_result { NO, POINT, SEGMENT, RAY };
  typedef typename R::FT FT;
  typedef typename R::Point_d Point_d;
  typedef typename R::Segment_d Segment_d;
  typedef typename R::Ray_d Ray_d;

protected:
  Ray_d _r1, _r2;
  bool _known;
  Intersection_result _result;
  Point_d _ip; Segment_d _is; Ray_d _ir;
public:
  Ray_d_Ray_d_pairH() : _known(false) {}
  Ray_d_Ray_d_pairH(const Ray_d& r1, const Ray_d& r2)
    : _r1(r1), _r2(r2), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Segment_d& result);
  bool intersection(Ray_d& result);
};

template <class R>
Ray_d_Ray_d_pairH<R>::Intersection_result
Ray_d_Ray_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  //CGAL_assertion((!_r1.is_degenerate()&&!_r2.is_degenerate());
  // none of the lines should be trivial

  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_r1.point(0),_r1.point(1),
              _r2.point(0),_r2.point(1),_ip,l1,l2);

  if (res == Int_obj_type::LINE)   
  {
    if ( _r1.direction() == _r2.direction() ) {
      if ( _r1.has_on(_r2.source()) ) _ir = _r2;
      else _ir = _r1; 
      return _result = RAY;
    }
    // now oppositely directed:
    if ( _r1.has_on(_r2.source()) ) {
      if ( _r1.source() != _r2.source() ) 
      { _is = Segment_d(_r1.source(),_r2.source()); 
        return _result = SEGMENT; }
      else
      { _ip = _r1.source(); return _result = POINT; }
    }
    return _result = NO;
  }


  if (res == Int_obj_type::POINT && 
      FT(0) <= l1 && FT(0) <= l2) 
  { return _result = POINT; }
  return _result = NO; 
  // now not parallel
}

template <class R>
bool
Ray_d_Ray_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool
Ray_d_Ray_d_pairH<R>::intersection(Segment_d& s)
{ if (!_known) intersection_type();
  if (_result != SEGMENT) return false;
  s = _is; return true;
}

template <class R>
bool
Ray_d_Ray_d_pairH<R>::intersection(Ray_d& r)
{ if (!_known) intersection_type();
  if (_result != RAY) return false;
  r = _ir; return true;
}

template <class R>
class Ray_d_Segment_d_pairH {
public:
  enum Intersection_result { NO, POINT, SEGMENT };
  typedef typename R::FT FT;
  typedef typename R::Point_d Point_d;
  typedef typename R::Segment_d Segment_d;
  typedef typename R::Ray_d Ray_d;

protected:
  Ray_d _r; Segment_d _s;
  bool _known;
  Intersection_result _result;
  Point_d _ip; Segment_d _is;  
public:
  Ray_d_Segment_d_pairH() : _known(false) {}
  Ray_d_Segment_d_pairH(const Ray_d& r, const Segment_d& s)
    : _r(r), _s(s), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Segment_d& result);
};

template <class R>
Ray_d_Segment_d_pairH<R>::Intersection_result
Ray_d_Segment_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  //CGAL_assertion(!_r.is_degenerate());
  if ( _s.is_degenerate() ) {
    if ( _r.has_on(_s.point(0)) ) 
    { _ip = _s.point(0); return _result = POINT; }
    return _result = NO; 
  }
  // now s is not degenerate
  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_r.point(0),_r.point(1),
              _s.point(0),_s.point(1),_ip,l1,l2);

  if (res == Int_obj_type::LINE)   
  {
    Point_d p1 = _s.point(0), p2 = _s.point(1);
    if ( _s.direction() != _r.direction() ) std::swap(p1,p2);
    // now order p2 after p1 on r underlying line 
    typename R::Position_on_line_d pos;
    pos(p1, _r.point(0),_r.point(1), l1);
    pos(p2, _r.point(0),_r.point(1), l2);
    if ( l1 < FT(0) ) p1 = _r.point(0);
    if ( l2 < FT(0) ) { return _result = NO; }
    if ( p1 == p2 ) { _ip = p1; return _result = NO; }
    _is = Segment_d(p1,p2);
    return _result = SEGMENT;
  }


  if (res == Int_obj_type::POINT && 
      FT(0) <= l1 && FT(0) <= l2 && l2 <= FT(1)) 
  { return _result = POINT; }
  return _result = NO; 
  // now not parallel
}

template <class R>
bool
Ray_d_Segment_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool
Ray_d_Segment_d_pairH<R>::intersection(Segment_d& s)
{ if (!_known) intersection_type();
  if (_result != SEGMENT) return false;
  s = _is; return true;
}


template <class R>
class Segment_d_Segment_d_pairH {
public:
  enum Intersection_result { NO, POINT, SEGMENT };
  typedef typename R::FT FT;
  typedef typename R::Point_d Point_d;
  typedef typename R::Segment_d Segment_d;

protected:
  Segment_d _s1,_s2;
  bool _known;
  Intersection_result _result;
  Point_d _ip; Segment_d _is;  
public:
  Segment_d_Segment_d_pairH() : _known(false) {}
  Segment_d_Segment_d_pairH(const Segment_d& s1, const Segment_d& s2)
    : _s1(s1), _s2(s2), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Segment_d& result);
};

template <class R>
Segment_d_Segment_d_pairH<R>::Intersection_result
Segment_d_Segment_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  if ( _s1.is_degenerate() ) {
    if ( _s2.is_degenerate() ) {
      if ( _s1.point(0) == _s2.point(0) ) 
      { _ip = _s1.point(0); return _result = POINT; }
      else { return _result = NO; }
    } else {
      if ( _s2.has_on(_s1.point(0)) ) 
      { _ip = _s1.point(0); return _result = POINT; }
      else { return _result = NO; }
    }
  }
  if ( _s2.is_degenerate() ) { 
    CGAL_assertion( !_s1.is_degenerate() );
    if ( _s1.has_on(_s2.point(0)) ) 
    { _ip = _s2.point(0); return _result = POINT; }
    else { return _result = NO; }
  }  

  // now s1,s2 not degenerate
  typedef Line_line_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT l1,l2; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_s1.point(0),_s1.point(1),
              _s2.point(0),_s2.point(1),_ip,l1,l2);

  if (res == Int_obj_type::LINE)   
  { Point_d p1 = _s1.min(), p2 = _s1.max();
    Point_d q1 = _s2.min(), q2 = _s2.max();
    Point_d s,t;
    // now order the for points along the line
    typename R::Position_on_line_d pos;
    pos(p1, q1, q2, l1); pos(p2, q1, q2, l2);
    if ( l1 < FT(0) ) {
      if ( l2 < FT(0) ) { return _result = NO; }
      else { s = q1; }
    } else { s = p1; } // l1 >= 0
    if ( l2 > FT(1) ) {
      if ( l1 > FT(1) ) { return _result = NO; }
      else { t = q2; }
    } else { t = p2; } // l2 <= 1
    if ( s == t ) { _ip = s; return _result = POINT; }
    _is = Segment_d(s,t); return _result = SEGMENT;
  }


  if ( res == Int_obj_type::POINT && 
       FT(0) <= l1 && l1 <= FT(1) && 
       FT(0) <= l2 && l2 <= FT(1) ) 
  { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool
Segment_d_Segment_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool
Segment_d_Segment_d_pairH<R>::intersection(Segment_d& s)
{ if (!_known) intersection_type();
  if (_result != SEGMENT) return false;
  s = _is; return true;
}


template <class R>
class Line_d_Hyperplane_d_pairH {
public:
  enum Intersection_result {NO, POINT, LINE};
  typedef typename R::Point_d      Point_d;
  typedef typename R::Line_d       Line_d;
  typedef typename R::Hyperplane_d Hyperplane_d;
  typedef typename R::FT           FT;

protected:
  Line_d _l; Hyperplane_d _h;
  bool _known;
  Intersection_result _result;
  Point_d _ip;
public:
  Line_d_Hyperplane_d_pairH() : _known(false) {}
  Line_d_Hyperplane_d_pairH(const Line_d& l, const Hyperplane_d& h) 
   : _l(l), _h(h), _known(false) {}

  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Line_d& result);
};

template <class R>
Line_d_Hyperplane_d_pairH<R>::Intersection_result
Line_d_Hyperplane_d_pairH<R>::intersection_type()
{ if (_known) return _result;
  _known = true;

  // CGAL_assertion_msg((!_l.is_degenerate())); 
  typedef Line_hyperplane_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT lambda; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_l.point(0),_l.point(1),_h,_ip,lambda);

  if (res == Int_obj_type::LINE)  { return _result = LINE; }
  if (res == Int_obj_type::POINT) { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool
Line_d_Hyperplane_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool
Line_d_Hyperplane_d_pairH<R>::intersection(Line_d& l)
{ if (!_known) intersection_type();
  if (_result != LINE) return false;
  l = _l; return true;
}

template <class R>
class Ray_d_Hyperplane_d_pairH {
public:
  enum Intersection_result {NO, POINT, RAY};
  typedef typename R::FT FT;
  typedef typename R::Point_d Point_d;
  typedef typename R::Ray_d Ray_d;
  typedef typename R::Hyperplane_d Hyperplane_d;
protected:
  Ray_d _r; Hyperplane_d _h;
  bool                    _known;
  Intersection_result     _result;
  Point_d                 _ip;  
public:
  Ray_d_Hyperplane_d_pairH() : _known(false) {}
  Ray_d_Hyperplane_d_pairH(const Ray_d& r, const Hyperplane_d& h) 
   : _r(r), _h(h), _known(false) {}
  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Ray_d& result);
};

template <class R>
Ray_d_Hyperplane_d_pairH<R>::Intersection_result
Ray_d_Hyperplane_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  // CGAL_assertion( !_r.is_degenerate() );
  typedef Line_hyperplane_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT lambda; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_r.point(0),_r.point(1),_h,_ip,lambda);

  if ( res == Int_obj_type::POINT && FT(0) <= lambda ) 
  { return _result = POINT; }
  if ( res == Int_obj_type::LINE )  
  { return _result = RAY; }
  return _result = NO; 
}

template <class R>
bool
Ray_d_Hyperplane_d_pairH<R>::intersection(Point_d& p)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  p = _ip; return true;
}

template <class R>
bool
Ray_d_Hyperplane_d_pairH<R>::intersection(Ray_d& r)
{ if (!_known) intersection_type();
  if (_result != RAY) return false;
  r = _r; return true;
}

template <class R>
class Segment_d_Hyperplane_d_pairH {
public:
  enum Intersection_result {NO, POINT, SEGMENT};
  typedef typename R::FT           FT;
  typedef typename R::Point_d Point_d;
  typedef typename R::Segment_d Segment_d;
  typedef typename R::Hyperplane_d Hyperplane_d;

protected:
  Segment_d    _s;
  Hyperplane_d _h;
  bool                    _known;
  Intersection_result     _result;
  Point_d                 _ip;  
public:
  /*{\Mcreation 4}*/
  Segment_d_Hyperplane_d_pairH() : _known(false) {}
  Segment_d_Hyperplane_d_pairH(
    const Segment_d& s, const Hyperplane_d& h) :
    _s(s), _h(h), _known(false) {}

  Intersection_result intersection_type();
  bool intersection(Point_d& result);
  bool intersection(Segment_d& result);
};

template <class R>
Segment_d_Hyperplane_d_pairH<R>::Intersection_result
Segment_d_Hyperplane_d_pairH<R>::intersection_type()
{ 
  if (_known) return _result;
  _known = true;

  if ( _s.is_degenerate() ) 
    if ( _h.has_on(_s.point(0)) ) 
    { _ip = _s.point(0); return _result = POINT; } 
    else { return _result = NO; }

  typedef Line_hyperplane_intersectionHd<R> Int_obj_type;
  Int_obj_type Intersect;
  FT lambda; 
  typename Int_obj_type::Intersection_result res = 
    Intersect(_s.point(0),_s.point(1),_h,_ip,lambda);

  if ( res == Int_obj_type::LINE )
  { return _result = SEGMENT; }
  if ( res == Int_obj_type::POINT && 
       FT(0) <= lambda && lambda <= FT(1) ) 
  { return _result = POINT; }
  return _result = NO; 
}

template <class R>
bool
Segment_d_Hyperplane_d_pairH<R>::intersection(Point_d& pt)
{ if (!_known) intersection_type();
  if (_result != POINT) return false;
  pt = _ip; return true;
}

template <class R>
bool
Segment_d_Hyperplane_d_pairH<R>::intersection(Segment_d& s)
{ if (!_known) intersection_type();
  if (_result != LINE) return false;
  s = _s; return true;
}



CGAL_END_NAMESPACE
#endif //CGAL_INTERSECTION_OBJECTSHD_H

