#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H

CGAL_BEGIN_NAMESPACE

template <class pR>
class Point_d : public pR::Point_d_base
{ public:
  typedef typename pR::Point_d_base Base;
  typedef Point_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Point_d(int d=0) : Base(d) {}
  Point_d(int d, const Origin &o) : Base(d,o) {}

  Point_d(int a, int b, int c = 1) : Base(a,b,c) {} 
  Point_d(const RT& a, const RT& b, const RT& c = 1) :
    Base(a,b,c) {}  
  Point_d(int a, int b, int c, int d) : Base(a,b,c,d) {}
  Point_d(const RT& a, const RT& b, const RT& c, const RT& d) :
    Base(a,b,c,d) {}

  template <class InputIterator>
  Point_d (int d, InputIterator first, InputIterator last)
    : Base (d, first, last) {}
  template <class InputIterator>
  Point_d(int d, InputIterator first, InputIterator last, const RT& D)
    : Base (d, first, last, D) {}
  Point_d(const Self &p) : Base(p) {}
  Point_d(const Base& p) : Base(p) {}
 
  Vector_d<R> operator-(const Origin& o) const 
  { return Base::operator-(o); }
  Vector_d<R> operator-(const Self& q) const
  { return Base::operator-(q); }
  Self operator+(const Vector_d<R>& v) const
  { return Base::operator+(v); }
  Self operator-(const Vector_d<R>& v) const
  { return Base::operator-(v); }
  Self& operator+=(const Vector_d<R>& v) 
  { return static_cast<Self&>(Base::operator+=(v)); }
  Self& operator-=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator-=(v)); }
  
};

CGAL_END_NAMESPACE
#endif //CGAL_POINT_D_H
