#ifndef CHR_MOEBIUS_POINT_H
#define CHR_MOEBIUS_POINT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class Point, class Weight>
class Moebius_point : public Point
{
 public:
  Moebius_point (const Point &p = Point (),
		  const Weight &lambda = Weight (1),
		  const Weight &mu = Weight (0))
    : Point (p), _lambda (lambda), _mu (mu) {}

  Weight lambda () const {  return _lambda; }
  Weight mu () const {  return _mu; }

 private:
  Weight _lambda, _mu;
};

template <class P, class W>
std::ostream & operator<< (std::ostream & os,
			   const Moebius_point<P,W> &p)
{
  return os << (P) p << ' ' << p.lambda() << ' ' << p.mu();
}

template <class P, class W>
std::istream & operator>> (std::istream & is,
			   Moebius_point<P,W> &p)
{
  W lambda, mu;
  P _p;
  is >> _p >> lambda >> mu;
  p = Moebius_point<P,W> (_p, lambda, mu);
  return is;
}



CGAL_END_NAMESPACE



#endif// CHR_MOEBIUS_POINT_H
