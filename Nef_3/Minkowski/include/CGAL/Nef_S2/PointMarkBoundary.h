#ifndef CGAL_POINT_MARK_BOUNDARY_H
#define CGAL_POINT_MARK_BOUNDARY_H

namespace CGAL {

template <class K>
class PointMarkBoundary {

  typedef typename K::Point_3 Point_3;
  typedef PointMarkBoundary<K> Self;

  Point_3 p;
  bool b;
  bool boundary;

 public:
  PointMark() : p(0,0,0), b(true), boundary(true) {}
  PointMark(const Self& pm) { p = pm.p; b = pm.b; }
  PointMark(const Point_3& p_, bool b_, bool boundary_) 
    : p(p_), b(b_) boundary(boundary_) {}

  Self& operator=(const Self& pm) {
    p = pm.p;
    b = pm.b;
    boundary = pm.boundary;
    return *this;
  }

  Self& operator+=(const Self& pm) {
    p = p + (pm.p - CGAL::ORIGIN);
    b = b && pm.b;
    boundary = boundary && pm.boundary;
    return *this;
  }

  Point_3 point() const {
    return p;
  }

  bool boolean() const {
    return b;
  }

  bool on_boundary() const {
    return boundary;
  }

  void set_boolean(bool b_) {
    b = b_;
  }

  void set_boundary(bool b_) {
    boundary = b;
  }

};

template <typename Kernel>
std::ostream& operator<<(std::ostream& out, 
			 const PointMark<Kernel>& pm) {
  out << pm.point() << "/" << pm.boolean();
  return out;
}

template <typename Kernel>
bool operator==(const PointMark<Kernel>& pm1,
		const PointMark<Kernel>& pm2) {
  return 
    pm1.point() == pm2.point() &&
    pm1.boolean() == pm2.boolean();
}

template <typename Kernel>
const PointMark<Kernel> operator+(const PointMark<Kernel>& pm1,
				  const PointMark<Kernel>& pm2) {
  PointMark<Kernel> ret(pm1);
  ret += pm2;
  return ret;
}

} //namespace CGAL
#endif
