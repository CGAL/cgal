
class PVector;
class segment;
class triangle;
class dir;
class ray;


class point{
public:
  friend PVector;
  point()
    : _x(0), _y(0)
  {}

  point(const point& p)
    : _x(p._x), _y(p._y)
  {}

  point(double x, double y)
    : _x(x), _y(y)
  {}

  double x() const
  {
    return _x;
  }

  double y() const
  {
    return _y;
  }

  bool operator==(const point& p) const
  {
    return x() == p.x() && y() == p.y();
  }
  
  bool operator!=(const point& p) const
  {
    return x() != p.x() || y() != p.y();
  }
  
private:
  double _x, _y;
};



ostream& operator<<(ostream& os, const point& p)
{
  os << "(" << p.x() << ", " << p.y() << ")";
  return os;
}



class dir
{
protected:
  double _x, _y;
public:
  dir() {}
  dir(double x, double y) : _x(x), _y(y) {}
  dir(const point* p, const point* q)
    : _x( q->x()-p->x() ), _y( q->y()-p->y() ) {}

  double x() const { return _x; }
  double y() const { return _y; }
  
  dir perpendicular() const { return dir(y(),-x()); }
};


ostream& operator<< (ostream& os, const dir& d)
{
  os << "(" << d.x() << ", " << d.y() << ")";
  return os;
}



class segment{
public:
  segment(point* s, point* t)
    : _s(s), _t(t)
  {}

  segment()
    : _s(0), _t(0)
  {}

  segment(const segment& s)
    : _s(s._s), _t(s._t)
  {}

  point* source() const
  {
    return _s;
  }

  point* target() const
  {
    return _t;
  }

  bool operator==(const segment& s) const
  {
    return *source() == *s.source() && *target() == *s.target();
  }
  
  bool operator!=(const segment& s) const
  {
    return *source() != *s.source() || *target() != *s.target();
  }
  
  dir direction() const
  {
    return dir(source(),target());
  }

private:
  point *_s, *_t;
};

ostream& operator<<(ostream& os, const segment& s)
{
  os << "(" << *(s.source()) << ", " << *(s.target()) << ")";
  return os;
}




class triangle{
public:
  triangle(point* p, point* q, point* r)
    : _p(p), _q(q), _r(r)
  {}

  triangle()
    : _p(0), _q(0), _r(0)
  {}

  triangle(const triangle& t)
    : _p(t._p), _q(t._q), _r(t._r)
  {}

  point* vertex(int i) const
  {
    switch(i){
    case 0:
      return _p;
    case 1:
      return _q;
    }
    return _r;
  }


  bool operator==(const triangle& t) const
  {
    return vertex(0) == t.vertex(0) && vertex(1) == t.vertex(1) && vertex(2) == t.vertex(2) ;
  }
  
  bool operator!=(const triangle& t) const
  {
    return vertex(0) != t.vertex(0) || vertex(1) != t.vertex(1) || vertex(2) != t.vertex(2) ;
  }
  
private:
  point *_p, *_q, *_r;
};

ostream& operator<<(ostream& os, const triangle& t)
{
  os << "(" << t.vertex(0) << ", " << t.vertex(1) << ", " << t.vertex(2) << ")";
  return os;
}


class ray
{
friend ostream& operator<< (ostream&, const ray&);

protected:
  point _p;
  dir _d;
public:
  ray() {}
  ray(const point& p, const dir& d) : _p(p), _d(d) {}

  point p() const { return _p; }
  dir   d() const { return _d; }
};

ostream& operator<<(ostream& os, const ray& r)
{
  os << r.p() << "+" << r.d();
  return os;
}



class PVector {
public:
  PVector()
    : V(NULL), _size(0)
  {}

  ~PVector()
  {
    if(V == NULL){
      delete []V;
    }
  }

  int size()
  {
    return _size;
  }

  point* operator[](int i) const
  {
    return V+i;
  }

  void reserve(int n)
  {
    if(V == NULL){
      delete []V;
    }
    V = new point[n];
    _size = n;
  }

  void read(istream& is)
  {
    int n;
    is >> n;

    reserve(n);
    double x,y;
    point* ptr = V;
    for(;n>0;n--){
      is >> x >> y;
      ptr->_x = x;
      ptr->_y = y;
      ++ptr;
    }
  }

private:
  point* V;
  int _size;
};

istream& operator>>(istream& is, PVector &pv)
{
  pv.read(is);
  return is;
}


#include <CGAL/Cartesian.h>
#include <CGAL/basic_constructionsC2.h>
#include <CGAL/predicates_on_pointsC2.h>



class Euclidean_2 {
public:
  typedef point*  Point;
  typedef segment Segment;
  typedef triangle Triangle;
  typedef dir Direction;
  typedef ray Ray;

 
  bool compare(const Point &p, const Point &q) const
  {
    return (p == q);
  }

  CGAL::Comparison_result compare_x(const Point &p, const Point &q) const
  {
    return CGAL::compare(p->x(), q->x());
  }

  CGAL::Comparison_result compare_y(const Point &p, const Point &q) const
  {
    return CGAL::compare(p->y(), q->y());
  }

  CGAL::Orientation orientation(const Point &p,
                               const Point &q,
                               const Point &r) const
  {
    if(*p==*q || *p == *r || *q == *r){
      cout << "coll" << endl;
      return CGAL::COLLINEAR;
    }
    return CGAL::orientationC2(p->x(), p->y(), q->x(), q->y(), r->x(), r->y());
  }


  CGAL::Orientation extremal(const Point &p,
                            const Point &q,
			    const Point &r) const
  {
    if(*p==*q || *p == *r || *q == *r){
      cout << "coll" << endl;
      return CGAL::COLLINEAR;
    }
    return CGAL::orientationC2(p->x(), p->y(), q->x(), q->y(), r->x(), r->y());
  }

  point circumcenter(const point& p, const point& q, const point& r)
  {
    double px( p.x());
    double py( p.y());
    double qx( q.x());
    double qy( q.y());
    double rx( r.x());
    double ry( r.y());

    double px_qx( px - qx);
    double py_qy( py - qy);
    double qx_rx( qx - rx);
    double qy_ry( qy - ry);
    double rx_px( rx - px);
    double ry_py( ry - py);

    double px2_py2( px*px + py*py);
    double qx2_qy2( qx*qx + qy*qy);
    double rx2_ry2( rx*rx + ry*ry);

    double num_x( px2_py2*qy_ry + qx2_qy2*ry_py + rx2_ry2*py_qy);
    double num_y( px2_py2*qx_rx + qx2_qy2*rx_px + rx2_ry2*px_qx);

    double den_x( ( px*qy_ry + qx*ry_py + rx*py_qy) * 2.0);
    double den_y( ( py*qx_rx + qy*rx_px + ry*px_qx) * 2.0);

    return point(num_x/den_x, num_y/den_y);
  }
};

