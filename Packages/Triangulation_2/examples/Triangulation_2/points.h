#include <CGAL/Cartesian.h>
#include <CGAL/predicates/kernel_ftC2.h>


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



std::ostream& operator<<(std::ostream& os, const point& p)
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


std::ostream& operator<< (std::ostream& os, const dir& d)
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

std::ostream& operator<<(std::ostream& os, const segment& s)
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

std::ostream& operator<<(std::ostream& os, const triangle& t)
{
  os << "(" << t.vertex(0) << ", " << t.vertex(1) << ", " << t.vertex(2) << ")";
  return os;
}


class ray
{
friend std::ostream& operator<< (std::ostream&, const ray&);

protected:
  point _p;
  dir _d;
public:
  ray() {}
  ray(const point& p, const dir& d) : _p(p), _d(d) {}

  point p() const { return _p; }
  dir   d() const { return _d; }
};

std::ostream& operator<<(std::ostream& os, const ray& r)
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

  void read(std::istream& is)
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

std::istream& operator>>(std::istream& is, PVector &pv)
{
  pv.read(is);
  return is;
}

class compare_x_2{
public:
  typedef point*  Point;
  CGAL::Comparison_result operator()( const Point& p, const Point&  q) const 
    {
      if (p->x() < q->x()) return CGAL::SMALLER;
      else if (p->x() > q->x()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }
};


class compare_y_2{
public:
  typedef point*  Point;
  CGAL::Comparison_result operator()( const Point& p, const Point&  q) const 
    {
      if (p->y() < q->y()) return CGAL::SMALLER;
      else if (p->y() > q->y()) return CGAL::LARGER;
      else return CGAL::EQUAL;
    }
};

class orientation_2
{
public:
  typedef point*  Point;
  CGAL::Orientation
  operator()( const Point& p, const Point&  q, const Point& r) const 
    {
      if(*p==*q || *p == *r || *q == *r){
	std::cout << "coll" << std::endl;
	return CGAL::COLLINEAR;
      }
      return CGAL::orientationC2(p->x(), p->y(), 
				 q->x(), q->y(), 
				 r->x(), r->y());
    }
};








class Euclidean_2 {
public:
  typedef point*  Point_2;
  typedef segment Segment_2;
  typedef triangle Triangle_2;

  typedef compare_x_2    Compare_x_2;
  typedef compare_y_2    Compare_y_2;
  typedef orientation_2  Orientation_2; 

  Compare_x_2
  compare_x_2_object() const
    { return Compare_x_2();}

  Compare_y_2
  compare_y_2_object() const
    { return Compare_y_2();}

  Orientation_2
  orientation_2_object() const
    { return Orientation_2();}

};

