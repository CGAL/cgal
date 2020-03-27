struct Point {
  double vec[2];

  Point() { vec[0]= vec[1] = 0; }
  Point (double x, double y) { vec[0]=x; vec[1]=y; }

  double x() const { return vec[ 0 ]; }
  double y() const { return vec[ 1 ]; }

  double& x() { return vec[ 0 ]; }
  double& y() { return vec[ 1 ]; }

  bool operator==(const Point& p) const
  {
    return (x() == p.x()) && (y() == p.y())  ;
  }

  bool  operator!=(const Point& p) const { return ! (*this == p); }
}; //end of class



namespace CGAL {

  template <>
  struct Kernel_traits<Point> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}


struct Construct_coord_iterator {
  typedef const double* result_type;
  const double* operator()(const Point& p) const
  { return static_cast<const double*>(p.vec); }

  const double* operator()(const Point& p, int) const
  { return static_cast<const double*>(p.vec+2); }


};
