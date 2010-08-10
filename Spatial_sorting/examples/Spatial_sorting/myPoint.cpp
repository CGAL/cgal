#include <CGAL/spatial_sort.h>

struct MyPoint {
  double x,y;

  MyPoint()
    : x(0), y(0)
  {}

  MyPoint(double x, double y)
    : x(x), y(y)
  {}
};


struct MyLessX {

  bool operator()(const MyPoint& p, const MyPoint& q) const
  {
    return p.x < q.x;
  }

};

struct MyLessY {

  bool operator()(const MyPoint& p, const MyPoint& q) const
  {
    return p.y < q.y;
  }

};

struct MySpatialSortingTraits {

  typedef MyPoint Point_2;

  typedef MyLessX Less_x_2;
  typedef MyLessY Less_y_2;
  
  Less_x_2 less_x_2_object() const
  {
    return Less_x_2();
  }

  Less_y_2 less_y_2_object() const
  {
    return Less_y_2();
  }
};

int main()
{
  MyPoint points[2];

  points[0] = MyPoint(78,12);
  points[1] = MyPoint(3,121);
  MySpatialSortingTraits sst;
  CGAL::spatial_sort(points, points+2, sst);
  std::cerr << "done" << std::endl;
  return 0;
}
