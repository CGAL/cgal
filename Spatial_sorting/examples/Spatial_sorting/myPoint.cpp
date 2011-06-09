#include <CGAL/spatial_sort.h>

struct MyPoint {
  double x,y;
  int color;
  MyPoint()
    : x(0), y(0),color(0)
  {}

  MyPoint(double x, double y, int color=0)
    : x(x), y(y), color(color)
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
  std::vector< MyPoint > points;

  points.push_back(MyPoint(14,12, 3));
  points.push_back(MyPoint(1,2  , 0));
  points.push_back(MyPoint(414,2, 5));
  points.push_back(MyPoint(4,21 , 1));
  points.push_back(MyPoint(7,74 , 2));
  points.push_back(MyPoint(74,4 , 4));  
  
  MySpatialSortingTraits sst;
  CGAL::spatial_sort(points.begin(), points.end(), sst);

  for (std::vector< MyPoint >::iterator it=points.begin();it!=points.end();++it)
    std::cout << it->color << " ";
  std::cout << "\n";  
  
  std::cerr << "done" << std::endl;
  return 0;
}
