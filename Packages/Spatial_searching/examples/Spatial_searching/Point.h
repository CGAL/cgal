
class Construct_coord_iterator;

struct Point {
  double vec[3];

  Point() { vec[0]= vec[1] = vec[2] = 0; }
  Point (double x, double y, double z) { vec[0]=x; vec[1]=y; vec[2]=z;  }
 
  double x() const { return vec[ 0 ]; }
  double y() const { return vec[ 1 ]; }
  double z() const { return vec[ 2 ]; }

  double& x() { return vec[ 0 ]; }
  double& y() { return vec[ 1 ]; }
  double& z() { return vec[ 2 ]; }
  
  bool operator==(const Point& p) const 
  {
    return (x() == p.x()) && (y() == p.y()) && (z() == p.z())  ;
  }

  bool  operator!=(const Point& p) const { return ! (*this == p); }
}; //end of class


namespace CGAL {

  template <>
  struct Kernel_traits<::Point> {
    struct Kernel {
      typedef double FT;
      typedef double RT;
    };
  };
}

class Construct_coord_iterator {
public:
  const double* operator()(const Point& p)
  {
    return static_cast<const double*>(p.vec);
  }

  const double* operator()(const Point& p, int)
  {
    return static_cast<const double*>(p.vec+3);
  }
};





class Distance
{
public:

  typedef Point Query_item;

  double distance(const Point& p1, const Point& p2) const
  {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  template <class TreeTraits>
  double min_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const 
  {   
    double distance(0.0);
    double h;
    h=p.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }
  
  template <class TreeTraits>
  double max_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const
  {   
    double distance(0.0);
    double h;
    h=p.x();
    if (h >= (b.min_coord(0)+b.max_coord(0))/2.0) 
      distance += (h-b.min_coord(0))*(h-b.min_coord(0)); 
    else
      distance += (b.max_coord(0)-h)*(b.max_coord(0)-h);
    h=p.y();
    if (h >= (b.min_coord(1)+b.max_coord(1))/2.0) 
      distance += (h-b.min_coord(1))*(h-b.min_coord(1)); 
    else
      distance += (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    if (h >= (b.min_coord(2)+b.max_coord(2))/2.0) 
      distance += (h-b.min_coord(2))*(h-b.min_coord(2)); 
    else
      distance += (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return distance;
  }
  
  double new_distance(double& dist, double old_off, double new_off,
		      int cutting_dimension)  const {
    return dist + new_off*new_off - old_off*old_off;
  }
  
  double transformed_distance(double d) const {
    return d*d;
  }
  
  double inverse_of_transformed_distance(double d) const {
    return sqrt(d);
  }
  
}; // end of class
