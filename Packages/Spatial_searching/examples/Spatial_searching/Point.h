
class Construct_coord_iterator;

class Point
{
public:

  friend Construct_coord_iterator;

  class R
  { 
  public:
    typedef double FT;
  };

  private:
  double   vec[ 3 ];
  
public: 
  
  Point()
  { 
    for  ( int ind = 0; ind < 3; ind++ )
      vec[ ind ] = 0;
  }

  Point (double& x, double& y, double& z)
  {
    vec[0]=x;
    vec[1]=y;
    vec[2]=z;
  }

  inline
  int dimension() const
  {
    return  3;
  }
 
  inline
  double x() const
  { 
	return vec[ 0 ];
  }

  inline
  double y() const
  { 
 	return vec[ 1 ];
  }
  
  inline
  double z() const
  { 
	return vec[ 2 ];
  }

  inline
  void set_coord(int k, double x)
  {
    vec[ k ] = x;
  }

  /*
  inline
  double  & operator[](int k)  
  {
    return  vec[ k ];
  }

  inline
  double  operator[](int k) const
  {
    return  vec[ k ];
  }
  */
}; //end of class


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


inline
bool
operator!=(const Point& p, const Point& q)
{
  return ( (p.x() != q.x()) || (p.y() != q.y()) || (p.z() != q.z()) ); 
}

inline
bool
operator==(const Point& p, const Point& q)
{
  return ( (p.x() == q.x()) && (p.y() == q.y()) && (p.z() == q.z()) ) ;
}


class Point3D_distance
{
public:

double distance(const Point& p1, const Point& p2) const
{
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
}

template <class TreeTraits>
double min_distance_to_queryitem(const Point& p,
                                        const CGAL::Kd_tree_rectangle<TreeTraits>& b) const 
{   double distance(0.0);
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
double max_distance_to_queryitem(const Point& p,
                                        const CGAL::Kd_tree_rectangle<TreeTraits>& b) const
{   double distance(0.0);
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
