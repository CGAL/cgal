#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Fuzzy_iso_box_d.h>

#include <vector>  
#include <iostream>

// create own Point type

template <class R_> 
class Point
{
public:
    typedef R_ R;
    typedef typename R::FT FT;

private:
  FT   vec[ 3 ];
  
public: 
  
  Point()
  { 
    for  ( int ind = 0; ind < 3; ind++ )
      vec[ ind ] = 0;
  }

  Point (FT& x, FT& y, FT& z)
  {
    vec[0]=x;
    vec[1]=y;
    vec[2]=z;
  }

  int dimension() const
  {
    return  3;
  }
 
  FT x() const
  { 
    return vec[ 0 ];
  }

  FT y() const
  { 
    return vec[ 1 ];
  }
  
  FT z() const
  { 
    return vec[ 2 ];
  }

  void set_coord(int k, FT x)
  {
    vec[ k ] = x;
  }
  

  FT  & operator[](int k)  
  {
    return  vec[ k ];
  }

  FT  operator[](int k) const
  {
    return  vec[ k ];
  }
}; //end of class

class R { // define representation class
public:
	typedef double FT;
	typedef double RT;
	typedef Point<R> Point_3;
}; // end of class

inline
bool
operator!=(const Point<R>& p, const Point<R>& q)
{
  return ( (p[0] != q[0]) || (p[1] != q[1]) || (p[2] != q[2]) ); 
}

inline
bool
operator==(const Point<R>& p, const Point<R>& q)
{
  return ( (p[0] == q[0]) && (p[1] == q[1]) && (p[2] == q[2]) ) ;
}

// not essential by specification but nice to have
std::ostream &operator<<(std::ostream &os, const Point<R> &p)
{
  std::cout << "(";
  for(int i = 0; i < 3; i++)
    {
      std::cout << p[i] ;
      if (i < p.dimension() - 1) std::cout << ", ";
    }
  std::cout << ")";
  return os;
}

// create own distance class

template <class Point>
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

double min_distance_to_queryitem(const Point& q,
                                        const CGAL::Kd_tree_rectangle<double>& b) const 
{   double distance(0.0);
    double h;
    h=q.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=q.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=q.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
}

double max_distance_to_queryitem(const Point& q,
                                        const CGAL::Kd_tree_rectangle<double>& b) const 
{   double distance(0.0);
    double h;
    h=q.x();
    if (h >= (b.min_coord(0)+b.max_coord(0))/2.0) 
                distance += (h-b.min_coord(0))*(h-b.min_coord(0)); 
	  else
                distance += (b.max_coord(0)-h)*(b.max_coord(0)-h);
    h=q.y();
    if (h >= (b.min_coord(1)+b.max_coord(1))/2.0) 
                distance += (h-b.min_coord(1))*(h-b.min_coord(1)); 
	  else
                distance += (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=q.z();
    if (h >= (b.min_coord(2)+b.max_coord(2))/2.0) 
                distance += (h-b.min_coord(2))*(h-b.min_coord(2)); 
	  else
                distance += (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return distance;
}

double new_distance(double& dist, double old_off, double new_off,
                int cutting_dimension) const {
                return dist + new_off*new_off - old_off*old_off;
}

double transformed_distance(double d) const {
        return d*d;
}

double inverse_of_transformed_distance(double d) const {
        return sqrt(d);
}

}; // end of class

typedef R::Point_3 My_point;
typedef CGAL::Creator_uniform_3<double,My_point> Creator;
typedef CGAL::Kd_tree_traits_point<My_point> TreeTraits;
typedef CGAL::Orthogonal_standard_search<TreeTraits, Point3D_distance<My_point> > 
NN_orthogonal_search;

typedef std::vector<TreeTraits::Point> Vector;

int main() {

  int bucket_size=10;
  
  const int data_point_number=1000;
  
  typedef std::list<My_point> Point_list;
  Point_list data_points, res;

  // generate random data points  
  CGAL::Random_points_in_cube_3<My_point,Creator> g( 1.0);
  CGAL::copy_n( g, data_point_number, std::back_inserter(data_points));
  
  TreeTraits tr(bucket_size);
  typedef CGAL::Kd_tree<TreeTraits> Tree;
  Tree d(data_points.begin(), data_points.end(), tr);

  // generate random query points
  const int query_point_number=5;
  CGAL::Random_points_in_cube_3<My_point,Creator> h( 1.0);
  Vector query_points;
  CGAL::copy_n(h, query_point_number, std::back_inserter(query_points));

  Point3D_distance<My_point> tr_dist;

  // nearest neighbor searching 
  std::vector<NN_orthogonal_search::Point_with_distance> 
  the_nearest_neighbors;
  
  // furthest neighbor searching 
  std::vector<NN_orthogonal_search::Point_with_distance> 
  the_furthest_neighbors;
  
  for (int i=0; i < query_point_number; i++) { 
     // nearest neighbour searching
     NN_orthogonal_search NN1(d, query_points[i], tr_dist);
     NN1.the_k_neighbors(std::back_inserter(the_nearest_neighbors));
     
     // furthest neighbour searching
     NN_orthogonal_search NN2(d, query_points[i], tr_dist, 1, 0.0, false);
     NN2.the_k_neighbors(std::back_inserter(the_furthest_neighbors));
  }
  
  std::cout << "results neighbor searching:" << std::endl;

  for (int j=0; j < query_point_number; j++) { 
     std::cout << " d(q, nearest neighbor)=  " << 
     tr_dist.inverse_of_transformed_distance(the_nearest_neighbors[j].second) << 
     "    d(q, furthest neighbor)= " << 
     tr_dist.inverse_of_transformed_distance(the_furthest_neighbors[j].second) << std::endl; 
  } 

  return 0;
};

