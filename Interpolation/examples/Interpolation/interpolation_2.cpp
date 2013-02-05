// compares the result of several interpolation methods

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Origin.h>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                         Coord_type;
typedef K::Vector_2                                   Vector;
typedef K::Point_2                                    Point;

typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
typedef CGAL::Interpolation_traits_2<K>               Traits;

typedef std::vector< std::pair<Point, Coord_type> >   Coordinate_vector;
typedef std::map<Point, Coord_type, K::Less_xy_2>     Point_value_map;
typedef std::map<Point,  Vector, K::Less_xy_2 >       Point_vector_map;

////////////////////////////////

int main()
{
  //number of sample points:
  int n = 24;
  //number of interpolation points:
  int m = 20;

  std::vector<Point> points;
  points.reserve(n+m);

  //put four bounding box points:
  points.push_back(Point(-3,-3));
  points.push_back(Point(3,-3));
  points.push_back(Point(-3,3));
  points.push_back(Point(3,3));

  // Create n+m-4 points within a disc of radius 2
  double r_d = 3;
  CGAL::Random_points_in_disc_2<Point> g(r_d );
  CGAL::cpp11::copy_n( g, n+m, std::back_inserter(points));

  Delaunay_triangulation T;


  Point_value_map values;
  Point_vector_map gradients;


  //parameters for quadratic function:
  Coord_type alpha = Coord_type(1.0),
    beta1 = Coord_type(2.0),
    beta2 = Coord_type(1.0),
    gamma1 = Coord_type(0.3),
    gamma2 = Coord_type(0.0),
    gamma3 = Coord_type(0.0),
    gamma4 = Coord_type(0.3);

  for(int j=0; j<n ; j++){
    T.insert(points[j]);

    //determine function value/gradient:
    Coord_type x(points[j].x());
    Coord_type y(points[j].y());

    Coord_type value = alpha + beta1*x + beta2*y + gamma1*(x*x) +
      gamma4*(y*y) + (gamma2+ gamma3) *(x*y);
    Vector gradient(beta1+ (gamma2+ gamma3)*y + Coord_type(2)*(gamma1*x),
		    beta2+ (gamma2+ gamma3)*x + Coord_type(2)*(gamma4*y));
    values.insert(std::make_pair(points[j], value));
    gradients.insert(std::make_pair(points[j], gradient));
  }

  //variables for statistics:
  std::pair<Coord_type, bool> res;
  Coord_type error, l_total = Coord_type(0),
    q_total(l_total), f_total(l_total), s_total(l_total),
    ssquare_total(l_total), l_max(l_total),
    q_max(l_total), f_max(l_total), s_max(l_total),
    ssquare_max(l_total),
    total_value(l_total), l_value(l_total);
  int failure(0);

  //interpolation + error statistics
  for(int i=n;i<n+m;i++){
    Coord_type x(points[i].x());
    Coord_type y(points[i].y());

    Coord_type exact_value = alpha + beta1*x + beta2*y + gamma1*(x*x) +
      gamma4*(y*y) + (gamma2+ gamma3) *(x*y);

    total_value += exact_value;

    //Coordinate_vector:
    std::vector< std::pair< Point, Coord_type > > coords;
    Coord_type norm =
      CGAL::natural_neighbor_coordinates_2(T, points[i],
					   std::back_inserter(coords)).second;

    assert(norm>0);

    //linear interpolant:
    l_value =  CGAL::linear_interpolation(coords.begin(), coords.end(),
					  norm,
					  CGAL::Data_access<Point_value_map>(values));

    error = CGAL_NTS abs(l_value - exact_value);
    l_total += error;
    if (error > l_max) l_max = error;

    //Farin interpolant:
    res =  CGAL::farin_c1_interpolation(coords.begin(),
					coords.end(), norm,points[i],
					CGAL::Data_access<Point_value_map>(values),
					CGAL::Data_access<Point_vector_map>
					(gradients),
					Traits());
    if(res.second){
      error = CGAL_NTS abs(res.first - exact_value);
      f_total += error;
      if (error > f_max) f_max = error;
    }else ++failure;


    //quadratic interpolant:
    res =  CGAL::quadratic_interpolation(coords.begin(), coords.end(),
					 norm,points[i],
					 CGAL::Data_access<Point_value_map>
					 (values),
					 CGAL::Data_access<Point_vector_map>
					 (gradients),
					 Traits());
    if(res.second){
      error = CGAL_NTS abs(res.first - exact_value);
      q_total += error;
      if (error > q_max) q_max = error;
    }else ++failure;

    //Sibson interpolant: version without sqrt:
    res =  CGAL::sibson_c1_interpolation_square(coords.begin(),
     						  coords.end(), norm,
     						  points[i],
     						  CGAL::Data_access<Point_value_map>
     						  (values),
     						  CGAL::Data_access<Point_vector_map>
     						  (gradients),
						Traits());
     //error statistics
    if(res.second){
      error = CGAL_NTS abs(res.first - exact_value);
      ssquare_total += error;
      if (error > ssquare_max) ssquare_max = error;
    } else ++failure;

    //with sqrt(the traditional):
    res =  CGAL::sibson_c1_interpolation(coords.begin(),
					 coords.end(), norm,
					 points[i],
					 CGAL::Data_access<Point_value_map>
					 (values),
					 CGAL::Data_access<Point_vector_map>
					 (gradients),
					   Traits());

    //error statistics
    if(res.second){
      error = CGAL_NTS abs(res.first - exact_value);
      s_total += error;
      if (error > s_max) s_max = error;
    } else ++failure;
  }

  /************** end of Interpolation: dump statistics **************/
  std::cout << "Result: -----------------------------------" << std::endl;
  std::cout <<  "Interpolation of '" << alpha <<" + "
	    << beta1<<" x + "
	    << beta2 << " y + " << gamma1 <<" x^2 + "  << gamma2+ gamma3
	    <<" xy + "  << gamma4 << " y^2'" << std::endl;
  std::cout << "Knowing " << m << " sample points. Interpolation on "
	    << n <<" random points. "<< std::endl;
  std::cout <<"Average function value "
	    << (1.0/n) * CGAL::to_double(total_value)
	    << ", nb of failures "<< failure << std::endl;

  std::cout << "linear interpolant mean error  "
	    << CGAL::to_double(l_total)/n << "  max "
	    << CGAL::to_double(l_max) <<std::endl;
  std::cout << "quadratic interpolant  mean error  "
	    << CGAL::to_double(q_total)/n << "  max "
	    << CGAL::to_double(q_max) << std::endl;
  std::cout << "Farin interpolant  mean error  "
	    << CGAL::to_double(f_total)/n << "  max "
	    << CGAL::to_double(f_max)  << std::endl;
  std::cout << "Sibson interpolant(classic) mean error  "
	    << CGAL::to_double(s_total)/n << "  max "
	    << CGAL::to_double(s_max)  << std::endl;
  std::cout << "Sibson interpolant(square_dist) mean error  "
	    << CGAL::to_double(ssquare_total)/n << "  max "
	    << CGAL::to_double(ssquare_max)  << std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
