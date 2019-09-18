#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Origin.h>
#include <CGAL/tuple.h>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               Coord_type;
typedef K::Vector_2                                         Vector;
typedef K::Point_2                                          Point;

template <typename V, typename G>
struct Value_and_gradient
{
  Value_and_gradient() : value(), gradient(CGAL::NULL_VECTOR) {}

  V value;
  G gradient;
};

typedef CGAL::Triangulation_vertex_base_with_info_2<
                Value_and_gradient<Coord_type, Vector>, K>  Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>            Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>              Delaunay_triangulation;
typedef Delaunay_triangulation::Vertex_handle               Vertex_handle;
typedef CGAL::Interpolation_traits_2<K>                     Traits;

typedef std::vector<std::pair<Vertex_handle, Coord_type> >  Coordinate_vector;

template <typename V, typename T>
struct Value_function
{
  typedef V                                                 argument_type;
  typedef std::pair<T, bool>                                result_type;

  result_type operator()(const argument_type& a) const {
    return result_type(a->info().value, true);
  }
};

template <typename V, typename G>
struct Gradient_function
  : public std::iterator<std::output_iterator_tag, void, void, void, void>
{

  typedef V                                                 argument_type;
  typedef std::pair<G, bool>                                result_type;

  result_type operator()(const argument_type& a) const {
    return std::make_pair(a->info().gradient, a->info().gradient != CGAL::NULL_VECTOR);
  }

  const Gradient_function& operator=(const std::pair<V, G>& p) const {
    p.first->info().gradient = p.second;
    return *this;
  }

  const Gradient_function& operator++(int) const { return *this; }
  const Gradient_function& operator*() const { return *this; }
};

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
  CGAL::Random rng(1513114263);
  CGAL::Random_points_in_disc_2<Point> g(r_d, rng);
  std::copy_n( g, n+m, std::back_inserter(points));

  Delaunay_triangulation T;

  // Note that a supported alternative to creating the functors below is to use lambdas
  Value_function<Vertex_handle, Coord_type> value_function;
  Gradient_function<Vertex_handle, Vector> gradient_function;

  //parameters for quadratic function:
  Coord_type alpha = Coord_type(1.0),
             beta1 = Coord_type(2.0),
             beta2 = Coord_type(1.0),
             gamma1 = Coord_type(0.3),
             gamma2 = Coord_type(0.0),
             gamma3 = Coord_type(0.0),
             gamma4 = Coord_type(0.3);

  for(int j=0; j<n ; j++)
  {
    Vertex_handle vh = T.insert(points[j]);

    //determine function value/gradient:
    Coord_type x(points[j].x());
    Coord_type y(points[j].y());

    Coord_type value = alpha + beta1*x + beta2*y + gamma1*(x*x) +
                       gamma4*(y*y) + (gamma2+ gamma3) *(x*y);
    Vector gradient(beta1+ (gamma2 + gamma3)*y + Coord_type(2)*(gamma1*x),
                    beta2+ (gamma2 + gamma3)*x + Coord_type(2)*(gamma4*y));
    vh->info().value = value;
    vh->info().gradient = gradient;
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
  for(int i=n; i<n+m; i++)
  {
    Coord_type x(points[i].x());
    Coord_type y(points[i].y());

    Coord_type exact_value = alpha + beta1*x + beta2*y + gamma1*(x*x) +
                             gamma4*(y*y) + (gamma2+ gamma3) *(x*y);

    total_value += exact_value;

    //Coordinate_vector:
    Coordinate_vector coords;
    typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> > Identity;
    Coord_type norm = CGAL::natural_neighbor_coordinates_2(T,
                                                           points[i],
                                                           std::back_inserter(coords),
                                                           Identity()).second;
    assert(norm > 0);

    //linear interpolant:
    l_value = CGAL::linear_interpolation(coords.begin(), coords.end(),
                                         norm, value_function);

    error = CGAL_NTS abs(l_value - exact_value);
    l_total += error;
    if (error > l_max) l_max = error;

    //Farin interpolant:
    res = CGAL::farin_c1_interpolation(coords.begin(),
                                       coords.end(), norm,points[i],
                                       value_function,
                                       gradient_function,
                                       Traits());


    if(res.second)
    {
      error = CGAL_NTS abs(res.first - exact_value);
      f_total += error;
      if (error > f_max)
        f_max = error;
    }
    else
    {
      ++failure;
    }

    //quadratic interpolant:
    res = CGAL::quadratic_interpolation(coords.begin(), coords.end(),
                                        norm, points[i],
                                        value_function,
                                        gradient_function,
                                        Traits());

    if(res.second)
    {
      error = CGAL_NTS abs(res.first - exact_value);
      q_total += error;
      if (error > q_max)
        q_max = error;
    }
    else
    {
      ++failure;
    }

    //Sibson interpolant: version without sqrt:
    res = CGAL::sibson_c1_interpolation_square(coords.begin(),
                                               coords.end(), norm,
                                               points[i],
                                               value_function,
                                               gradient_function,
                                               Traits());
    //error statistics
    if(res.second)
    {
      error = CGAL_NTS abs(res.first - exact_value);
      ssquare_total += error;
      if (error > ssquare_max)
        ssquare_max = error;
    }
    else
    {
      ++failure;
    }

    //with sqrt(the traditional):
    res = CGAL::sibson_c1_interpolation(coords.begin(),
                                        coords.end(), norm,
                                        points[i],
                                        value_function,
                                        gradient_function,
                                        Traits());

    //error statistics
    if(res.second)
    {
      error = CGAL_NTS abs(res.first - exact_value);
      s_total += error;
      if (error > s_max)
        s_max = error;
    }
    else
    {
      ++failure;
    }
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

  return EXIT_SUCCESS;
}
