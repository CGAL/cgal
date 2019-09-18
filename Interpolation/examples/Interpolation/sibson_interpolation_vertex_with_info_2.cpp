#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>
#include <CGAL/interpolation_functions.h>

#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Interpolation_gradient_fitting_traits_2<K>    Traits;

typedef K::FT                                               Coord_type;
typedef K::Point_2                                          Point;
typedef K::Vector_2                                         Vector;

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
typedef CGAL::Delaunay_triangulation_2<K,Tds>               Delaunay_triangulation;
typedef Delaunay_triangulation::Vertex_handle               Vertex_handle;

template <typename V, typename T>
struct Value_function
{
  typedef V                                                 argument_type;
  typedef std::pair<T, bool>                                result_type;

  // read operation
  result_type operator()(const argument_type& a) const {
    return result_type(a->info().value, true);
  }
};

template <typename V, typename G>
struct Gradient_function
  : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
  typedef V                                                argument_type;
  typedef std::pair<G, bool>                               result_type;

  // read operation
  result_type operator()(const argument_type& a) const {
    return std::make_pair(a->info().gradient, a->info().gradient != CGAL::NULL_VECTOR);
  }

  // write operation
  const Gradient_function& operator=(const std::pair<V, G>& p) const {
    p.first->info().gradient = p.second;
    return *this;
  }

  const Gradient_function& operator++(int) const { return *this; }
  const Gradient_function& operator*() const { return *this; }
};

int main()
{
  Delaunay_triangulation dt;

  // Note that a supported alternative to creating the functors below is to use lambdas
  Value_function<Vertex_handle, Coord_type> value_function;
  Gradient_function<Vertex_handle, Vector> gradient_function;

  // parameters for spherical function:
  Coord_type a(0.25), bx(1.3), by(-0.7), c(0.2);
  for(int y=0 ; y<4 ; ++y) {
    for(int x=0 ; x<4 ; ++x) {
      K::Point_2 p(x,y);
      Vertex_handle vh = dt.insert(p);
      Coord_type value = a + bx* x+ by*y + c*(x*x+y*y);
      vh->info().value = value; // store the value directly in the vertex
    }
  }

  CGAL::sibson_gradient_fitting_nn_2(dt,
                                     gradient_function,
                                     CGAL::Identity<std::pair<Vertex_handle, Vector> >(),
                                     value_function,
                                     Traits());

  // coordinate computation
  K::Point_2 p(1.6, 1.4);
  std::vector<std::pair<Vertex_handle, Coord_type> > coords;
  typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> > Identity;
  Coord_type norm = CGAL::natural_neighbor_coordinates_2(dt,
                                                         p,
                                                         std::back_inserter(coords),
                                                         Identity()).second;

  // Sibson interpolant: version without sqrt:
  std::pair<Coord_type, bool> res = CGAL::sibson_c1_interpolation_square(coords.begin(),
                                                                         coords.end(),
                                                                         norm,
                                                                         p,
                                                                         value_function,
                                                                         gradient_function,
                                                                         Traits());

  if(res.second)
    std::cout << "Tested interpolation on " << p
              << " interpolation: " << res.first << " exact: "
              << a + bx*p.x() + by*p.y() + c*(p.x()*p.x() + p.y()*p.y())
              << std::endl;
  else
    std::cout << "C^1 Interpolation not successful." << std::endl
              << " not all function_gradients are provided." << std::endl
              << " You may resort to linear interpolation." << std::endl;

  return EXIT_SUCCESS;
}
