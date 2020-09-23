// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto
//                 Mael Rouxel-Labbé

#include <CGAL/Interpolation/internal/helpers.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>

#include <CGAL/algorithm.h>
#include <CGAL/double.h>
#include <CGAL/function_objects.h>
#include <functional>
#include <CGAL/Origin.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>

#include <boost/utility/result_of.hpp>

#include <iostream>
#include <cassert>
#include <utility>

template < typename Traits_ >
struct Extract_point
{
  typedef Traits_                               Traits;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Weighted_point_2     Weighted_point_2;

  typedef typename Traits::Construct_point_2    Construct_point_2;

  Extract_point(const Traits& traits = Traits()) : traits(traits) {}

  const Point_2& operator()(const Point_2& p) const { return p; }

  Point_2 operator()(const Weighted_point_2& wp) const {
    return traits.construct_point_2_object()(wp);
  }

  template <typename VH>
  typename boost::result_of<const Construct_point_2(const Point_2&)>::type
  operator()(const VH& vh) const {
    return traits.construct_point_2_object()(vh->point());
  }

private:
  Traits traits;
};

template <typename V, typename T>
struct Value_function
{
  typedef V                                                 argument_type;
  typedef std::pair<T, bool>                                result_type;

  Value_function(std::size_t i) : index(i) { }

  result_type operator()(const argument_type& a) const {
    return result_type(a->info()[index].value, true);
  }

private:
  std::size_t index;
};

template <typename V, typename G>
struct Gradient_function
  : public CGAL::cpp98::iterator<std::output_iterator_tag, void, void, void, void>
{

  typedef V                                                 argument_type;
  typedef std::pair<G, bool>                                result_type;

  Gradient_function(std::size_t i) : index(i) { }

  result_type operator()(const argument_type& a) const {
    return std::make_pair(a->info()[index].gradient,
                          a->info()[index].gradient != CGAL::NULL_VECTOR);
  }

private:
  std::size_t index;
};


template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
               typename std::iterator_traits<ForwardIterator>::value_type::second_type norm)
{
  typename std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
  for(; first!=beyond; ++first)
    sum += first->second;

  if(norm != sum)
  {
    std::cerr << "norm / sum differs: " << norm << " " << sum << std::endl;
    return false;
  }
  else
  {
    return true;
  }
}

template < class Tr, class ForwardIterator, class Point >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
                     typename std::iterator_traits<ForwardIterator>::value_type::second_type norm,
                     const Point& p,
                     const typename std::iterator_traits<ForwardIterator>::value_type::second_type& tolerance)
{
  typedef typename Tr::Geom_traits                                               Gt;
  typedef typename Gt::Point_2                                                   Bare_point;
  Extract_point<Gt> cp;

  Bare_point b(CGAL::ORIGIN);
  for(; first!=beyond; ++first)
    b = b + (first->second / norm) * (cp(first->first) - CGAL::ORIGIN);

  return CGAL::squared_distance(cp(p), b) <= tolerance;
}

/////////////////////////////////////////////////////////////////
// Accessory function testing functions that require sqrt().
// Doesn't instantiate anything if RT doesn't support sqrt().
template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Gt, class Point>
bool
_test_sibson_c1_interpolation_sqrt(ForwardIterator , ForwardIterator ,
                                   const typename std::iterator_traits<ForwardIterator>::value_type::second_type& /* norm */,
                                   const Point& /* p */,
                                   ValueFunctor /* f */,
                                   GradFunctor /* grad_f */,
                                   const Gt& /*geom_traits */,
                                   const typename std::iterator_traits<ForwardIterator>::value_type::second_type& /* tolerance */,
                                   const typename std::iterator_traits<ForwardIterator>::value_type::second_type& /* exact_value */,
                                   CGAL::Integral_domain_without_division_tag t)
{
  std::cout << std::endl
            << "NOTE: "
            << typeid(t).name()
            << " doesn't support sqrt() and sibson_c1_interpolation is thus not tested." << std::endl;
  return true;
}

template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Gt, class Point>
bool _test_sibson_c1_interpolation_sqrt(ForwardIterator first, ForwardIterator beyond,
                                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                                        const Point& p,
                                        ValueFunctor f,
                                        GradFunctor grad_f,
                                        const Gt& geom_traits,
                                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& tolerance,
                                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& exact_value,
                                        CGAL::Field_with_sqrt_tag)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  value_functor_result_type res = CGAL::sibson_c1_interpolation(first, beyond,
                                                                norm, p, f,
                                                                grad_f, geom_traits);

  return res.second && (CGAL_NTS abs(res.first-exact_value) <= tolerance);
}

template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Gt, class Point>
bool test_interpolation_with_value(ForwardIterator first, ForwardIterator beyond,
                                   const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                                   const Point& p,
                                   const typename boost::result_of<
                                                    ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
                                                    ::type::first_type exact_value,
                                   ValueFunctor f,
                                   GradFunctor grad_f,
                                   const Gt& geom_traits,
                                   const int& i,
                                   const typename std::iterator_traits<ForwardIterator>::value_type::second_type& tolerance)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;
  typedef typename value_functor_result_type::first_type                          Value_type;

  if(i == 0)
  {
    Value_type val = CGAL::linear_interpolation(first, beyond, norm, f);
    assert(CGAL_NTS abs(val - exact_value) <= tolerance);
  }

  value_functor_result_type res = CGAL::quadratic_interpolation(first, beyond, norm, p, f,
                                                                grad_f, geom_traits);
  assert(res.second && (CGAL_NTS abs(res.first - exact_value) <= tolerance));

  if(i<2)
  {
    //without sqrt:
    res = CGAL::sibson_c1_interpolation_square(first, beyond,
                                               norm, p, f,
                                               grad_f, geom_traits);
    assert(res.second && (CGAL_NTS abs(res.first - exact_value) <= tolerance));

    //with sqrt:
    typedef CGAL::Algebraic_structure_traits<Value_type> AST;

    assert(_test_sibson_c1_interpolation_sqrt(first, beyond, norm, p, f, grad_f,
                                              geom_traits, tolerance, exact_value,
                                              typename AST::Algebraic_category()));
  }

  res = CGAL::farin_c1_interpolation(first, beyond, norm, p, f, grad_f, geom_traits);
  assert(res.second);
  assert(CGAL_NTS abs(res.first - exact_value) <= tolerance);

  return true;
}

template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Gt, class Point, class Value_type>
bool test_interpolation(ForwardIterator first, ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const Point& p,
                        const Value_type exact_value,
                        ValueFunctor f,
                        GradFunctor grad_f,
                        const Gt& geom_traits,
                        const int& i,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& tolerance)
{
  return test_interpolation_with_value(first, beyond, norm, p, exact_value, f, grad_f, geom_traits, i, tolerance);
}

template <class Dt>
void _test_interpolation_functions_2_Delaunay_without_OutputFunctor(const Dt&, const typename Dt::Geom_traits::FT& tolerance)
{
  std::cout << "Testing backward compatibility..." << std::endl;

  CGAL::Set_ieee_double_precision pfr;
  Dt T;

  int n=20, m=20;
  double r = 3;
  double max_value = 5;

  typedef typename Dt::Geom_traits                             Gt;
  typedef CGAL::Interpolation_traits_2<Gt>                     Traits;

  typedef typename Dt::Point                                   Point;

  typedef typename Gt::FT                                      Coord_type;
  typedef typename Gt::Vector_2                                Vector;

  typedef std::map<Point, Coord_type, typename Gt::Less_xy_2>  Point_value_map;
  typedef std::map<Point, Vector, typename Gt::Less_xy_2>      Point_vector_map;

  typedef std::vector<std::pair<Point, Coord_type> >           Point_coordinate_vector;

  std::cout << "NN2: Testing random points." << std::endl;

  //test random points in a square of length r:
  std::vector<Point> points;
  points.reserve(n+m);

  //put four bounding box points:
  points.push_back(Point(-r, -r));
  points.push_back(Point(r, -r));
  points.push_back(Point(-r, r));
  points.push_back(Point(r, r));

  // Create n+m-4 points within a disc of radius 2
  CGAL::Random_points_in_square_2<Point> g(r);
  std::copy_n(g, n+m, std::back_inserter(points));

  CGAL::Random random;

  Point_value_map values[3];
  Point_vector_map gradients[3];

  Coord_type alpha = Coord_type(random.get_double(-max_value, max_value)),
             beta1 = Coord_type(random.get_double(-max_value, max_value)),
             beta2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma1 = Coord_type(random.get_double(-max_value, max_value)),
             gamma2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma3 = Coord_type(random.get_double(-max_value, max_value));

  //INSERTION + DET. of GRADIENT for n DATA POINTS :
  for(int j=0; j<n; ++j)
  {
    T.insert(points[j]);

    gradients[0].insert(std::make_pair(points[j], Vector(beta1, beta2)));

    gradients[1].insert(std::make_pair(points[j],
                                       Vector(beta1 + Coord_type(2)*gamma1*points[j].x(),
                                              beta2 + Coord_type(2)*gamma1*points[j].y())));
    gradients[2].insert(std::make_pair(points[j],
                                       Vector(beta1 + Coord_type(2)*gamma1*points[j].x() + gamma3*points[j].y(),
                                              beta2 + Coord_type(2)*gamma2*points[j].y() + gamma3*points[j].x())));
  }

  //DETERMINE VALUES FOR n DATA POINTS AND m RANDOM TEST POINTS:
  for(int j=0; j<n+m; j++)
  {
    //linear function
    values[0].insert(std::make_pair(points[j], alpha + beta1*points[j].x() + beta2*points[j].y()));

    //spherical function:
    values[1].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                       beta2*points[j].y() +
                                                       gamma1*points[j].x()*points[j].x()+
                                                       gamma1*points[j].y()*points[j].y()));

    //quadratic function
    values[2].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                       beta2*points[j].y() +
                                                       gamma1*points[j].x()*points[j].x() +
                                                       gamma2*points[j].y()*points[j].y() +
                                                       gamma3*points[j].x()*points[j].y()));
  }

  //INTERPOLATION OF RANDOM POINTS:
  Coord_type norm;
  Point_coordinate_vector coords;
  for(int j=n; j<n+m; ++j)
  {
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
        CGAL::natural_neighbor_coordinates_2(T, points[j], std::back_inserter(coords));
    assert(coordinate_result.third);
    norm = coordinate_result.second;

    bool is_equal = test_norm(coords.begin(), coords.end(), norm);
    assert(norm > 0);
    assert(is_equal);
    is_equal = test_barycenter<Dt>(coords.begin(), coords.end(), norm, points[j], tolerance);
    assert(is_equal);

    for(int i=0; i<3; ++i)
    {
      assert(test_interpolation(coords.begin(), coords.end(), norm,
                                points[j], values[i][points[j]],
                                CGAL::Data_access< Point_value_map >(values[i]),
                                CGAL::Data_access< Point_vector_map >(gradients[i]),
                                Traits(), i, tolerance));
    }
    coords.clear();
  }

  //TESTING THE GRADIENT APPRXIMATION METHOD:
    // sibson_gradient_fitting is tested when calling nn_2
  std::cout << "Testing gradient estimation method on random points." << std::endl;

  typedef CGAL::Interpolation_gradient_fitting_traits_2<Gt> GradTraits;
  Point_vector_map approx_gradients[2];
  CGAL::sibson_gradient_fitting_nn_2(T,
                                     std::inserter(approx_gradients[0], approx_gradients[0].begin()),
                                     CGAL::Data_access<Point_value_map>(values[0]), GradTraits());

  CGAL::sibson_gradient_fitting_nn_2(T,
                                     std::inserter(approx_gradients[1], approx_gradients[1].begin()),
                                     CGAL::Data_access<Point_value_map>(values[1]), GradTraits());

  for(int j=0; j<n; ++j)
  {
    std::pair<Vector, bool> res = CGAL::Data_access<Point_vector_map>(approx_gradients[0])(points[j]);

    if(res.second)
    {
      // if it is the exact computation kernel: test the equality:
      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Point_vector_map>(gradients[0])(points[j]).first);
      res = CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]);

      // if one exists->the other must also exist
      assert(res.second);

      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Point_vector_map>(gradients[1])(points[j]).first);
    }
    else
    {
      assert(!CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]).second);
    }
  }

  //TESTING A POINT == A DATA POINT:
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
      CGAL::natural_neighbor_coordinates_2(T, points[n/2], std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));

  typename std::vector<std::pair<Point, Coord_type> >::iterator ci = coords.begin();
  assert(ci->first == points[n/2]);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci==coords.end());

  for(int j=0; j<3; ++j)
  {
    assert(test_interpolation(coords.begin(), coords.end(), norm,
                              points[n/2], values[j][points[n/2]],
                              CGAL::Data_access<Point_value_map>(values[j]),
                              CGAL::Data_access<Point_vector_map>(gradients[j]),
                              Traits(), j, tolerance));
  }
  coords.clear();
}

template <class Dt>
void _test_interpolation_functions_2_Delaunay_with_OutputFunctor(const Dt&, const typename Dt::Geom_traits::FT& tolerance)
{
  std::cout << "Testing with OutputFunctor..." << std::endl;

  CGAL::Set_ieee_double_precision pfr;
  Dt T;

  int n=20, m=20;
  double r = 3;
  double max_value = 5;

  typedef typename Dt::Geom_traits                               Gt;
  typedef CGAL::Interpolation_traits_2<Gt>                       Traits;

  typedef typename Dt::Vertex_handle                             Vertex_handle;

  typedef typename Gt::FT                                        Coord_type;
  typedef typename Dt::Point                                     Point;
  typedef typename Gt::Vector_2                                  Vector;

  typedef std::vector<std::pair<Vertex_handle, Coord_type> >     Coordinate_vector;
  typedef typename Coordinate_vector::const_iterator             CV_cit;
  typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> >  Output_functor;

  typedef std::map<Point, Coord_type>                            Point_value_map;
  typedef std::map<Point, Vector>                                Point_vector_map;

  std::cout << "NN2: Testing random points." << std::endl;

  //test random points in a square of length r:
  std::vector<Point> points;
  points.reserve(n+m);

  //put four bounding box points:
  points.push_back(Point(-r, -r));
  points.push_back(Point(r, -r));
  points.push_back(Point(-r, r));
  points.push_back(Point(r, r));

  // Create n+m-4 points within a disc of radius 2
  CGAL::Random_points_in_square_2<Point> g(r);
  std::copy_n(g, n+m, std::back_inserter(points));

  CGAL::Random random;

  Coord_type alpha = Coord_type(random.get_double(-max_value, max_value)),
             beta1 = Coord_type(random.get_double(-max_value, max_value)),
             beta2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma1 = Coord_type(random.get_double(-max_value, max_value)),
             gamma2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma3 = Coord_type(random.get_double(-max_value, max_value));

  //INSERTION + DET. of GRADIENT for n DATA POINTS :
  //DETERMINE VALUES FOR n DATA POINTS AND m RANDOM TEST POINTS:
  Point_value_map exact_values[3];

  std::map<Point, Vertex_handle> p_to_vh;

  for(int j=0; j<n+m; ++j)
  {
    Vector gradient0(beta1, beta2);
    Vector gradient1(beta1 + Coord_type(2)*gamma1*points[j].x(),
                     beta2 + Coord_type(2)*gamma1*points[j].y());
    Vector gradient2(beta1 + Coord_type(2)*gamma1*points[j].x() + gamma3*points[j].y(),
                     beta2 + Coord_type(2)*gamma2*points[j].y() + gamma3*points[j].x());

    Coord_type value0 = alpha + beta1*points[j].x() + beta2*points[j].y();
    Coord_type value1 = alpha + beta1*points[j].x()
                              + beta2*points[j].y()
                              + gamma1*points[j].x()*points[j].x()
                              + gamma1*points[j].y()*points[j].y();
    Coord_type value2 = alpha + beta1*points[j].x()
                              + beta2*points[j].y()
                              + gamma1*points[j].x()*points[j].x()
                              + gamma2*points[j].y()*points[j].y()
                              + gamma3*points[j].x()*points[j].y();

    if(j<n) // only insert n points
    {
      Vertex_handle vh = T.insert(points[j]);
      p_to_vh[points[j]] = vh;

      vh->info()[0].gradient = gradient0;
      vh->info()[1].gradient = gradient1;
      vh->info()[2].gradient = gradient2;


      vh->info()[0].value = value0;
      vh->info()[1].value = value1;
      vh->info()[2].value = value2;
    }
    else
    {
      exact_values[0][points[j]] = value0;
      exact_values[1][points[j]] = value1;
      exact_values[2][points[j]] = value2;
    }
  }

  //INTERPOLATION OF RANDOM POINTS:
  Coord_type norm;

  Coordinate_vector coords;
  Output_functor out_fct;

  for(int j=n; j<n+m; ++j)
  {
    CGAL::Triple<std::back_insert_iterator<Coordinate_vector>, Coord_type, bool> coordinate_result =
        CGAL::natural_neighbor_coordinates_2(T, points[j], std::back_inserter(coords), out_fct);
    assert(coordinate_result.third);
    norm = coordinate_result.second;

    bool is_equal = test_norm(coords.begin(), coords.end(), norm);
    assert(norm > 0);
    assert(is_equal);
    is_equal = test_barycenter<Dt>(coords.begin(), coords.end(), norm, points[j], tolerance);
    assert(is_equal);

    assert(test_interpolation(coords.begin(), coords.end(), norm,
                              points[j], exact_values[0][points[j]],
                              [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[0].value, true); },
                              [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[0].gradient, true); },
                              Traits(), 0, tolerance));


    // wrapping the lambda in a std function
    std::function<std::pair<Coord_type, bool>(const Vertex_handle)> value_function_1 =
      [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[1].value, true); };

    std::function<std::pair<Vector, bool>(const Vertex_handle)> gradient_function_1 =
      [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[1].gradient, true); };

    assert(test_interpolation(coords.begin(), coords.end(), norm,
                              points[j], exact_values[1][points[j]],
                              value_function_1, gradient_function_1,
                              Traits(), 1, tolerance));

    assert(test_interpolation(coords.begin(), coords.end(), norm,
                              points[j], exact_values[2][points[j]],
                              [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[2].value, true); },
                              [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[2].gradient, true); },
                              Traits(), 2, tolerance));

    coords.clear();
  }

  //TESTING THE GRADIENT APPRXIMATION METHOD:
    // sibson_gradient_fitting is tested when calling nn_2
  std::cout << "Testing gradient estimation method on random points." << std::endl;

  typedef CGAL::Interpolation_gradient_fitting_traits_2<Gt> GradTraits;

  Point_vector_map approx_gradients[2];

  {
    CGAL::sibson_gradient_fitting_nn_2(T,
                                       std::inserter(approx_gradients[0], approx_gradients[0].begin()), // OutputIterator
                                       CGAL::Interpolation::internal::Extract_point_in_pair<Dt, Vector>(), // OutputFunctor
                                       [](const Vertex_handle vh)
                                          -> std::pair<Coord_type, bool>
                                         { return std::make_pair(vh->info()[0].value, true); },
                                       GradTraits());

    std::function<std::pair<Coord_type, bool>(const Vertex_handle)> value_function_1 =
      [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[1].value, true); };

    CGAL::sibson_gradient_fitting_nn_2(T,
                                       std::inserter(approx_gradients[1], approx_gradients[1].begin()),
                                       CGAL::Interpolation::internal::Extract_point_in_pair<Dt, Vector>(),
                                       value_function_1,
                                       GradTraits());
  }

  for(int j=0; j<n; ++j)
  {
    std::pair<Vector, bool> res = CGAL::Data_access<Point_vector_map>(approx_gradients[0])(points[j]);

    if(res.second)
    {
      Gradient_function<Vertex_handle, Vector> gradient_function_0(0);
      Gradient_function<Vertex_handle, Vector> gradient_function_1(1);

      // if it is the exact computation kernel: test the equality:
      assert(tolerance > Coord_type(0) ||
             res.first == (gradient_function_0(p_to_vh[points[j]])).first);
      res = CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]);

      // if one exists->the other must also exist
      assert(res.second);

      assert(tolerance > Coord_type(0) ||
             res.first == gradient_function_1(p_to_vh[points[j]]).first);
    }
    else
    {
      assert(!CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]).second);
    }
  }

  //TESTING A POINT == A DATA POINT:
  CGAL::Triple<std::back_insert_iterator<Coordinate_vector>, Coord_type, bool> coordinate_result =
      CGAL::natural_neighbor_coordinates_2(T, points[n/2], std::back_inserter(coords), out_fct);
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));

  CV_cit ci = coords.begin();
  assert(ci->first == p_to_vh[points[n/2]]);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci == coords.end());

  Value_function<Vertex_handle, Coord_type> value_function_0(0);
  Value_function<Vertex_handle, Coord_type> value_function_2(2);

  assert(test_interpolation(coords.begin(), coords.end(), norm,
                            points[n/2], value_function_0(p_to_vh[points[n/2]]).first,
                            [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[0].value, true); },
                            [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[0].gradient, true); },
                            Traits(), 0, tolerance));

  // wrapping the lambda in a std function
  std::function<std::pair<Coord_type, bool>(const Vertex_handle)> value_function_1 =
      [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[1].value, true); };

  std::function<std::pair<Vector, bool>(const Vertex_handle)> gradient_function_1 =
      [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[1].gradient, true); };

  assert(test_interpolation(coords.begin(), coords.end(), norm,
                            points[n/2], value_function_1(p_to_vh[points[n/2]]).first,
                            value_function_1, gradient_function_1,
                            Traits(), 1, tolerance));

  assert(test_interpolation(coords.begin(), coords.end(), norm,
                            points[n/2], value_function_2(p_to_vh[points[n/2]]).first,
                            [](const Vertex_handle vh) -> std::pair<Coord_type, bool> { return std::make_pair(vh->info()[2].value, true); },
                            [](const Vertex_handle vh) -> std::pair<Vector, bool> { return std::make_pair(vh->info()[2].gradient, true); },
                            Traits(), 2, tolerance));

  coords.clear();
}

template <class Rt>
void _test_interpolation_functions_2_regular_without_OutputFunctor(const Rt&, const typename Rt::Geom_traits::FT& tolerance)
{
  std::cout << "Testing backward compatibility..." << std::endl;

  CGAL::Set_ieee_double_precision pfr;
  Rt T;

  std::size_t n=20, m=20;
  double r = 3.;
  double max_value = 5.;

  typedef typename Rt::Geom_traits                             Gt;
  typedef CGAL::Interpolation_traits_2<Gt>                     Traits;

  typedef typename Rt::Bare_point                              Bare_point;
  typedef typename Rt::Weighted_point                          Weighted_point;

  typedef typename Gt::FT                                      Coord_type;
  typedef typename Gt::Vector_2                                Vector;

  typedef std::map<Weighted_point, Coord_type>                 Point_value_map;
  typedef std::map<Weighted_point, Vector>                     Point_vector_map;

  typedef std::vector<std::pair<Weighted_point, Coord_type> >  Point_coordinate_vector;

  std::cout << "RN2: Testing random points." << std::endl;

  // test random points in a square of length r:
  std::vector<Weighted_point> points;
  points.reserve(n+m);

  // put four bounding box points to never have infinite vertices
  // move them away from the cube of side 'r' to be sure they are not hidden
  double rr = 10 * r;
  points.push_back(Weighted_point(Bare_point(-rr, -rr), 0));
  points.push_back(Weighted_point(Bare_point(rr, -rr), 0));
  points.push_back(Weighted_point(Bare_point(-rr, rr), 0));
  points.push_back(Weighted_point(Bare_point(rr, rr), 0));

  // Create n+m-4 points within a disc of radius 2, with random weights
  CGAL::Random random;
  CGAL::Random_points_in_square_2<Bare_point> g(r);

  while(points.size() != n+m)
  {
    Weighted_point p(*g++, random.get_double(0., 0.5));
    points.push_back(p);
  }

  Point_value_map values[3];
  Point_vector_map gradients[3];

  Coord_type alpha = Coord_type(random.get_double(-max_value, max_value)),
             beta1 = Coord_type(random.get_double(-max_value, max_value)),
             beta2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma1 = Coord_type(random.get_double(-max_value, max_value)),
             gamma2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma3 = Coord_type(random.get_double(-max_value, max_value));

  //INSERTION + DET. of GRADIENT for n DATA POINTS :
  for(std::size_t j=0; j<n; ++j)
  {
    T.insert(points[j]);

    gradients[0].insert(std::make_pair(points[j], Vector(beta1, beta2)));

    gradients[1].insert(std::make_pair(points[j],
                                       Vector(beta1 + Coord_type(2)*gamma1*points[j].x(),
                                              beta2 + Coord_type(2)*gamma1*points[j].y())));
    gradients[2].insert(std::make_pair(points[j],
                                       Vector(beta1 + Coord_type(2)*gamma1*points[j].x() + gamma3*points[j].y(),
                                              beta2 + Coord_type(2)*gamma2*points[j].y() + gamma3*points[j].x())));
  }

  //DETERMINE VALUES FOR n DATA POINTS AND m RANDOM TEST POINTS:
  for(std::size_t j=0; j<n+m; j++)
  {
    //linear function
    values[0].insert(std::make_pair(points[j], alpha + beta1*points[j].x() + beta2*points[j].y()));

    //spherical function:
    values[1].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                       beta2*points[j].y() +
                                                       gamma1*points[j].x()*points[j].x()+
                                                       gamma1*points[j].y()*points[j].y()));

    //quadratic function
    values[2].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                       beta2*points[j].y() +
                                                       gamma1*points[j].x()*points[j].x() +
                                                       gamma2*points[j].y()*points[j].y() +
                                                       gamma3*points[j].x()*points[j].y()));
  }

  //INTERPOLATION OF RANDOM POINTS:
  Coord_type norm;
  Point_coordinate_vector coords;
  for(std::size_t j=n; j<n+m; ++j)
  {
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
        CGAL::regular_neighbor_coordinates_2(T, points[j], std::back_inserter(coords));
    assert(coordinate_result.third);

    norm = coordinate_result.second;
    assert(norm >= 0); // can be null if the weight was too small
    bool is_equal = test_norm(coords.begin(), coords.end(), norm);
    assert(is_equal);

    // further tests are pointless if there are no coordinates
    if(norm == 0)
      continue;

    is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, points[j], tolerance);
    assert(is_equal);

    for(int i=0; i<3; ++i)
    {
      assert(test_interpolation(coords.begin(), coords.end(), norm,
                                points[j], values[i][points[j]],
                                CGAL::Data_access< Point_value_map >(values[i]),
                                CGAL::Data_access< Point_vector_map >(gradients[i]),
                                Traits(), i, tolerance));
    }
    coords.clear();
  }

  //TESTING THE GRADIENT APPRXIMATION METHOD:
    // sibson_gradient_fitting is tested when calling rn_2
  std::cout << "Testing gradient estimation method on random points." << std::endl;
  typedef CGAL::Interpolation_gradient_fitting_traits_2<Gt> GradTraits;
  Point_vector_map approx_gradients[2];
  CGAL::sibson_gradient_fitting_rn_2(T,
                                     std::inserter(approx_gradients[0], approx_gradients[0].begin()),
                                     CGAL::Data_access<Point_value_map>(values[0]), GradTraits());

  CGAL::sibson_gradient_fitting_rn_2(T,
                                     std::inserter(approx_gradients[1], approx_gradients[1].begin()),
                                     CGAL::Data_access<Point_value_map>(values[1]), GradTraits());

  for(std::size_t j=0; j<n; ++j)
  {
    std::pair<Vector, bool> res = CGAL::Data_access<Point_vector_map>(approx_gradients[0])(points[j]);

    if(res.second)
    {
      // if it is the exact computation kernel: test the equality:
      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Point_vector_map>(gradients[0])(points[j]).first);
      res = CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]);

      // if one exists->the other must also exist
      assert(res.second);
      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Point_vector_map>(gradients[1])(points[j]).first);
    }
    else
    {
      assert(!CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]).second);
    }
  }

  //TESTING A POINT == A DATA POINT:
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
      CGAL::regular_neighbor_coordinates_2(T, points[n/2], std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  // if the point is hidden, the norm will be 0
  assert(norm == Coord_type(0) || norm == Coord_type(1));

  if(norm == Coord_type(1)) // nothing to do if the point is hidden
  {
    typename std::vector<std::pair<Weighted_point, Coord_type> >::iterator ci = coords.begin();
    assert(ci->first == points[n/2]);
    assert(ci->second == Coord_type(1));
    ci++;
    assert(ci == coords.end());

    for(int j=0; j<3; ++j)
    {
      assert(test_interpolation(coords.begin(), coords.end(), norm,
                                points[n/2], values[j][points[n/2]],
                                CGAL::Data_access<Point_value_map>(values[j]),
                                CGAL::Data_access<Point_vector_map>(gradients[j]),
                                Traits(), j, tolerance));
    }
    coords.clear();
  }
  else
  {
    assert(coords.empty()); // point was hidden
  }
}

template <class Rt>
void _test_interpolation_functions_2_regular_with_OutputFunctor(const Rt&, const typename Rt::Geom_traits::FT& tolerance)
{
  CGAL::Set_ieee_double_precision pfr;
  Rt T;

  std::size_t n=20, m=20;
  double r = 3.;
  double max_value = 5.;

  typedef typename Rt::Geom_traits                               Gt;
  typedef CGAL::Interpolation_traits_2<Gt>                       Traits;

  typedef typename Rt::Bare_point                                Bare_point;
  typedef typename Rt::Weighted_point                            Weighted_point;

  typedef typename Gt::FT                                        FT;
  typedef typename Gt::FT                                        Coord_type;
  typedef typename Gt::Vector_2                                  Vector;

  typedef typename Rt::Vertex_handle                             Vertex_handle;

  // These are the values at points which won't be inserted in the triangulation
  typedef std::map<Weighted_point, Coord_type>                   Point_value_map;

  typedef std::map<Vertex_handle, Coord_type>                    Vertex_value_map;
  typedef std::map<Vertex_handle, Vector>                        Vertex_vector_map;

  typedef std::vector<std::pair<Vertex_handle, Coord_type> >     Vertex_coordinate_vector;
  typedef typename Vertex_coordinate_vector::const_iterator      VCV_cit;
  typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> >  Identity_output_functor;

  Identity_output_functor vh_fct;

  std::cout << "RN2: Testing random points." << std::endl;

  // test random points in a square of length r:
  std::vector<Weighted_point> points;
  points.reserve(n+m);

  // put four bounding box points to never have infinite vertices
  // move them away from the cube of side 'r' to be sure they are not hidden
  double rr = 10 * r;
  points.push_back(Weighted_point(Bare_point(-rr, -rr), 0));
  points.push_back(Weighted_point(Bare_point(rr, -rr), 0));
  points.push_back(Weighted_point(Bare_point(-rr, rr), 0));
  points.push_back(Weighted_point(Bare_point(rr, rr), 0));

  // Create n+m-4 points within a disc of radius 2, with random weights
  CGAL::Random random;
  CGAL::Random_points_in_square_2<Bare_point> g(r);

  while(points.size() != n+m)
  {
    Weighted_point p(*g++, FT(random.get_double(0., 0.5)));
    points.push_back(p);
  }

  Vertex_value_map values[3];
  Vertex_vector_map gradients[3];
  Point_value_map pt_values[3];

  Coord_type alpha = Coord_type(random.get_double(-max_value, max_value)),
             beta1 = Coord_type(random.get_double(-max_value, max_value)),
             beta2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma1 = Coord_type(random.get_double(-max_value, max_value)),
             gamma2 = Coord_type(random.get_double(-max_value, max_value)),
             gamma3 = Coord_type(random.get_double(-max_value, max_value));

  //INSERTION + DET. of GRADIENT for n DATA POINTS :
  for(std::size_t j=0; j<n; ++j)
  {
    Vertex_handle vh = T.insert(points[j]);

    if(vh != Vertex_handle())
    {
      gradients[0].insert(std::make_pair(vh, Vector(beta1, beta2)));

      gradients[1].insert(std::make_pair(vh,
                                         Vector(beta1 + Coord_type(2)*gamma1*points[j].x(),
                                                beta2 + Coord_type(2)*gamma1*points[j].y())));
      gradients[2].insert(std::make_pair(vh,
                                         Vector(beta1 + Coord_type(2)*gamma1*points[j].x() + gamma3*points[j].y(),
                                                beta2 + Coord_type(2)*gamma2*points[j].y() + gamma3*points[j].x())));

      //linear function
      values[0].insert(std::make_pair(vh, alpha + beta1*points[j].x() + beta2*points[j].y()));

      //spherical function:
      values[1].insert(std::make_pair(vh, alpha + beta1*points[j].x() +
                                                  beta2*points[j].y() +
                                                  gamma1*points[j].x()*points[j].x()+
                                                  gamma1*points[j].y()*points[j].y()));

      //quadratic function
      values[2].insert(std::make_pair(vh, alpha + beta1*points[j].x() +
                                                  beta2*points[j].y() +
                                                  gamma1*points[j].x()*points[j].x() +
                                                  gamma2*points[j].y()*points[j].y() +
                                                  gamma3*points[j].x()*points[j].y()));
    }

    // Fill the values of points that are not in the triangulation (these are exact
    // values that are used to compare with approximations obtained through interpolation functions)
    for(std::size_t j=0; j<n+m; j++)
    {
      // linear function
      pt_values[0].insert(std::make_pair(points[j], alpha + beta1*points[j].x() + beta2*points[j].y()));

      // spherical function:
      pt_values[1].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                            beta2*points[j].y() +
                                                            gamma1*points[j].x()*points[j].x()+
                                                            gamma1*points[j].y()*points[j].y()));

      // quadratic function
      pt_values[2].insert(std::make_pair(points[j], alpha + beta1*points[j].x() +
                                                            beta2*points[j].y() +
                                                            gamma1*points[j].x()*points[j].x() +
                                                            gamma2*points[j].y()*points[j].y() +
                                                            gamma3*points[j].x()*points[j].y()));
    }
  }

  //INTERPOLATION OF RANDOM POINTS:
  Coord_type norm;
  Vertex_coordinate_vector vh_coords;
  for(std::size_t j=n; j<n+m; ++j)
  {
    CGAL::Triple<std::back_insert_iterator<Vertex_coordinate_vector>, Coord_type, bool> coordinate_result =
        CGAL::regular_neighbor_coordinates_2(T, points[j], std::back_inserter(vh_coords), vh_fct);
    assert(coordinate_result.third);

    norm = coordinate_result.second;
    assert(norm >= 0); // can be null if the weight was too small
    bool is_equal = test_norm(vh_coords.begin(), vh_coords.end(), norm);
    assert(is_equal);

    // further tests are pointless if there are no coordinates
    if(norm == 0)
      continue;

    // equivalent to a call to 'test_barycenter'
    Extract_point<Gt> cp;
    Bare_point b(CGAL::ORIGIN);
    VCV_cit ci = vh_coords.begin();
    for(; ci!=vh_coords.end(); ++ci)
    {
      Vertex_handle vh = ci->first;
      b = b + (ci->second / norm) * (cp(vh->point()) - CGAL::ORIGIN);
    }
    assert(CGAL::squared_distance(cp(points[j]), b) <= tolerance);

    for(int i=0; i<3; ++i)
    {
      std::pair<FT, bool> ev = CGAL::Data_access<Point_value_map>(pt_values[i])(points[j]);
      assert(ev.second);

      assert(test_interpolation_with_value(vh_coords.begin(), vh_coords.end(), norm, points[j],
                                           ev.first /*exact value*/,
                                           CGAL::Data_access<Vertex_value_map>(values[i]),
                                           CGAL::Data_access<Vertex_vector_map>(gradients[i]),
                                           Traits(), i, tolerance));
    }
    vh_coords.clear();
  }

  //TESTING THE GRADIENT APPRXIMATION METHOD:
    // sibson_gradient_fitting is tested when calling rn_2
  std::cout << "Testing gradient estimation method on random points." << std::endl;
  typedef CGAL::Interpolation_gradient_fitting_traits_2<Gt> GradTraits;
  Vertex_vector_map approx_gradients[2];
  CGAL::sibson_gradient_fitting_rn_2(T,
                                     std::inserter(approx_gradients[0], approx_gradients[0].begin()),
                                     CGAL::Identity<std::pair<Vertex_handle, Vector> >(),
                                     CGAL::Data_access<Vertex_value_map>(values[0]), GradTraits());

  CGAL::sibson_gradient_fitting_rn_2(T,
                                     std::inserter(approx_gradients[1], approx_gradients[1].begin()),
                                     CGAL::Identity<std::pair<Vertex_handle, Vector> >(),
                                     CGAL::Data_access<Vertex_value_map>(values[1]), GradTraits());

  typename Rt::Finite_vertices_iterator vit = T.finite_vertices_begin();
  for(; vit!=T.finite_vertices_end(); ++vit)
  {
    Vertex_handle vh = vit;
    assert(vh != Vertex_handle());
    std::pair<Vector, bool> res = CGAL::Data_access<Vertex_vector_map>(approx_gradients[0])(vh);

    if(res.second)
    {
      // if it is the exact computation kernel: test the equality:
      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Vertex_vector_map>(gradients[0])(vh).first);
      res = CGAL::Data_access<Vertex_vector_map>(approx_gradients[1])(vh);

      // if one exists->the other must also exist
      assert(res.second);
      assert(tolerance > Coord_type(0) ||
             res.first == CGAL::Data_access<Vertex_vector_map>(gradients[1])(vh).first);
    }
    else
    {
      assert(!CGAL::Data_access<Vertex_vector_map>(approx_gradients[1])(vh).second);
    }
  }

  //TESTING A POINT == A DATA POINT:
  Vertex_handle vh = T.finite_vertices_begin();
  std::cout << "testing at data point: " << vh->point() << std::endl;
  CGAL::Triple<std::back_insert_iterator<Vertex_coordinate_vector>, Coord_type, bool> coordinate_result =
      CGAL::regular_neighbor_coordinates_2(T, vh->point(), std::back_inserter(vh_coords), vh_fct);
  assert(coordinate_result.third);
  norm = coordinate_result.second;

  // if the point is hidden, the norm will be 0
  assert(norm == Coord_type(0) || norm == Coord_type(1));

  if(norm == Coord_type(1)) // nothing to do if the point is hidden
  {
    VCV_cit ci = vh_coords.begin();
    assert(ci->first == vh);
    assert(ci->second == Coord_type(1));
    ci++;
    assert(ci == vh_coords.end());

    for(int j=0; j<3; ++j)
    {
      std::pair<FT, bool> ev = CGAL::Data_access<Vertex_value_map>(values[j])(vh);
      assert(ev.second);

      assert(test_interpolation_with_value(vh_coords.begin(), vh_coords.end(), norm,
                                           vh->point(), ev.first /*exact value*/,
                                           CGAL::Data_access<Vertex_value_map>(values[j]),
                                           CGAL::Data_access<Vertex_vector_map>(gradients[j]),
                                           Traits(), j, tolerance));
    }
    vh_coords.clear();
  }
  else
  {
    assert(vh_coords.empty()); // point was hidden
  }
}

template <class Dt>
void _test_interpolation_functions_2_Delaunay(const Dt& dt, const typename Dt::Geom_traits::FT& tolerance)
{
  _test_interpolation_functions_2_Delaunay_without_OutputFunctor(dt, tolerance);
  _test_interpolation_functions_2_Delaunay_with_OutputFunctor(dt, tolerance);
}

template <class Rt>
void _test_interpolation_functions_2_regular(const Rt& rt, const typename Rt::Geom_traits::FT& tolerance)
{
  _test_interpolation_functions_2_regular_without_OutputFunctor(rt, tolerance);
  _test_interpolation_functions_2_regular_with_OutputFunctor(rt, tolerance);
}
