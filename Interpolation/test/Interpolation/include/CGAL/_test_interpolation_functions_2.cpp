// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Julia Floetotto

#include <iostream>
#include <cassert>
#include <utility>

#include <CGAL/double.h>
#include <CGAL/Random.h>
#include <CGAL/Origin.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

#include <CGAL/natural_neighbor_coordinates_2.h>

#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>


template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
	       typename std::iterator_traits<ForwardIterator>::
	       value_type::second_type norm)
{
  typename std::iterator_traits<ForwardIterator>::value_type::second_type
           sum(0);
   for(; first !=beyond; first++)
     sum+= first->second;

   return norm==sum;
}

template < class ForwardIterator >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
		     typename std::iterator_traits<ForwardIterator>::
		     value_type::second_type norm,
		     const typename std::iterator_traits<ForwardIterator>::
		     value_type::first_type& p, const typename
		     std::iterator_traits<ForwardIterator>::
		     value_type::second_type& tolerance)
{
  typedef typename
    std::iterator_traits<ForwardIterator>::value_type::first_type Point;
  Point b(CGAL::ORIGIN);
  for(; first !=beyond; first++)
    b = b+ (first->second/norm) * (first->first - CGAL::ORIGIN);
  return CGAL::squared_distance(p,b) <= tolerance;
}

/////////////////////////////////////////////////////////////////
// Accessory function testing functions that require sqrt().
// Doesn't instantiate anything if RT doesn't support sqrt().
template < class ForwardIterator, class Functor, class GradFunctor,
           class Gt>
bool
_test_sibson_c1_interpolation_sqrt(ForwardIterator ,
				   ForwardIterator ,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type&
				   /* norm */, 
                                   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::first_type& /* p */,
				   Functor /* f */, 
                                   GradFunctor /* grad_f */,
				   const Gt& /*geom_traits */,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type& /* tolerance */,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type& /* exact_value */,
                                   CGAL::Integral_domain_without_division_tag)
{
//  bool UNTESTED_STUFF_BECAUSE_SQRT_IS_NOT_SUPPORTED;
  std::cout << std::endl
            << "NOTE: FT doesn't support sqrt(),"
    " hence sibson_c1_interpolation is not tested." << std::endl;
  return true;
}

template < class ForwardIterator, class Functor, class GradFunctor,
           class Gt>
bool
_test_sibson_c1_interpolation_sqrt(ForwardIterator first,
				   ForwardIterator beyond,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type&
				   norm, const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type
				   ::first_type& p,
				   Functor f, GradFunctor grad_f,
				   const Gt& geom_traits,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type&
				   tolerance,
				   const typename
				   std::iterator_traits<ForwardIterator>
				   ::value_type::second_type& exact_value,
				   CGAL::Field_with_sqrt_tag)
{
  typename Functor::result_type res =
    CGAL::sibson_c1_interpolation(first, beyond,
				  norm,p,f,
				  grad_f, geom_traits);
  return res.second && (CGAL_NTS abs(res.first-exact_value)<= tolerance);
}

template < class ForwardIterator, class Functor, class GradFunctor,
           class Gt>
bool test_interpolation(ForwardIterator first, ForwardIterator beyond,
			const typename
			std::iterator_traits<ForwardIterator>
			::value_type::second_type&
			norm, const typename
			std::iterator_traits<ForwardIterator>::value_type
			::first_type& p,
			Functor f, GradFunctor grad_f,
			const Gt& geom_traits,
			const int& i,
			const typename std::iterator_traits<ForwardIterator>
			::value_type::second_type& tolerance)
{
  typedef typename Functor::result_type::first_type Value_type;
  assert(f(p).second);
  Value_type exact_value = f(p).first;

  if(i==0){
    Value_type val =  CGAL::linear_interpolation(first, beyond, norm,f);
    assert(CGAL_NTS abs(val-exact_value)<= tolerance);
  }

  typename Functor::result_type
    res =  CGAL::quadratic_interpolation(first, beyond, norm,p,f,
					 grad_f, geom_traits);
  assert(res.second && (CGAL_NTS abs(res.first-exact_value)<=tolerance));

  if(i<2){
    //without sqrt:
    res =  CGAL::sibson_c1_interpolation_square(first, beyond,
						norm,p,f,
						grad_f, geom_traits);
    assert(res.second && (CGAL_NTS abs(res.first-exact_value)<= tolerance));

    //with sqrt:
    typedef CGAL::Algebraic_structure_traits<Value_type> AST;
    
    assert(_test_sibson_c1_interpolation_sqrt(
                   first, beyond, norm,p,f, grad_f, 
                   geom_traits, tolerance, exact_value,
                   typename AST::Algebraic_category()));
  }
  res =  CGAL::farin_c1_interpolation(first, beyond, norm,p,f,grad_f,
				      geom_traits);
  assert(res.second);
  assert(CGAL_NTS abs(res.first-exact_value)<= tolerance);
  return true;
}

template <class Triangul>
void
_test_interpolation_functions_2_delaunay( const Triangul &,
		  const typename Triangul::Geom_traits::FT& tolerance)
{
  CGAL::Set_ieee_double_precision pfr;
  Triangul T;

  int n=20, m=20;
  double r = 3;
  double max_value = 5;

  typedef typename Triangul::Geom_traits          Gt;
  typedef CGAL::Interpolation_traits_2<Gt>        Traits;

  typedef typename Triangul::Face_handle         Face_handle;

  typedef typename Gt::Point_2                    Point;
  typedef typename Gt::FT                         Coord_type;
  typedef typename Gt::Vector_2                   Vector;

  typedef std::map<Point, Coord_type, typename Gt::Less_xy_2> Point_value_map ;
  typedef std::map<Point, Vector, typename Gt::Less_xy_2 > Point_vector_map;

  typedef std::vector< std::pair<Point, Coord_type> > Point_coordinate_vector;

  std::cout << "NN2: Testing random points." << std::endl;
  //test random points in a square of length r:
  std::vector<Point> points;
  points.reserve(n+m);

  //put four bounding box points:
  points.push_back(Point(-r,-r));
  points.push_back(Point(r,-r));
  points.push_back(Point(-r,r));
  points.push_back(Point(r,r));

  // Create n+m-4 points within a disc of radius 2
  CGAL::Random_points_in_square_2<Point> g(r);
  CGAL::cpp11::copy_n( g, n+m, std::back_inserter(points));

  CGAL::Random random;

  Point_value_map values[3];
  Point_vector_map gradients[3];

  Coord_type alpha = Coord_type(random.get_double(-max_value,max_value)),
    beta1 = Coord_type(random.get_double(-max_value,max_value)),
    beta2 = Coord_type(random.get_double(-max_value,max_value)),
    gamma1 = Coord_type(random.get_double(-max_value,max_value)),
    gamma2 = Coord_type(random.get_double(-max_value,max_value)),
    gamma3 = Coord_type(random.get_double(-max_value,max_value));

  //INSERTION + DET. of GRADIENT for n DATA POINTS :
  for(int j=0; j<n; j++){
    T.insert(points[j]);

    gradients[0].insert(std::make_pair(points[j], Vector( beta1, beta2)));

    gradients[1].insert(std::make_pair
			(points[j],
			 Vector( beta1+Coord_type(2)*gamma1*points[j].x(),
				 beta2+Coord_type(2)*gamma1* points[j].y())));
    gradients[2].insert
      (std::make_pair(points[j],Vector(beta1+Coord_type(2)*gamma1
				       *points[j].x()+gamma3* points[j].y(),
				       beta2 +
				       Coord_type(2) * gamma2* points[j].y()+
				       gamma3* points[j].x())));
  }

  //DETERMINE FUNCTION VALUE FOR n DATA POINTS AND m RANDOM TEST POINTS:
  for(int j=0; j<n+m; j++){
    //linear function
    values[0].insert(std::make_pair(points[j],alpha +
				    beta1*points[j].x()
				    + beta2*points[j].y()));

    //spherical function:
    values[1].insert(std::make_pair(points[j],alpha +
				    beta1*points[j].x() +
				    beta2*points[j].y() +
				    gamma1*points[j].x()*points[j].x()+
				    gamma1*points[j].y()*points[j].y()));

    //quadratic function
    values[2].insert(std::make_pair(points[j],alpha +
				    beta1*points[j].x() +
				    beta2*points[j].y() +
				    gamma1*points[j].x()*points[j].x()+
				    gamma2*points[j].y()*points[j].y()
				    + gamma3*points[j].x()*points[j].y()));
  }

  //INTERPOLATION OF RANDOM POINTS:
  Coord_type norm;
  Point_coordinate_vector coords;
  for(int j=n;j<n+m;j++){
    CGAL::Triple<
      std::back_insert_iterator<Point_coordinate_vector>,
      Coord_type, bool> coordinate_result =
      CGAL::natural_neighbor_coordinates_2(T, points[j],
					   std::back_inserter(coords));
    assert(coordinate_result.third);
    norm = coordinate_result.second;

    assert(norm>0);
    assert(test_norm( coords.begin(), coords.end(),norm));
    assert(test_barycenter( coords.begin(), coords.end(),norm,
			    points[j], tolerance));

    for(int i=0; i<3; i++)
      assert(test_interpolation(coords.begin(), coords.end(),norm,points[j],
       				CGAL::Data_access< Point_value_map >
				(values[i]),
       				CGAL::Data_access< Point_vector_map >
				(gradients[i]),
       				Traits(), i, tolerance ));
    coords.clear();
  }

  //TESTING THE GRADIENT APPRXIMATION METHOD:
  std::cout << "Testing gradient estimation method on random points."
	    <<std::endl;
  typedef CGAL::Interpolation_gradient_fitting_traits_2<Gt> GradTraits;
  Point_vector_map approx_gradients[2];
  CGAL::sibson_gradient_fitting_nn_2(T,std::inserter(approx_gradients[0],
						     approx_gradients[0].
						     begin()),
				     CGAL::Data_access<Point_value_map>
				     (values[0]),
				     GradTraits());
  CGAL::sibson_gradient_fitting_nn_2
    (T,std::inserter(approx_gradients[1],approx_gradients[1].begin()),
     CGAL::Data_access<Point_value_map>(values[1]),GradTraits());
  for(int j=0; j<n; j++){
    std::pair<Vector, bool> res =
      CGAL::Data_access<Point_vector_map>(approx_gradients[0])(points[j]);
    if(res.second){

      //if it is the exact computation kernel: test the equality:
      assert(tolerance > Coord_type(0) || res.first ==
	     CGAL::Data_access<Point_vector_map>(gradients[0])
	     (points[j]).first);
      res =
	CGAL::Data_access<Point_vector_map>(approx_gradients[1])(points[j]);
      //if one exists->the other must also exist
      assert(res.second);


      assert(tolerance > Coord_type(0) || res.first ==
	     CGAL::Data_access<Point_vector_map>(gradients[1])
	     (points[j]).first);
    }else
      assert(!CGAL::Data_access<Point_vector_map>(approx_gradients[1])
	     (points[j]).second);
  }


  //TESTING A POINT == A DATA POINT:
  CGAL::Triple<
    std::back_insert_iterator<Point_coordinate_vector>,
    Coord_type, bool> coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T,points[n/2],std::back_inserter
					 (coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  typename std::vector< std::pair< Point, Coord_type > >::iterator
    ci= coords.begin();
  assert(ci->first == points[n/2]);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci==coords.end());
  for(int j=0; j<3; j++) {
    assert(test_interpolation(coords.begin(), coords.end(),norm,points[n/2],
			      CGAL::Data_access< Point_value_map >(values[j]),
			      CGAL::Data_access< Point_vector_map >
			      (gradients[j]),
			      Traits(),j, tolerance));
  }
  coords.clear();

  //FURTHER TESTS FOR NATURAL NEIGHBOR COORDINATES 2:
  //TESTING A POINT on an EDGE of the triangulation:
  Face_handle fh = T.finite_faces_begin();
  int i =0;
  if(T.is_infinite(fh->neighbor(i)))
    i++;
  assert(!T.is_infinite(fh->neighbor(i)));
  Point p = fh->vertex(T.ccw(i))->point() + Coord_type(0.5)*
    (fh->vertex(T.cw(i))->point()-fh->vertex(T.ccw(i))->point());

  coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T, p,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm =  coordinate_result.second;
  assert(test_norm( coords.begin(), coords.end(),norm));
  assert(test_barycenter( coords.begin(), coords.end(),norm,p, tolerance));
  coords.clear();
  //END OF TEST WITH EDGE


  //TESTING a GRID POINT SET
  std::cout << "NN2: Testing grid points." << std::endl;

  Triangul T2;
  //Grid points:
  Point p1_2(-2, -2);
  Point p2_2(-2,2);
  Point p3_2(2,-2);
  Point p4_2(2,2);
  Point p1(-1, -1);
  Point p2(-1,1);
  Point p3(1,-1);
  Point p4(1,1);
  Point p12(-1, 0);
  Point p23(0,1);
  Point p34(1,0);
  Point p41(0,-1);

  T2.insert(p1_2);T2.insert(p2_2);T2.insert(p3_2);T2.insert(p4_2);
  T2.insert(p1);T2.insert(p2);T2.insert(p3);T2.insert(p4);
  T2.insert(p12);T2.insert(p23);T2.insert(p34);T2.insert(p41);

  coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T2,Point(0,0),
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  ci= coords.begin();
  for(; ci!= coords.end(); ci++)
    assert(ci->second == Coord_type(0.25));
  assert(test_barycenter( coords.begin(), coords.end(),norm,
			  Point(0,0), tolerance));
  coords.clear();

  //add the middle point of the grid
  T2.insert(Point(0,0));

  //point on a vertex;
  coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T2,p34,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  ci= coords.begin();
  assert(ci->first == p34);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci==coords.end());
  coords.clear();

  //point on an edge:
  p= Point(0,0.5);
  coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T2,p,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm, p, tolerance));
  coords.clear();

  //Point outside convex hull:
  p= Point(3,0.5);
  coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T2,p,
					 std::back_inserter(coords));
  assert(!coordinate_result.third);
  //Point on a convex hull edge:
  p= Point(2,1);
   coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T2,p,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
}
