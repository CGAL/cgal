#ifndef CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H
#define CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H

#include <CGAL/aff_transformation_tags.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Interpolation_gradient_fitting_traits_2
//-----------------------------------------------------------------------//
// The class meets the requirement of the concept InterpolationTraits2
// and GradientFittingTraits2: it defines the geometric
// operations used by the interpolation methods and the gradient
// fitting function.

template<class Aff_2>
class Construct_sum_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Construct_sum_matrix_2(){};
  
  Aff_transformation_2
  operator()(const Aff_transformation_2& tr1,
	     const Aff_transformation_2& tr2){
    
    return(Aff_transformation_2(tr1.m(0,0) + tr2.m(0,0),
				 tr1.m(0,1) + tr2.m(0,1),
				 tr1.m(0,2) + tr2.m(0,2),
				 tr1.m(1,0) + tr2.m(1,0),
				 tr1.m(1,1) + tr2.m(1,1),
				 tr1.m(1,2) + tr2.m(1,2)));
  }

};

template<class Aff_2>
class Construct_null_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Construct_null_matrix_2(){};
  
  Aff_transformation_2
  operator()(){
    
    return Aff_transformation_2(0,0,0,0,0,0);
  }

};

template<class Aff_2>
class Construct_scaling_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Construct_scaling_matrix_2(){};
  
  Aff_transformation_2
  operator()(typename Aff_transformation_2::R::FT scale){
    
    return Aff_transformation_2(SCALING, scale);
  }

};

template< class  R >
class Construct_outer_product_2
{
public:
  typedef typename R::Aff_transformation_2       Aff_transformation_2;
  typedef typename R::Vector_2                   Vector_2;
  
  Construct_outer_product_2(){};
  
  Aff_transformation_2
  operator()(const Vector_2& v){  
    return(Aff_transformation_2(v.x()*v.x(),v.x()*v.y(),v.x()*v.y(),
				v.y()*v.y()));
  }

};


template <class R>
class Interpolation_gradient_fitting_traits_2 
{
public:
  typedef R                                          Rep;

  typedef typename Rep::FT                           FT;
  typedef typename Rep::Point_2                      Point;
  typedef typename Rep::Vector_2                     Vector;

  typedef typename Rep::Construct_vector_2           Construct_vector;
  typedef typename Rep::Construct_scaled_vector_2    Construct_scaled_vector;
  //only one not needed by gradient fitting:
  typedef typename Rep::Compute_squared_distance_2   Compute_squared_distance;
  
  
  //additional types for gradient computation:
  typedef typename Rep::Aff_transformation_2         Aff_transformation;

  typedef Construct_null_matrix_2<Aff_transformation>   
                                                     Construct_null_matrix;
  typedef Construct_scaling_matrix_2<Aff_transformation> 
                                                     Construct_scaling_matrix;
  typedef Construct_sum_matrix_2<Aff_transformation> Construct_sum_matrix;
  typedef Construct_outer_product_2<Rep>             Construct_outer_product;
  
  
  Construct_outer_product
  construct_outer_product_object() const
    {return Construct_outer_product();}
  
  Construct_sum_matrix
  construct_sum_matrix_object() const
    {return Construct_sum_matrix();}
  
  Construct_scaling_matrix
  construct_scaling_matrix_object() const
    {return Construct_scaling_matrix();}

  Construct_null_matrix
  construct_null_matrix_object() const
    {return Construct_null_matrix();}

  //also in the traits without gradient computation:
  Construct_scaled_vector
  construct_scaled_vector_object()const
    {return Construct_scaled_vector();}
  
  Construct_vector
  construct_vector_object()const
    {return Construct_vector();}
  
  Compute_squared_distance
  compute_squared_distance_object()const
    {return Compute_squared_distance();}

};
CGAL_END_NAMESPACE

#endif // CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H

