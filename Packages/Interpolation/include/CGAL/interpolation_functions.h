// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source:
// /CVSROOT/CGAL/Packages/Interpolation/include/CGAL/interpolation_functions.h,v $
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Julia Floetotto
#ifndef CGAL_INTERPOLATION_FUNCTIONS_H
#define CGAL_INTERPOLATION_FUNCTIONS_H

#include <utility>
#include <CGAL/double.h>

CGAL_BEGIN_NAMESPACE 


template < class ForwardIterator, class Functor>
typename Functor::result_type 
linear_interpolation(ForwardIterator first, ForwardIterator beyond,
		     const typename
		     std::iterator_traits<ForwardIterator>::value_type::
		     second_type& norm, Functor f)
{  
  CGAL_precondition(norm>0);
  typedef typename Functor::result_type Value_type;
  Value_type result(0);
  for(; first !=beyond; first++)
    result += (first->second/norm) * f(first->first);
  
  return result;
};


template < class ForwardIterator, class Functor, class GradFunctor,class Gt>
typename Functor::result_type 
quadratic_interpolation(ForwardIterator first, ForwardIterator beyond,  
			const typename
			std::iterator_traits<ForwardIterator>::
			value_type::second_type& norm, const typename
			std::iterator_traits<ForwardIterator>::value_type::
			first_type& p,
			Functor f, GradFunctor grad_f,
			const Gt& geom_traits)
{  
  CGAL_precondition(norm >0);
  typedef typename Functor::result_type Value_type;
  Value_type result(0);
  
  for(; first !=beyond; first++)
    result += (first->second/norm) 
      *( f(first->first) + grad_f(first->first)*
	 geom_traits.construct_scaled_vector_object()
	 (geom_traits.construct_vector_object()(first->first, p),0.5));
  
  return result;
};


template < class ForwardIterator, class Functor, class GradFunctor, class Gt>
typename Functor::result_type 
sibson_c1_interpolation(ForwardIterator first, ForwardIterator beyond,
			const typename
			std::iterator_traits<ForwardIterator>::
			value_type::second_type&
			norm, const typename
			std::iterator_traits<ForwardIterator>::value_type::
			first_type& p, 
			Functor f, GradFunctor grad_f, 
			const Gt& geom_traits)
{  
  CGAL_precondition(norm >0);
  typedef typename Functor::result_type Value_type;
  typedef typename Gt::FT   Coord_type;

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0),gradient_int(0);
  
  for(; first !=beyond; first++){
    
    Coord_type coeff = first->second/norm;
    Coord_type squared_dist = geom_traits.
      compute_squared_distance_object()(first->first, p);
    Coord_type dist = CGAL_NTS sqrt(squared_dist);
    
    if(squared_dist ==0){
      ForwardIterator it = first;
      CGAL_assertion(it++==beyond);
      return f(first->first);
    }
    //three different terms to mix linear and gradient
    //interpolation
    term1 +=  coeff/dist;
    term2 +=  coeff * squared_dist;
    term3 +=  coeff * dist;
    
    linear_int += coeff * f(first->first);
    
    gradient_int += (coeff/dist)
      *(f(first->first) + grad_f(first->first)*
	geom_traits.construct_vector_object()(first->first, p));
  }
      
  term4 = term3/ term1;
  gradient_int = gradient_int / term1;
 
  return (term4* linear_int + term2 * gradient_int)/
    (term4 + term2);
};

//this method works with rational number types:
//modification of Sibson's interpolant without sqrt
//following a proposition by Gunther Rote:
//
// the general scheme:
//  Coord_type inv_weight = f(dist); //i.e. dist^2
//   	term1 +=  coeff/inv_weight;
//	term2 +=  coeff * squared_dist;
//	term3 +=  coeff*(squared_dist/inv_weight);  
// 	gradient_int +=  (coeff/inv_weight)*
// 	  (vh->get_value()+  vh->get_gradient()
// 	   *(p - vh->point()));

template < class ForwardIterator, class Functor, class GradFunctor,
  class Gt>
typename Functor::result_type 
sibson_c1_interpolation_square(ForwardIterator first, ForwardIterator beyond,
			       const typename
			       std::iterator_traits<ForwardIterator>::
			       value_type::second_type& norm, const typename
			       std::iterator_traits<ForwardIterator>::
			       value_type::first_type& p, 
			       Functor f, GradFunctor grad_f, 
			       const Gt& geom_traits)
{  
  CGAL_precondition(norm >0);
  typedef typename Functor::result_type Value_type;
  typedef typename Gt::FT               Coord_type;

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0),gradient_int(0);

  for(; first !=beyond; first++){
    
    Coord_type coeff = first->second/norm;
    Coord_type squared_dist = geom_traits.
      compute_squared_distance_object()(first->first, p);
    
    if(squared_dist ==0){
      ForwardIterator it = first;
      CGAL_assertion(++it==beyond);
      return f(first->first);
    }
    //three different terms to mix linear and gradient
    //interpolation
    term1 +=  coeff/squared_dist;
    term2 +=  coeff * squared_dist;
    term3 +=  coeff;
    
    linear_int += coeff * f(first->first);
    
    gradient_int += (coeff/squared_dist)
      *(f(first->first) + grad_f(first->first)*
	geom_traits.construct_vector_object()(first->first, p));
  }
      
  term4 = term3/ term1;
  gradient_int = gradient_int / term1;
  
  
  return (term4* linear_int + term2 * gradient_int)/
    (term4 + term2);
};


template < class RandomAccessIterator, class Functor, class
GradFunctor, class Gt>
typename Functor::result_type 
farin_c1_interpolation(RandomAccessIterator first,
		       RandomAccessIterator beyond,
		       const typename 
		       std::iterator_traits<RandomAccessIterator>::
		       value_type::second_type& norm, const typename
		       std::iterator_traits<RandomAccessIterator>::
		       value_type::first_type& p, 
		       Functor f, GradFunctor grad_f, const Gt& geom_traits)
{ 
  CGAL_precondition(norm >0);
  typedef typename Functor::result_type  Value_type;
  typedef typename Gt::FT                Coord_type;
  
 
  int n= beyond - first;
  if( n==1)
    return f(first->first);

  //there must be one or at least three NN-neighbors:
  CGAL_assertion(n > 2);
   
  RandomAccessIterator it2, it;
  
  Value_type result(0);
  const Coord_type fac3(3);

  std::vector< std::vector<Value_type> > 
    ordinates(n,std::vector<Value_type>(n, Value_type(0)));
  
  for(int i =0; i<n; i++){
    it = first+i;
    Coord_type coord_i_square = CGAL_NTS square(it->second);
    
    //for later:
    ordinates[i][i] =  f(it->first);
    
    //control point = data point
    result += coord_i_square * it->second* ordinates[i][i];
    
    
    //compute tangent plane control point (one 2, one 1 entry)
    Value_type res_i(0);
    for(int j =0; j<n; j++){
      if(i!=j){
	it2 = first+j;
	//  ordinates[i][j] = (p_j - p_i) * g_i
	ordinates[i][j] = grad_f(it->first) * 
	  geom_traits.construct_vector_object()(it->first,it2->first);
	
	// a point in the tangent plane: 
	// 3( f(p_i) + (1/3)(p_j - p_i) * g_i) 
	// => 3*f(p_i) + (p_j - p_i) * g_i
	res_i += (fac3 * ordinates[i][i] + ordinates[i][j])* it2->second;
      }
    }
    //res_i already multiplied by three
    result += coord_i_square *res_i;
  }
  
  //the third type of control points: three 1 entries i,j,k
  for(int i=0; i< n; i++)
    for(int j=i+1; j< n; j++)
      for(int k=j+1; k<n; k++){
	// add 6* (u_i*u_j*u_k) * b_ijk
	//  b_ijk = 1.5 * a - 0.5*c 
	//where 
	//c : average of the three data control points
	//a : 1.5*a = 1/12 * (ord[i][j] + ord[i][k] + ord[j][i] +
	//            ord[j][k] + ord[k][i]+ ord[k][j])
	// =>  6 * b_ijk = 3*(f_i + f_j + f_k) + 0.5*a 
	result += (Coord_type(2.0)*( ordinates[i][i]+ ordinates[j][j]+
				     ordinates[k][k])
		   + Coord_type(0.5)*(ordinates[i][j] + ordinates[i][k] 
				      + ordinates[j][i] +
				      ordinates[j][k] + ordinates[k][i]+
				      ordinates[k][j]))
	  *(first+i)->second *(first+j)->second *(first+k)->second ;
      }
    
 
  
  return(result/(CGAL_NTS square(norm)*norm)); 
}

CGAL_END_NAMESPACE

#endif // CGAL_INTERPOLATION_FUNCTIONS_H
