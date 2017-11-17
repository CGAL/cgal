// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de> 
//              
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Polynomial_traits_d.h>

#ifndef CGAL_POLYNOMIAL_UTILS_H
#define CGAL_POLYNOMIAL_UTILS_H

#define CGAL_UNARY_POLY_FUNCTION(functor,function)                      \
  template <typename Polynomial_d>  inline                              \
  typename Polynomial_traits_d<Polynomial_d>::functor::result_type      \
  function(const Polynomial_d& p){                                      \
    typedef Polynomial_traits_d<Polynomial_d> PT;                       \
    return typename PT::functor()(p);                                   \
  }   

#define CGAL_UNARY_POLY_FUNCTION_INDEX(functor,function)                \
  CGAL_UNARY_POLY_FUNCTION(functor,function)                            \
  template <typename Polynomial_d>  inline                              \
  typename Polynomial_traits_d<Polynomial_d>::functor::result_type      \
  function(const Polynomial_d& p, int index ){                          \
    typedef Polynomial_traits_d<Polynomial_d> PT;                       \
    typename PT::functor fo;                                            \
    return fo(p,index);                                                 \
  }      

#define CGAL_BINARY_POLY_FUNCTION(functor,function)                     \
  template <typename Polynomial_d>  inline                              \
  typename Polynomial_traits_d<Polynomial_d>::functor::result_type      \
  function(const Polynomial_d& p,                                       \
      const typename Polynomial_traits_d<Polynomial_d>::                \
      functor::second_argument_type& second                             \
  ){                                                                    \
    typedef Polynomial_traits_d<Polynomial_d> PT;                       \
    return typename PT::functor()(p,second);                            \
  }  

#define CGAL_BINARY_POLY_FUNCTION_INDEX(functor,function)               \
  CGAL_BINARY_POLY_FUNCTION(functor,function)                           \
  template <typename Polynomial_d>  inline                              \
  typename Polynomial_traits_d<Polynomial_d>::functor::result_type      \
  function(const Polynomial_d& p,                                       \
      const typename Polynomial_traits_d<Polynomial_d>::                \
      functor::second_argument_type& second,                            \
      int index                                                         \
  ){                                                                    \
    typedef Polynomial_traits_d<Polynomial_d> PT;                       \
    return typename PT::functor()(p,second,index);                      \
  }  

namespace CGAL {

// GetCoefficient
template <typename Polynomial_d> inline  
typename Polynomial_traits_d<Polynomial_d>::Get_coefficient::result_type
get_coefficient(const Polynomial_d& p, int i){
  typename Polynomial_traits_d<Polynomial_d>::Get_coefficient get_coefficient;
  return get_coefficient(p,i);
}
// GetInnermostCoefficient
template <typename Polynomial_d> inline  
typename Polynomial_traits_d<Polynomial_d>
::Get_innermost_coefficient::result_type
get_innermost_coefficient(const Polynomial_d& p, Exponent_vector ev){
  typename Polynomial_traits_d<Polynomial_d>::Get_innermost_coefficient gic;
  return gic(p,ev);
}
// ConstructCoefficientConstIteratorRange
// ConstructInnermostCoefficientConstIteratorRange
// Swap
template <typename Polynomial_d> inline  
typename Polynomial_traits_d<Polynomial_d>::Swap::result_type
swap(const Polynomial_d& p, int i, int j){
  typename Polynomial_traits_d<Polynomial_d>::Swap swap;
  return swap(p,i,j);
}
// Move
template <typename Polynomial_d> inline  
typename Polynomial_traits_d<Polynomial_d>::Move::result_type
move(const Polynomial_d& p, int i, int j){
  typename Polynomial_traits_d<Polynomial_d>::Move move;
  return move(p,i,j);
}
// Permute
template <typename Polynomial_d, typename Input_iterator> inline  
typename Polynomial_traits_d<Polynomial_d>::Permute::result_type
permute(const Polynomial_d& p, Input_iterator begin, Input_iterator end){
  typename Polynomial_traits_d<Polynomial_d>::Permute permute;
  return permute(p,begin,end);
}

// Degree
CGAL_UNARY_POLY_FUNCTION_INDEX(Degree,degree)
// TotalDegree
CGAL_UNARY_POLY_FUNCTION(Total_degree,total_degree)
// DegreeVector
CGAL_UNARY_POLY_FUNCTION(Degree_vector,degree_vector)
// LeadingCoefficient
CGAL_UNARY_POLY_FUNCTION(Leading_coefficient,leading_coefficient)
// InnermostLeadingCoefficient
CGAL_UNARY_POLY_FUNCTION(
    Innermost_leading_coefficient,
    innermost_leading_coefficient)

// Canonicalize
CGAL_UNARY_POLY_FUNCTION(Canonicalize, canonicalize)
// Differentiate
CGAL_UNARY_POLY_FUNCTION_INDEX(Differentiate, differentiate)

// Evaluate
CGAL_BINARY_POLY_FUNCTION(Evaluate,evaluate)
// EvaluateHomogeneous
template <typename Polynomial_d, typename T>  inline                                 
typename Polynomial_traits_d<Polynomial_d>::Evaluate_homogeneous::result_type
evaluate_homogeneous(const Polynomial_d& p,const T& num, const T& den){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Evaluate_homogeneous()(p,num,den);                      
}  

// Substitute
template <typename Polynomial_d, typename Input_iterator>  inline     
typename CGAL::Coercion_traits<
    typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type,
    typename std::iterator_traits<Input_iterator>::value_type>              
  ::Type  
substitute(const Polynomial_d& p,Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Substitute()(p,begin, end);                      
}  
// IsZeroAt
template <typename Polynomial_d, typename Input_iterator>  inline     
typename Polynomial_traits_d<Polynomial_d>::Is_zero_at::result_type
is_zero_at(const Polynomial_d& p, Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Is_zero_at()(p,begin, end);                      
}  
// SignAt
template <typename Polynomial_d, typename Input_iterator>  inline     
typename Polynomial_traits_d<Polynomial_d>::Sign_at::result_type
sign_at(const Polynomial_d& p, Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Sign_at()(p,begin, end);                      
}  

// SubstituteHomogeneous
template <typename Polynomial_d, typename Input_iterator>  inline     
typename CGAL::Coercion_traits<
    typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type,
    typename std::iterator_traits<Input_iterator>::value_type>              
  ::Type  
substitute_homogeneous(
    const Polynomial_d& p,Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Substitute_homogeneous()(p,begin, end);
}  
// IsZeroAtHomogeneous
template <typename Polynomial_d, typename Input_iterator>  inline     
typename Polynomial_traits_d<Polynomial_d>::Is_zero_at_homogeneous::result_type
is_zero_at_homogeneous(
    const Polynomial_d& p, Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Is_zero_at_homogeneous()(p,begin, end);
}  
// SignAtHomogeneous
template <typename Polynomial_d, typename Input_iterator>  inline     
typename Polynomial_traits_d<Polynomial_d>::Sign_at_homogeneous::result_type
sign_at_homogeneous(
    const Polynomial_d& p, Input_iterator begin, Input_iterator end){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Sign_at_homogeneous()(p,begin, end);                      
}  

// Compare // provided by number_type utils 
// CGAL_BINARY_POLY_FUNCTION(Compare,compare);

// UnivariateContent
CGAL_UNARY_POLY_FUNCTION(Univariate_content, univariate_content)
// MultivariateContent
CGAL_UNARY_POLY_FUNCTION(Multivariate_content, multivariate_content)

// SquareFreeFactorize
template <typename Polynomial_d, typename OutputIterator>  inline     
OutputIterator
square_free_factorize(const Polynomial_d& p, OutputIterator oi){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Square_free_factorize()(p,oi);                      
}  
// MakeSquareFree
CGAL_UNARY_POLY_FUNCTION(Make_square_free, make_square_free)
// IsSquareFree
CGAL_UNARY_POLY_FUNCTION(Is_square_free, is_square_free)

// PseudoDivision
// PseudoDivisionQuotient
// PseudoDivisionRemainder
template <typename Polynomial_d>  inline     
void 
pseudo_division(
    const Polynomial_d& f, const Polynomial_d& g, 
    Polynomial_d& q, Polynomial_d& r, 
    typename Polynomial_traits_d<Polynomial_d>::Coefficient_type& D){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  typename PT::Pseudo_division()(f,g,q,r,D);                      
  return;
}  
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Pseudo_division_quotient::result_type
pseudo_division_quotient(const Polynomial_d& f, const Polynomial_d& g){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Pseudo_division_quotient()(f,g);                      
}  
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Pseudo_division_remainder::result_type
pseudo_division_remainder(const Polynomial_d& f, const Polynomial_d& g){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Pseudo_division_remainder()(f,g);                      
}  

// GcdUpToConstantFactor
CGAL_BINARY_POLY_FUNCTION(
    Gcd_up_to_constant_factor, 
    gcd_up_to_constant_factor)
// IntegralDivisionUpToConstantFactor
CGAL_BINARY_POLY_FUNCTION(
    Integral_division_up_to_constant_factor, 
    integral_division_up_to_constant_factor)
// UnivariateContentUpToConstantFactor
CGAL_UNARY_POLY_FUNCTION(
    Univariate_content_up_to_constant_factor, 
    univariate_content_up_to_constant_factor)
// SquareFreeFactorizeUpToConstantFactor
template <typename Polynomial_d, typename OutputIterator>  inline     
OutputIterator
square_free_factorize_up_to_constant_factor(
    const Polynomial_d& p, OutputIterator oi){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Square_free_factorize_up_to_constant_factor()(p,oi);
}  

// Shift
CGAL_BINARY_POLY_FUNCTION_INDEX(Shift,shift)
// Negate
CGAL_UNARY_POLY_FUNCTION_INDEX(Negate,negate)
// Invert
CGAL_UNARY_POLY_FUNCTION_INDEX(Invert,invert)
// Translate
CGAL_BINARY_POLY_FUNCTION_INDEX(Translate,translate)
// TranslateHomogeneous
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Translate_homogeneous::result_type
translate_homogeneous(const Polynomial_d& f,
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& num, 
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& den){
      
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Translate_homogeneous()(f,num,den);                      
}  
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Translate_homogeneous::result_type
translate_homogeneous(const Polynomial_d& f,
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& num, 
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& den,
    int index ){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Translate_homogeneous()(f,num,den,index);
}      
// Scale
CGAL_BINARY_POLY_FUNCTION_INDEX(Scale,scale)
// ScaleHomogeneous
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Scale_homogeneous::result_type
scale_homogeneous(const Polynomial_d& f,
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& num, 
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& den){
      
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Scale_homogeneous()(f,num,den);                      
}  
template <typename Polynomial_d>  inline     
typename Polynomial_traits_d<Polynomial_d>::Scale_homogeneous::result_type
scale_homogeneous(const Polynomial_d& f,
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& num, 
    const typename Polynomial_traits_d<Polynomial_d>::Innermost_coefficient_type& den,
    int index ){
  typedef Polynomial_traits_d<Polynomial_d> PT;                       
  return typename PT::Scale_homogeneous()(f,num,den,index);
}  
// Resultant
CGAL_BINARY_POLY_FUNCTION(Resultant,resultant)

template <typename Polynomial_d,typename OutputIterator> inline
OutputIterator polynomial_subresultants
(Polynomial_d p, Polynomial_d q, OutputIterator out) {
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Polynomial_subresultants() (p, q, out);
}   

template <typename Polynomial_d,typename OutputIterator> inline
OutputIterator principal_subresultants
(Polynomial_d p, Polynomial_d q, OutputIterator out) {
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Principal_subresultants() (p, q, out);
}   

template<typename Polynomial_d,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3> inline
OutputIterator1 polynomial_subresultants_with_cofactors
(Polynomial_d p,
 Polynomial_d q,
 OutputIterator1 sres_out,
 OutputIterator2 coP_out,
 OutputIterator3 coQ_out) {
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Polynomial_subresultants_with_cofactors() 
        (p, q, sres_out, coP_out, coQ_out);
}


template <typename Polynomial_d,typename OutputIterator> inline
OutputIterator
principal_sturm_habicht_sequence
(Polynomial_d f, OutputIterator out){
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Principal_sturm_habicht_sequence() (f, out);
}
  
template<typename Polynomial_d,typename OutputIterator> OutputIterator
sturm_habicht_sequence(Polynomial_d f,OutputIterator out) {
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Sturm_habicht_sequence() (f, out);
}

template<typename Polynomial_d,
    typename OutputIterator1,
    typename OutputIterator2,
    typename OutputIterator3> 
OutputIterator1
sturm_habicht_sequence_with_cofactors
(Polynomial_d f,
 OutputIterator1 stha_out,
 OutputIterator2 cof_out,
 OutputIterator3 cofx_out) {
    typedef Polynomial_traits_d<Polynomial_d> PT;
    return typename PT::Sturm_habicht_sequence_with_cofactors() 
        (f, stha_out, cof_out, cofx_out);
}


// TODO: REMOVE function below ?

template<typename NT> inline
Polynomial<NT> scale_up(const Polynomial<NT>& p, const NT& a)
{ Polynomial<NT> q(p); q.scale_up(a); return q; }

template<typename NT> inline
Polynomial<NT> scale_down(const Polynomial<NT>& p, const NT& b)
{ Polynomial<NT> q(p); q.scale_down(b); return q; }


template<typename NT> inline
Polynomial<NT> translate_by_one(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.translate_by_one(); return q; }


template<typename NT> inline
Polynomial<NT> reversal(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.reversal(); return q; }


} //namespace CGAL

#undef CGAL_UNARY_POLY_FUNCTION
#undef CGAL_UNARY_POLY_FUNCTION_INDEX
#undef CGAL_BINARY_POLY_FUNCTION
#undef CGAL_BINARY_POLY_FUNCTION_INDEX

#endif // CGAL_POLYNOMIAL_UTILS_H
