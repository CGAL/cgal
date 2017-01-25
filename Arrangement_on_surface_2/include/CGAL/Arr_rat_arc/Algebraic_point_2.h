// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_ALBERAIC_POINT_D_1_H
#define CGAL_ALBERAIC_POINT_D_1_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <ostream>

#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Cache.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Arr_rat_arc/Rational_function_canonicalized_pair.h>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Polynomial.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Bigfloat_interval_traits.h>

namespace CGAL {
namespace Arr_rational_arc {

//-------------------
//Algebraic_point_2_rep
//-------------------
template <typename AlgebraicKernelD_1>
class Algebraic_point_2_rep : public Base_rational_arc_ds_1<AlgebraicKernelD_1>
{
public:
  typedef AlgebraicKernelD_1                        Algebraic_kernel_d_1;
  typedef Base_rational_arc_ds_1<Algebraic_kernel_d_1>
                                                    Base;
 
  typedef CGAL::Arr_rational_arc::Rational_function<Algebraic_kernel_d_1>
                                                    Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_pair<Algebraic_kernel_d_1>
                                                    Rational_function_pair;
  
  typedef typename Base::Algebraic_real_1           Algebraic_real_1;
  typedef typename Algebraic_kernel_d_1::Bound      Bound;
  typedef typename Base::Integer                    Integer ;
  typedef typename Base::Rational                   Rational ;
  typedef typename Base::Polynomial_1               Polynomial_1;

  typedef typename Base::Root_multiplicity_vector   Root_multiplicity_vector;

  typedef typename Get_arithmetic_kernel<Algebraic_real_1>::Arithmetic_kernel
                                                    AK;
  typedef typename AK::Bigfloat_interval            BFI; 
  typedef Bigfloat_interval_traits<BFI>             BFI_traits;
  typedef CGAL::Polynomial<BFI>                     BFI_polynomial;
  typedef CGAL::Polynomial_traits_d<BFI_polynomial> BFI_polynomial_traits;

  
  typedef typename Base::FT_rat_1                   FT_rat_1; 
  typedef typename Base::Polynomial_traits_1        Polynomial_traits_1;

  typedef CGAL::Arr_rational_arc::Cache<Algebraic_kernel_d_1>
                                                    Cache;

public:
  Algebraic_point_2_rep(){}
  Algebraic_point_2_rep(const Rational_function& rational_function,
                        const Algebraic_real_1& x_coordinate) :
    _rational_function(rational_function),
    _x_coordinate(x_coordinate) {}

  //assignment oparator
  Algebraic_point_2_rep& operator=(const Algebraic_point_2_rep& other)
  {
    if (this != &other) // protect against invalid self-assignment
    {
      _rational_function = other._rational_function;
      _x_coordinate = other._x_coordinate;
    }
    return *this;
  }

  Comparison_result compare_xy_2(const Algebraic_point_2_rep& other,
                                 const Cache& cache) const
  {
    Comparison_result comp = CGAL::compare(_x_coordinate, other.x());
    if (comp != EQUAL)
      return comp;
    if (_rational_function == other.rational_function())
      return EQUAL;
    Rational_function_pair rat_func_pair =
      cache.get_rational_pair(_rational_function, other.rational_function());
    return rat_func_pair.compare_f_g_at(_x_coordinate);
  }

  Algebraic_real_1& x() 
  {
    return _x_coordinate;
  }

  const Algebraic_real_1& x() const 
  {
    return _x_coordinate;
  }

  const Rational_function& rational_function() const 
  {
    return _rational_function;
  }

  const Polynomial_1& numerator() const { return _rational_function.numer(); }
  const Polynomial_1& denominator() const { return _rational_function.denom(); }

  //new functions...
  Algebraic_real_1 y() const
  {
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    //converting the defining polynomial of x and the rational function to
    //bivariate polynomials
    Polynomial_2 f(_algebraic_kernel.compute_polynomial_1_object()(_x_coordinate));
   
    Polynomial_2 y(CGAL::shift(Polynomial_2(1),1));
    Polynomial_2 g(_rational_function.numer() - y * _rational_function.denom());

    f=CGAL::swap(f, 0, 1);  //swap x and y in the polynomial f
    g=CGAL::swap(g, 0, 1);  //swap x and y in the polynomial g
    //compute the resultant in x (polynomial in y)
    Polynomial_1 r(CGAL::resultant(f,g));
  
    //solve for all roots of resultant
    std::list<Algebraic_real_1> roots;
    _algebraic_kernel.solve_1_object()(r, false, std::back_inserter(roots));
    //isolate the right root
    unsigned int initial_precision = 16;
    int error_bound = 2;

    while (roots.size() > 1)
    {
      std::pair<Bound, Bound>
        y_bounds(this->approximate_absolute_y(error_bound,initial_precision));
      while (CGAL::compare(roots.front(),y_bounds.first) == SMALLER)
        roots.pop_front();
  
      while (CGAL::compare(y_bounds.second,roots.back()) == SMALLER)
        roots.pop_back();

      error_bound *= 2;
    }

    CGAL_postcondition (roots.size() == 1);
    return roots.front();
  }
  std::pair<double,double> to_double() const
  {
    double x = CGAL::to_double(_x_coordinate);
    double numer_val = evaluate_at(_rational_function.numer(),x);
    double denom_val = evaluate_at(_rational_function.denom(),x);
    return std::make_pair(x,numer_val/denom_val);
  }
  std::pair<Bound,Bound> approximate_absolute_x( int a) const
  {
    return _algebraic_kernel.approximate_absolute_1_object()(_x_coordinate,a);
  }
  std::pair<Bound,Bound> approximate_absolute_y( int a) const
  {
    unsigned int precision = 16;
    return approximate_absolute_y(a,precision);
  }
  std::pair<Bound,Bound> approximate_relative_x( int r) const
  {    
    return _algebraic_kernel.approximate_relative_1_object()(_x_coordinate,r);
  }
  std::pair<Bound,Bound> approximate_relative_y( int r) const
  {
    // return zero if y is zero (the code below would not terminate)
    // moreover approx relative is actually not well defined 
    if(_rational_function.sign_at(_x_coordinate)==CGAL::ZERO)
      return std::make_pair(Bound(0),Bound(0));
    
    typename BFI_traits::Set_precision       set_precision;
    typename BFI_polynomial_traits::Evaluate evaluate;
    
    typedef typename BFI_traits::Bound    BF;

    long precision = 16;
    set_precision(precision);
    BF eps = CGAL::ipower(BF(1)/2,r);
    
    while (true){      
      set_precision(precision);
      BFI x_bfi(convert_to_bfi(_x_coordinate));
      
      BFI_polynomial
        numer_bfi(convert_to_bfi_extended(_rational_function.numer()));
      BFI_polynomial
        denom_bfi(convert_to_bfi_extended(_rational_function.denom()));

      BFI y_numer_bfi(evaluate(numer_bfi,x_bfi));
      BFI y_denom_bfi(evaluate(denom_bfi,x_bfi));
      
      if (CGAL::zero_in(y_denom_bfi) == false)
      {
        BFI y_bfi(y_numer_bfi/y_denom_bfi);
       
        if (CGAL::compare( 
                CGAL::width(y_bfi),
                CGAL::lower(CGAL::abs(y_bfi)) * eps)
            == SMALLER)
          return std::make_pair(
              Bound(CGAL::lower(y_bfi)),
              Bound(CGAL::upper(y_bfi)));
      }
      else precision*=2;
    }
  }

  std::ostream& print (std::ostream& os) const
  {
    std::pair<double,double> double_p;
    switch(::CGAL::get_mode(os))
    {
     case ::CGAL::IO::PRETTY:
      double_p = this->to_double();
      os <<"(" ;
      os << double_p.first;
      os <<" , " ;
      os << double_p.second;
      os << ")";
      break;
          
     case ::CGAL::IO::BINARY:
      std::cerr << "BINARY format not yet implemented" << std::endl;
      break;
          
     default:
      // ASCII
      os <<"x = " << _x_coordinate<<" ";
      os <<"rational function = ( " ;
      Base::print_polynomial(os, _rational_function.numer(), 'x');
      os << ") / (";
      Base::print_polynomial(os, _rational_function.denom(), 'x');
      os << ")";
    }
     
    return os;
 }
private:
  std::pair<Bound,Bound> approximate_absolute_y( int a,unsigned int& precision ) const
  {
    typename BFI_traits::Set_precision       set_precision;
    typename BFI_polynomial_traits::Evaluate evaluate;
    
    typedef typename BFI_traits::Bound             BF;
    
    BF eps = CGAL::ipower(BF(1)/2,a);
    while (true)
    {
      set_precision(precision);
      BFI x_bfi(convert_to_bfi(_x_coordinate));

      BFI_polynomial numer_bfi(convert_to_bfi_extended(_rational_function.numer()));
      BFI_polynomial denom_bfi(convert_to_bfi_extended(_rational_function.denom()));

      BFI y_numer_bfi(evaluate(numer_bfi,x_bfi));
      BFI y_denom_bfi(evaluate(denom_bfi,x_bfi));
      
      if (CGAL::zero_in(y_denom_bfi) == false)
      {
        BFI y_bfi(y_numer_bfi/y_denom_bfi);
        if (CGAL::width(y_bfi) < eps )
          return std::make_pair(
              Bound(CGAL::lower(y_bfi)),
              Bound(CGAL::upper(y_bfi)));
       
      }
      else precision*=2;
    }
  }

  template <typename NTX>
  static typename
  CGAL::Coercion_traits<NTX,
                        typename Get_arithmetic_kernel<NTX>::
                          Arithmetic_kernel::Bigfloat_interval>::Type
  convert_to_bfi_extended(const NTX& x) 
  {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval BFI; 
    typedef CGAL::Coercion_traits<NTX, BFI> CT;
    return typename CT::Cast()(x);
  }

  double evaluate_at(const Polynomial_1& poly, const double x) const
  {
    double x_val = 1;
    double ret_val(0);
    for (int i(0); i <= poly.degree(); ++i)
    {
      ret_val = ret_val + x_val*CGAL::to_double(poly[i]);
      x_val = x_val*x;
    }
    return ret_val;
  }

private:
  Rational_function _rational_function;  //supporting rational function
  Algebraic_real_1  _x_coordinate;

  Algebraic_kernel_d_1 _algebraic_kernel;
};

template <typename Algebraic_kernel_ >
class Algebraic_point_2 :
    public Handle_with_policy<Algebraic_point_2_rep<Algebraic_kernel_> >
{
 
public:
  typedef Algebraic_kernel_                               Algebraic_kernel_d_1;
  typedef Handle_with_policy<Algebraic_point_2_rep<Algebraic_kernel_> >
                                                          Base;
  typedef Algebraic_point_2<Algebraic_kernel_d_1>         Self;
  typedef Algebraic_point_2_rep<Algebraic_kernel_>        Rep;
  typedef typename Rep::Rational                          Rational;
  typedef typename Rep::Algebraic_real_1                  Algebraic_real_1;
  typedef typename Rep::Rational_function                 Rational_function;
  typedef typename Rep::Bound                             Bound;
  typedef typename Rep::Cache                             Cache;
  typedef typename Rep::Polynomial_1                      Polynomial_1; 

private:
  static Self& get_default_instance()
  {
    static Algebraic_kernel_d_1 kernel;
    static typename Rational_function::Polynomial_1 numer(0);
    static typename Rational_function::Polynomial_1 denom(1);
    static Rational_function rational_function(numer, denom, &kernel);
    
    static Algebraic_real_1 x_coordinate =
      kernel.construct_algebraic_real_1_object()(Rational(0));
    
    static Self default_instance(rational_function,x_coordinate); 
    
    return default_instance;

    /*static Self x = Self(Rational(0),Rational(0),_cache); 
    return x; */
  }

public:
  explicit Algebraic_point_2(const Rational_function& rational_function,
                             const Algebraic_real_1& x_coordinate) :
    Base(rational_function,x_coordinate) {}
  
  Algebraic_point_2() :
    Base(static_cast<const Base &> (get_default_instance())) {}

  // explicit copy-constructor, required by VC9
  Algebraic_point_2 (const Self & p)
    : Base(static_cast<const Base &> (p)) {}

  Comparison_result compare_xy_2(const Algebraic_point_2& other,
                                 const Cache& cache) const
  {
    if (this->is_identical (other))
      return CGAL::EQUAL;
    return this->ptr()->compare_xy_2(*other.ptr(), cache);
  }

  const Polynomial_1& numerator() const { return this->ptr()->numerator(); }
  const Polynomial_1& denominator() const { return this->ptr()->denominator(); }
  
  Algebraic_real_1& x()
  {
    if (this->is_shared())
      this->copy_on_write();
    return this->ptr()->x();
  }

  const Algebraic_real_1& x() const 
  {
    return this->ptr()->x();
  }

  const Rational_function& rational_function() const 
  {
    return this->ptr()->rational_function();
  }

  Algebraic_real_1 y() const
  {
    return this->ptr()->y();
  }

  std::pair<double,double> to_double() const
  {
    return this->ptr()->to_double();
  }

  std::pair<Bound,Bound> approximate_absolute_x( int a) const
  {
    return this->ptr()->approximate_absolute_x(a);
  }

  std::pair<Bound,Bound> approximate_absolute_y( int a) const
  {
    return this->ptr()->approximate_absolute_y(a);
  }

  std::pair<Bound,Bound> approximate_relative_x( int r) const
  {    
    return this->ptr()->approximate_relative_x(r);
  }

  std::pair<Bound,Bound> approximate_relative_y( int r) const
  {
    return this->ptr()->approximate_relative_y(r);
  }

  std::ostream& print(std::ostream& os) const
  {
    return this->ptr()->print(os);
  }
};  //Algebraic_point_2


template < typename Algebraic_kernel_>
std::ostream& operator<<(std::ostream& os,
                         const Algebraic_point_2<Algebraic_kernel_> & p)
{
  return (p.print(os));
}

}   //namespace Arr_rational_arc
}   //namespace CGAL {   

#endif //CGAL_ALBERAIC_POINT_D_1_H
