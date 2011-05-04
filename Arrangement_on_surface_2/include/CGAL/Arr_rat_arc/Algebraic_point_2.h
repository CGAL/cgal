#ifndef CGAL_ALBERAIC_POINT_D_1_H
#define CGAL_ALBERAIC_POINT_D_1_H

#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Cache.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Arr_rat_arc/Rational_function_canonicalized_pair.h>
#include <CGAL/Arr_rat_arc/Singleton.h>

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
template < class Kernel_, 
   class Algebraic_kernel_ = Algebraic_kernel_d_1 <typename Fraction_traits <typename Kernel_::FT>::Numerator_type> >
class Algebraic_point_2_rep: public Base_rational_arc_ds_1<Kernel_, Algebraic_kernel_>
{
public:
  typedef Kernel_                                           Kernel;
  typedef Algebraic_kernel_                                 Algebraic_kernel;
  typedef Base_rational_arc_ds_1<Kernel_, Algebraic_kernel> Base;
 
  typedef CGAL::Arr_rational_arc::Rational_function<Kernel_, Algebraic_kernel>        Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_pair<Kernel_,Algebraic_kernel>    Rational_function_pair;
  
  typedef typename Base::Algebraic_real_1   Algebraic_real_1;
  typedef typename Algebraic_kernel::Bound  Bound;
  typedef typename Base::Integer            Integer ;
  typedef typename Base::Rational           Rational ;

  typedef typename Get_arithmetic_kernel<Algebraic_real_1>::Arithmetic_kernel AK;
  typedef typename AK::Bigfloat_interval                                      BFI; 
  typedef Bigfloat_interval_traits<BFI>                                       BFI_traits;
  typedef CGAL::Polynomial<BFI>                                                     BFI_polynomial;
  typedef CGAL::Polynomial_traits_d<BFI_polynomial>                           BFI_polynomial_traits;

  
  typedef typename Base::FT_rat_1               FT_rat_1; 
  typedef typename Base::Polynomial_traits_1    Polynomial_traits_1;

  typedef CGAL::Arr_rational_arc::Cache<Kernel_,Algebraic_kernel_>      Cache;
  typedef Singleton<Cache>                                              Rational_functions_cache;

public:
  Algebraic_point_2_rep(){}
  Algebraic_point_2_rep(const Rational_function& rational_function,const Algebraic_real_1& x_coordinate):
    _rational_function(rational_function),_x_coordinate(x_coordinate) {}
  Algebraic_point_2_rep(const Rational& x,const Rational& y)
  {
    Algebraic_kernel kernel;
    _x_coordinate = kernel.construct_algebraic_real_1_object()(x);

    Integer  y_numer,y_denom;
    typename FT_rat_1::Decompose()(y,y_numer,y_denom);
  
    Cache* cache = Rational_functions_cache::instance();
    _rational_function =  cache->get_rational_function (Rational(y_numer  , y_denom ));
  }
  Algebraic_point_2_rep(const Algebraic_real_1& x,const Rational& y)
    :_x_coordinate(x)
  {
    Integer  y_numer,y_denom;
    typename FT_rat_1::Decompose()(y,y_numer,y_denom);
  
    Cache* cache = Rational_functions_cache::instance();
    _rational_function =  cache->get_rational_function (Rational(y_numer  , y_denom ));
  }

  //assignment oparator
  Algebraic_point_2_rep & operator = (const Algebraic_point_2_rep & other)
  {
    if (this != &other) // protect against invalid self-assignment
      {
        _rational_function = other._rational_function;
        _x_coordinate = other._x_coordinate;
      }
    return *this;
  }
  Comparison_result compare_xy_2 (const Algebraic_point_2_rep & other) const
  {
    Comparison_result comp = CGAL::compare (_x_coordinate, other.x());
    if (comp != EQUAL)
      return comp;
    if (this->_rational_function == other.rational_function())
      return EQUAL;
    Cache* cache = Rational_functions_cache::instance();
    Rational_function_pair rat_func_pair = cache->get_rational_pair(this->_rational_function,
        other.rational_function()); 
    return rat_func_pair.compare_f_g_at(_x_coordinate);
  }
  bool operator == (const Algebraic_point_2_rep & other) const
  {
    return (compare_xy_2 (other) == EQUAL);
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


  //new functions...
  Algebraic_real_1 y() const
  {
  }
  std::pair<double,double> to_double() const
  {
  }
  std::pair<Bound,Bound> approximate_absolute_x( int a) const
  {
    return _algebraic_kernel.approximate_absolute_1_object()(_x_coordinate,a);
  }
  std::pair<Bound,Bound> approximate_absolute_y( int a) const
  {
    
    BFI_traits::Set_precision       set_precision;
    BFI_polynomial_traits::Evaluate evaluate;
    long precision = 16;
    Rational error_bound = CGAL::ipower(Rational(1,2),a);
    while (true)
    {
      set_precision(precision);
      BFI x_bfi(convert_to_bfi(_x_coordinate));

      BFI_polynomial numer_bfi(convert_to_bfi(_rational_function.numer()));
      BFI_polynomial denom_bfi(convert_to_bfi(_rational_function.denom()));

      BFI y_numer_bfi(evaluate(numer_bfi,x_bfi));
      BFI y_denom_bfi(evaluate(denom_bfi,x_bfi));
      
      if (CGAL::zero_in(y_denom_bfi) == false)
      {
        BFI y_bfi(y_numer_bfi/y_denom_bfi);
        if (CGAL::width(y_bfi) < error_bound )
          return std::make_pair(Bound(CGAL::lower(y_bfi)),
                                Bound(CGAL::upper(y_bfi)) );
      }
      else precision*=2;
    }
  }
  std::pair<Bound,Bound> approximate_relative_x( int r) const
  {    
    return _algebraic_kernel.approximate_relative_1_object()(_x_coordinate,a);
  }
  std::pair<Bound,Bound> approximate_relative_y( int r) const
  {
    BFI_traits::Set_precision       set_precision;
    BFI_polynomial_traits::Evaluate evaluate;
    long precision = 16;
    Rational error_bound = CGAL::ipower(Rational(1,2),r);
    while (true)
    {
      set_precision(precision);
      BFI x_bfi(convert_to_bfi(_x_coordinate));

      BFI_polynomial numer_bfi(convert_to_bfi(_rational_function.numer()));
      BFI_polynomial denom_bfi(convert_to_bfi(_rational_function.denom()));

      BFI y_numer_bfi(evaluate(numer_bfi,x_bfi));
      BFI y_denom_bfi(evaluate(denom_bfi,x_bfi));
      
      if (CGAL::zero_in(y_denom_bfi) == false)
      {
        BFI y_bfi(y_numer_bfi/y_denom_bfi);
        if (CGAL::width(y_bfi) < Rational(CGAL::lower(CGAL::absolute(y_bfi)))*error_bound )
          return std::make_pair(Bound(CGAL::lower(y_bfi)),
                                Bound(CGAL::upper(y_bfi)) );
      }
      else precision*=2;
    }
  }

private:
  Rational_function _rational_function;  //supporting rational function
  Algebraic_real_1  _x_coordinate;

  Algebraic_kernel _algebraic_kernel;
};



template < class Kernel_, 
   class Algebraic_kernel_ = Algebraic_kernel_d_1 <typename Fraction_traits <typename Kernel_::FT>::Numerator_type> >
class Algebraic_point_2: public Handle_with_policy<Algebraic_point_2_rep<Kernel_,Algebraic_kernel_> >
{
 
public:
  typedef Kernel_            Kernel;
  typedef Algebraic_kernel_         Algebraic_kernel;
  typedef Handle_with_policy<Algebraic_point_2_rep<Kernel_,Algebraic_kernel_> > Base;
  typedef Algebraic_point_2<Kernel,Algebraic_kernel>          Self;
  typedef Algebraic_point_2_rep<Kernel_,Algebraic_kernel_>    Rep;
  typedef typename Rep::Rational                              Rational;
  typedef typename Rep::Algebraic_real_1                      Algebraic_real_1;
  typedef typename Rep::Rational_function                     Rational_function;

  //  Algebraic_point_2 () :Base() {}
private:
  static Self& get_default_instance(){
    static Self x = Self(Rational(0),Rational(0)); 
    return x; 
  }
public:
  template <typename _T> 
  explicit Algebraic_point_2(const _T& t) : Base(t) {}

  /*template <typename _T1, typename _T2> 
  explicit Algebraic_point_2(const _T1& t1 ,const _T2& t2) : Base(t1,t2) {}*/

  explicit Algebraic_point_2( const Rational_function& rational_function,
                              const Algebraic_real_1& x_coordinate) 
                              : Base(rational_function,x_coordinate) {}
  explicit Algebraic_point_2( const Rational& x,
                              const Rational& y) 
                              : Base(x,y) {}
  explicit Algebraic_point_2( const Algebraic_real_1& x,
                              const Rational& y) 
                              : Base(x,y) {}

  //used to solve VS bug...
  Algebraic_point_2 (const Self & p = get_default_instance()) : Base(static_cast<const Base &> (p)) {}


  Comparison_result compare_xy_2 (const Self & other) const
  {
    if (this->is_identical (other))
      return CGAL::EQUAL;
    return this->ptr()->compare_xy_2(*other.ptr());
  }
    
  bool operator== (const Self & other) const
  {
    if (this->is_identical (other))
      return true;
    return (*(this->ptr()) == *(other.ptr()));
  }

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

};  //Algebraic_point_2

}   //namespace Arr_rational_arc
}   //namespace CGAL {   
#endif //CGAL_ALBERAIC_POINT_D_1_H
