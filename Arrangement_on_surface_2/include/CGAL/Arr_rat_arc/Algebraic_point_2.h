#ifndef CGAL_ALBERAIC_POINT_D_1_H
#define CGAL_ALBERAIC_POINT_D_1_H

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
  typedef typename Base::Polynomial_1       Polynomial_1;

  typedef typename Base::Root_multiplicity_vector Root_multiplicity_vector;

  typedef typename Get_arithmetic_kernel<Algebraic_real_1>::Arithmetic_kernel AK;
  typedef typename AK::Bigfloat_interval                                      BFI; 
  typedef Bigfloat_interval_traits<BFI>                                       BFI_traits;
  typedef CGAL::Polynomial<BFI>                                               BFI_polynomial;
  typedef CGAL::Polynomial_traits_d<BFI_polynomial>                           BFI_polynomial_traits;

  
  typedef typename Base::FT_rat_1               FT_rat_1; 
  typedef typename Base::Polynomial_traits_1    Polynomial_traits_1;

  typedef CGAL::Arr_rational_arc::Cache<Kernel_,Algebraic_kernel_>      Cache;

public:
  Algebraic_point_2_rep(){}
  Algebraic_point_2_rep(const Rational_function& rational_function,const Algebraic_real_1& x_coordinate):
    _rational_function(rational_function),_x_coordinate(x_coordinate) {}
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
  Comparison_result compare_xy_2 (const Algebraic_point_2_rep & other,Cache& cache) const
  {
    Comparison_result comp = CGAL::compare (_x_coordinate, other.x());
    if (comp != EQUAL)
      return comp;
    if (this->_rational_function == other.rational_function())
      return EQUAL;
    Rational_function_pair rat_func_pair = cache.get_rational_pair( this->_rational_function,
                                                                    other.rational_function()); 
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


  //new functions...
  Algebraic_real_1 y() const
  {
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    //converting the defining polynomial of x and the rational function to bivariate polynomials
    Polynomial_2 f(_algebraic_kernel.compute_polynomial_1_object()(_x_coordinate));
   
    Polynomial_2 y(CGAL::shift(Polynomial_2(1),1));
    Polynomial_2 g(_rational_function.numer() - y * _rational_function.denom());

    f=CGAL::swap(f,0,1);  //swap x and y in the polynomial f
    g=CGAL::swap(g,0,1);  //swap x and y in the polynomial g
    Polynomial_1 r(CGAL::resultant(f,g)); //compute the resultant in x (polynomial in y)
  
    //solve for all roots of resultant
    std::list<Algebraic_real_1> roots;
    _algebraic_kernel.solve_1_object()(r,false,std::back_inserter(roots));
    //isolate the right root
    unsigned int initial_precision = 16;
    int error_bound = 2;

    while (roots.size() > 1)
    {
      std::pair<Bound,Bound> y_bounds (this->approximate_absolute_y(error_bound,initial_precision));
      while (CGAL::compare(roots.front(),y_bounds.first) == SMALLER)
        roots.pop_front();
  
      while (CGAL::compare(y_bounds.second,roots.back()) == SMALLER)
        roots.pop_back();

      error_bound*=2;
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
    BFI_traits::Set_precision       set_precision;
    BFI_polynomial_traits::Evaluate evaluate;
    long precision = 16;
    Rational error_bound = CGAL::ipower(Rational(1,2),r);
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
        if (CGAL::compare (CGAL::width(y_bfi),Rational(CGAL::lower(CGAL::abs(y_bfi)))*error_bound ) == SMALLER)
          return std::make_pair(Bound(CGAL::lower(y_bfi)),
                                Bound(CGAL::upper(y_bfi)) );
      }
      else precision*=2;
    }
  }

  std::ostream& print (std::ostream& os) const
  {
    os <<"x = " << _x_coordinate<<" ";
    os <<"rational function = ( " ;
    Base::print_polynomial(os, _rational_function.numer(), 'x');
    os << ") / (";
    Base::print_polynomial(os, _rational_function.denom(), 'x');
    os << ")";

    return os;
  }
private:
  std::pair<Bound,Bound> approximate_absolute_y( int a,unsigned int& precision ) const
  {
    BFI_traits::Set_precision       set_precision;
    BFI_polynomial_traits::Evaluate evaluate;
    Rational error_bound = CGAL::ipower(Rational(1,2),a);
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
        if (Bound(CGAL::width(y_bfi)) < error_bound )
          return std::make_pair(Bound(CGAL::lower(y_bfi)),
                                Bound(CGAL::upper(y_bfi)) );
      }
      else precision*=2;
    }
  }
  template <class NTX>
  static typename CGAL::Coercion_traits<NTX,typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval>::Type
  convert_to_bfi_extended(const NTX& x) 
  {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval BFI; 
    typedef CGAL::Coercion_traits<NTX,BFI> CT;
    return typename CT::Cast()(x);
  }
  double evaluate_at(const Polynomial_1& poly,const double x) const
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

  Algebraic_kernel _algebraic_kernel;
};



template < class Kernel_, 
   class Algebraic_kernel_ = Algebraic_kernel_d_1 <typename Fraction_traits <typename Kernel_::FT>::Numerator_type> >
class Algebraic_point_2: public Handle_with_policy<Algebraic_point_2_rep<Kernel_,Algebraic_kernel_> >
{
 
public:
  typedef Kernel_                                             Kernel;
  typedef Algebraic_kernel_                                   Algebraic_kernel;
  typedef Handle_with_policy<Algebraic_point_2_rep<Kernel_,Algebraic_kernel_> > Base;
  typedef Algebraic_point_2<Kernel,Algebraic_kernel>          Self;
  typedef Algebraic_point_2_rep<Kernel_,Algebraic_kernel_>    Rep;
  typedef typename Rep::Rational                              Rational;
  typedef typename Rep::Algebraic_real_1                      Algebraic_real_1;
  typedef typename Rep::Rational_function                     Rational_function;
  typedef typename Rep::Bound                                 Bound;
  typedef typename Rep::Cache                                 Cache;

private:
  
  static Self& get_default_instance()
  {
    static Rational_function::Polynomial_1 numer(0);
    static Rational_function::Polynomial_1 denom(1);
    static Rational_function rational_function(numer,denom);

    static Algebraic_kernel kernel;
    static Algebraic_real_1 x_coordinate = kernel.construct_algebraic_real_1_object()(Rational(0));
    
    static Self default_instance(rational_function,x_coordinate); 
    
    return default_instance;

    /*static Self x = Self(Rational(0),Rational(0),_cache); 
    return x; */
  }
public:
  explicit Algebraic_point_2( const Rational_function& rational_function,
                              const Algebraic_real_1& x_coordinate) 
                              : Base(rational_function,x_coordinate) {}
  //used to solve VS bug...
  Algebraic_point_2 (const Self & p = get_default_instance()) : Base(static_cast<const Base &> (p)) {}


  Comparison_result compare_xy_2 (const Self & other,Cache& cache) const
  {
    if (this->is_identical (other))
      return CGAL::EQUAL;
    return this->ptr()->compare_xy_2(*other.ptr(),cache);
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

  std::ostream& print (std::ostream& os) const
  {
    return this->ptr()->print(os);
  }
};  //Algebraic_point_2


template < class Kernel_,
   class Algebraic_kernel_  >
std::ostream&
operator<< (std::ostream& os,
            const Algebraic_point_2<Kernel_, Algebraic_kernel_> & p)
{
  return (p.print (os));
}
}   //namespace Arr_rational_arc
}   //namespace CGAL {   
#endif //CGAL_ALBERAIC_POINT_D_1_H
