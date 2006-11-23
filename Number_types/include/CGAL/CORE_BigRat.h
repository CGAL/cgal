// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Hemmer  <mhemmer@uni-mainz.de>
#ifndef CGAL_CORE_BIGRAT_H
#define CGAL_CORE_BIGRAT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/CORE_coercion_traits.h>

#include <CGAL/CORE_Expr.h> // used for To_interval-functor
#include <CGAL/CORE_BigInt.h>


#define CORE_LEVEL 4
#include <CORE/CORE.h>


//#if defined(CGAL_CORE_BIGRAT_NUMER_DENOM_ARE_MEMBERS)
//  #define CGAL_CORE_NUMERATOR(X) ((X).numerator())
//  #define CGAL_CORE_DENOMINATOR(X) ((X).denominator())
//#elif defined(CGAL_CORE_BIGRAT_NUMER_DENOM_ARE_NONMEMBERS)
  #define CGAL_CORE_NUMERATOR(X) (numerator((X)))
  #define CGAL_CORE_DENOMINATOR(X) (denominator((X)))
//#else

CGAL_BEGIN_NAMESPACE

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigRat >
  : public Algebraic_structure_traits_base< CORE::BigRat, 
                                            Field_tag >  {
  public:
    typedef Tag_true            Is_exact;

    // BigRat are always normalized, so no special simplify-functor is needed
    
    // Nothing new...                
};




//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigRat > 
  : public Real_embeddable_traits_base< CORE::BigRat > {
  public:
      
    class Abs 
      : public Unary_function< Real_embeddable, Real_embeddable > {
      public:        
        Real_embeddable operator()( const Real_embeddable& x ) const {
          return CORE::abs( x );
        }
    };
    
    class Sign 
      : public Unary_function< Real_embeddable, ::CGAL::Sign > {
      public:        
        ::CGAL::Sign operator()( const Real_embeddable& x ) const {
          return (::CGAL::Sign) CORE::sign( x );
        }        
    };
    
    class Compare 
      : public Binary_function< Real_embeddable, Real_embeddable,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Real_embeddable& x, 
                                            const Real_embeddable& y ) const {
          typedef Real_embeddable_traits<int> Int_traits;
          return Int_traits::Sign()( ::CORE::cmp(x,y));
        }
    };
    
    class To_double 
      : public Unary_function< Real_embeddable, double > {
      public:
        double operator()( const Real_embeddable& x ) const {
          // this call is required to get reasonable values for the double
          // approximation 
          return x.doubleValue();
        }
    };
    
    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {

          Real_embeddable_traits<CORE::Expr>::To_interval to_interval;
          CORE::Expr temp(x);
            
          return to_interval(temp); 
        }          
    };
};

/*! \ingroup NiX_Fraction_traits_spec
 *  \brief Specialization of Fraction_traits for ::leda::rational
 */
template <>
class Fraction_traits< CORE::BigRat > {
public:
    typedef CORE::BigRat Fraction;
    typedef ::CGAL::Tag_true Is_fraction;
    typedef CORE::BigInt Numerator;
    typedef Numerator Denominator;

    typedef Algebraic_structure_traits< Numerator >::Gcd Common_factor;

    class Decompose {
    public:
        typedef Fraction first_argument_type;
        typedef Numerator& second_argument_type;
        typedef Numerator& third_argument_type;
        void operator () (
                const Fraction& rat,
                Numerator& num,
                Numerator& den) {
            num = CGAL_CORE_NUMERATOR(rat);
            den = CGAL_CORE_DENOMINATOR(rat);
        }
    };
    
    class Compose {
    public:
        typedef Numerator first_argument_type;
        typedef Numerator second_argument_type;
        typedef Fraction result_type;
        Fraction operator ()(
                const Numerator& num , 
                const Numerator& den ) {
            return Fraction(num, den);
        }
    };
};

template <class F>
class Output_rep< ::CORE::BigRat, F> {
    const ::CORE::BigRat& t;
public:
    //! initialize with a const reference to \a t.
    Output_rep( const ::CORE::BigRat& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const {
        switch (get_mode(out)) {
        case IO::PRETTY:{
            if(CGAL_CORE_DENOMINATOR(t) == ::CORE::BigRat(1))
                return out <<CGAL_CORE_NUMERATOR(t);
            else
                return out << CGAL_CORE_NUMERATOR(t)
                           << "/" 
                           << CGAL_CORE_DENOMINATOR(t);
            break;
        }
            
        default:
            return out << CGAL_CORE_NUMERATOR(t)
                       << "/" 
                       << CGAL_CORE_DENOMINATOR(t);
        }
    }
};

template <>
struct Needs_parens_as_product< ::CORE::BigRat >{
    bool operator()( ::CORE::BigRat t){
        if (CGAL_CORE_DENOMINATOR(t) != 1 ) 
            return true;
        else
            return needs_parens_as_product(CGAL_CORE_NUMERATOR(t)) ;
    }
};

template <>
class Output_rep< ::CORE::BigRat, Parens_as_product_tag > {
    const ::CORE::BigRat& t;
public:
    // Constructor 
    Output_rep( const ::CORE::BigRat& tt) : t(tt) {}
    // operator 
    std::ostream& operator()( std::ostream& out) const { 
        Needs_parens_as_product< ::CORE::BigRat > needs_parens_as_product;
        if (needs_parens_as_product(t))
            return out <<"("<< oformat(t) <<")";
        else
            return out << oformat(t);
    }
};




CGAL_END_NAMESPACE

#endif // CGAL_CORE_BIGRAT_H
