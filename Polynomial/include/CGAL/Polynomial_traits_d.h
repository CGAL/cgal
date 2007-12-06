// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//
// ============================================================================
#ifndef CGAL_POLYNOMIAL_TRAITS_D_H
#define CGAL_POLYNOMIAL_TRAITS_D_H

#include <CGAL/basic.h>


#include <CGAL/Polynomial/polynomial_utils.h>
#include <CGAL/Polynomial/resultant.h>
#include <CGAL/Polynomial/square_free_factorization.h>

#define CGAL_POLYNOMIAL_TRAITS_D_BASE_TYPEDEFS                          \
    typedef Polynomial_traits_d< Polynomial< Coefficient_ > > PT;       \
    typedef Polynomial_traits_d< Coefficient_ > PTC;                    \
                                                                        \
    public:                                                             \
    typedef Polynomial<Coefficient_>                  Polynomial_d;     \
    typedef Coefficient_                              Coefficient;      \
                                                                        \
    typedef typename Innermost_coefficient<Polynomial_d>::Type          \
    Innermost_coefficient;                                              \
    static const int d = Dimension<Polynomial_d>::value;                \
                                                                        \
                                                                        \
    private:                                                            \
    typedef std::pair< Exponent_vector, Innermost_coefficient >         \
    Exponents_coeff_pair;                                               \
    typedef std::vector< Exponents_coeff_pair > Monom_rep;              \
                                                                        \
    typedef CGAL::Recursive_const_flattening< d-1,                      \
    typename CGAL::Polynomial<Coefficient>::const_iterator >            \
    Coefficient_flattening;                                             \
                                                                        \
    public:                                                             \
    typedef typename Coefficient_flattening::Recursive_flattening_iterator \
    Innermost_coefficient_iterator;                                     \
    typedef typename  Polynomial_d::iterator Coefficient_iterator;      \
                                                                        \
    private:

CGAL_BEGIN_NAMESPACE;

namespace POLYNOMIAL {

// template meta function Innermost_coefficient
// returns the tpye of the innermost coefficient 
template <class T> struct Innermost_coefficient{ typedef T Type; };
template <class Coefficient> 
struct Innermost_coefficient<Polynomial<Coefficient> >{
    typedef typename Innermost_coefficient<Coefficient>::Type Type; 
};

// template meta function Dimension
// returns the number of variables 
template <class T> struct Dimension{ static const int value = 0;};
template <class Coefficient> 
struct Dimension<Polynomial<Coefficient> > {
    static const int value = Dimension<Coefficient>::value + 1 ; 
};

// Base class for functors depending on the algebraic category of the
// innermost coefficient
template< class Coefficient_, class ICoeffAlgebraicCategory >
class Polynomial_traits_d_base_icoeff_algebraic_category {    
    public:    
        typedef Null_functor    Multivariate_content;
};

// Specializations
template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Integral_domain_without_division_tag > 
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Null_tag > {};

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Integral_domain_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Integral_domain_without_division_tag > {}; 

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Unique_factorization_domain_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Integral_domain_tag > {
    CGAL_POLYNOMIAL_TRAITS_D_BASE_TYPEDEFS
    
    public:    
    
    //       Multivariate_content;
    struct Multivariate_content
        : public Unary_function< Polynomial_d , Innermost_coefficient >{
        Innermost_coefficient 
        operator()(const Polynomial_d& p) const {
            typedef Innermost_coefficient_iterator IT;
            Innermost_coefficient content(0);
            for (IT it = typename PT::Innermost_coefficient_begin()(p);
                 it != typename PT::Innermost_coefficient_end()(p);
                 it++){
                content = CGAL::gcd(content, *it);
                if(CGAL::is_one(content)) break;
            }
            return content;
        }
    };
}; 

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Euclidean_ring_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Unique_factorization_domain_tag > {}; 

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Field_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Integral_domain_tag > {
    CGAL_POLYNOMIAL_TRAITS_D_BASE_TYPEDEFS

    public:
    
    //       Multivariate_content;
    struct Multivariate_content
        : public Unary_function< Polynomial_d , Innermost_coefficient >{
        Innermost_coefficient 
        operator()(const Polynomial_d& p) const {
            typename PT::Compare compare;
            if( compare( p, Polynomial_d(0) ) == EQUAL )
                return Innermost_coefficient(0);
            else
                return Innermost_coefficient(1);
        }
    };
};

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Field_with_sqrt_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Field_tag > {}; 

template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Field_with_kth_root_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Field_with_sqrt_tag > {}; 
 
 template< class Coefficient_ >
class Polynomial_traits_d_base_icoeff_algebraic_category< 
            Polynomial< Coefficient_ >, Field_with_root_of_tag >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
                Polynomial< Coefficient_ >, Field_with_kth_root_tag > {}; 

// Base class for functors depending on the algebraic category of the
// Polynomial type
template< class Coefficient_, class PolynomialAlgebraicCategory >
class Polynomial_traits_d_base_polynomial_algebraic_category {        
    public:
        typedef Null_functor    Univariate_content;
        typedef Null_functor    Square_free_factorization;
};

// Specializations
template< class Coefficient_ >
class Polynomial_traits_d_base_polynomial_algebraic_category<
            Polynomial< Coefficient_ >, Integral_domain_without_division_tag >
    : public Polynomial_traits_d_base_polynomial_algebraic_category<
                Polynomial< Coefficient_ >, Null_tag > {};

template< class Coefficient_ >
class Polynomial_traits_d_base_polynomial_algebraic_category<
            Polynomial< Coefficient_ >, Integral_domain_tag >
    : public Polynomial_traits_d_base_polynomial_algebraic_category<
                Polynomial< Coefficient_ >, Integral_domain_without_division_tag > {};

template< class Coefficient_ >
class Polynomial_traits_d_base_polynomial_algebraic_category<
            Polynomial< Coefficient_ >, Unique_factorization_domain_tag >
    : public Polynomial_traits_d_base_polynomial_algebraic_category<
                Polynomial< Coefficient_ >, Integral_domain_tag > {
    CGAL_POLYNOMIAL_TRAITS_D_BASE_TYPEDEFS

    public:
    
    //       Univariate_content
    struct Univariate_content
        : public Unary_function< Polynomial_d , Coefficient>{
        Coefficient operator()(const Polynomial_d& p) const {
            return p.content();
        }
        Coefficient operator()(Polynomial_d p, int i) const {
            return typename PT::Swap()(p,i,PT::d-1).content();
        }
    };
    
    //       Square_free_factorization;
    struct Square_free_factorization{
        typedef int result_type;
        
        template < class OutputIterator1, class OutputIterator2 >
        int operator()(
                const Polynomial_d& p, 
                OutputIterator1 fit, 
                OutputIterator2 mit) const {
            return square_free_factorization( p, fit, mit );
        }
        
        template< class OutputIterator1, class OutputIterator2 >
        int operator()( const Polynomial_d& p, OutputIterator1 fit,
                        OutputIterator2 mit, Innermost_coefficient& a ) {
            if( p == Polynomial_d(0) ) {
                a = Innermost_coefficient(0);
                return 0;
            }
            a = CGAL::unit_part( typename Polynomial_traits_d< Polynomial_d >::Innermost_leading_coefficient()( p ) ) * 
                typename Polynomial_traits_d< Polynomial_d >::Multivariate_content()( p );            
            return square_free_factorization( p/Polynomial_d(a), fit, mit );
        }
          
    };
        
};

template< class Coefficient_ >
class Polynomial_traits_d_base_polynomial_algebraic_category<
            Polynomial< Coefficient_ >, Euclidean_ring_tag >
    : public Polynomial_traits_d_base_polynomial_algebraic_category<
                Polynomial< Coefficient_ >, Unique_factorization_domain_tag > {};


// Polynomial_traits_d_base class connecting the two base classes which depend
//  on the algebraic category of the innermost coefficient type and the poly-
//  nomial type.

// First the general base class for the innermost coefficient
template< class InnermostCoefficient, 
          class ICoeffAlgebraicCategory, class PolynomialAlgebraicCategory >
class Polynomial_traits_d_base {
    typedef InnermostCoefficient ICoeff;
  public:
    static const int d = 0;
    
    typedef ICoeff Polynomial_d;
    typedef ICoeff Coefficient;
    typedef ICoeff Innermost_coefficient;
    
    struct Degree 
        : public Unary_function< ICoeff , int > {
        int operator()(const ICoeff& c) const { return 0; }
    };
    struct Total_degree 
        : public Unary_function< ICoeff , int > {
        int operator()(const ICoeff& c) const { return 0; }
    };
    
    typedef Null_functor  Construct_polynomial;
    typedef Null_functor  Get_coefficient;
    typedef Null_functor  Degree;
    typedef Null_functor  Total_degree;
    typedef Null_functor  Leading_coefficient;
    typedef Null_functor  Univariate_content;
    typedef Null_functor  Multivariate_content;
    typedef Null_functor  Shift;
    typedef Null_functor  Negate;
    typedef Null_functor  Invert;
    typedef Null_functor  Translate;
    typedef Null_functor  Translate_homogeneous;
    typedef Null_functor  Scale_homogeneous;
    typedef Null_functor  Derivative;
    
    struct Is_square_free
        : public Unary_function< ICoeff, bool > {
        bool operator()( const ICoeff& ) const {
            return true;
        }
    };
    
    struct Make_square_free 
        : public Unary_function< ICoeff, ICoeff>{
        ICoeff operator()( const ICoeff& x ) const {
            if (CGAL::is_zero(x)) return x ;
            else  return ICoeff(1);
        }
    };
        
    typedef Null_functor  Square_free_factorization;
    typedef Null_functor  Pseudo_division;
    typedef Null_functor  Pseudo_division_remainder;
    typedef Null_functor  Pseudo_division_quotient;

    struct Gcd_up_to_constant_factor 
        : public Binary_function< ICoeff, ICoeff, ICoeff >{
        ICoeff operator()(const ICoeff& x, const ICoeff& y) const {
            if (CGAL::is_zero(x) && CGAL::is_zero(y)) 
                return ICoeff(0);
            else
                return ICoeff(1);
        }
    };
    typedef Null_functor  Integral_division_up_to_constant_factor;
    struct Univariate_content_up_to_constant_factor
        : public Unary_function< ICoeff, ICoeff >{
        ICoeff operator()(const ICoeff& ) const {
            return ICoeff(1);
        }
    };

    typedef Null_functor  Square_free_factorization_up_to_constant_factor;
    typedef Null_functor  Resultant;
    typedef Null_functor  Canonicalize;
    typedef Null_functor  Evaluate_homogeneous;
    
    struct Innermost_leading_coefficient
        :public Unary_function <ICoeff, ICoeff>{
        ICoeff operator()(const ICoeff& x){return x;}
    };
    struct Degree_vector{
        typedef Exponent_vector         result_type;
        typedef Coefficient             argument_type;
        
        // returns the exponent vector of inner_most_lcoeff. 
        result_type operator()(const Coefficient& constant){
            return Exponent_vector();
        }
    };
    
    struct Get_innermost_coefficient 
        : public Binary_function< ICoeff, Polynomial_d, Exponent_vector > {
        
        ICoeff operator()( const Polynomial_d& p, Exponent_vector ev ) {
            CGAL_precondition( ev.empty() );
            return p;
        }
        
    };
    
    struct Evaluate {        
        template< class Input_iterator >
        ICoeff operator()( const Polynomial_d& p, Input_iterator, Input_iterator ) {
            //std::cerr << p << std::endl;
            return p;
        } 
    };
    
};

// Evaluate_homogeneous_func for recursive homogeneous evaluation of a
//  polynomial, used by Polynomial_traits_d_base for polynomials.
    template< class Polynomial, int d = CGAL::Polynomial_traits_d< Polynomial>::d >
    struct Evaluate_homogeneous_func;

    template< class Polynomial >
    struct Evaluate_homogeneous_func< Polynomial, 1 > {
        typedef typename CGAL::Polynomial_traits_d< Polynomial > PT;
        typedef typename PT::Coefficient Coefficient;
        typedef typename PT::Innermost_coefficient ICoeff;
        typedef typename CGAL::Polynomial_traits_d< Coefficient > PTC;
        
        template< class Input_iterator >
        ICoeff operator()( const Polynomial& p,
                           Input_iterator begin,
                           Input_iterator end,
                           int total_degree,
                           const ICoeff& v ) const {
            --end;
            CGAL_precondition( begin == end );
/*            std::cerr << (*end) << ", " << v << ", " << total_degree << std::endl;
            std::cerr << p << std::endl;*/
            return p.evaluate_homogeneous( (*end), v, total_degree );
        }
        
    }; 

    template< class Polynomial, int d >
    struct Evaluate_homogeneous_func {
        typedef typename CGAL::Polynomial_traits_d< Polynomial > PT;
        typedef typename PT::Coefficient Coefficient;
        typedef typename PT::Innermost_coefficient ICoeff;
        typedef typename CGAL::Polynomial_traits_d< Coefficient > PTC;
        
        template< class Input_iterator >
        ICoeff operator()( const Polynomial& p,
                           Input_iterator begin,
                           Input_iterator end,
                           int total_degree,
                           const ICoeff& v ) const {
            CGAL_precondition( begin != end );
            //typename PT::Evaluate evaluate;
            typename PT::Degree degree;            
            Evaluate_homogeneous_func< Coefficient > eval_hom;
            --end;
            
            std::vector< ICoeff > cv;
            
            for( int i = 0; i <= degree(p); ++i ) {
                cv.push_back( eval_hom( p[i], begin, end, total_degree - i, v ) );                     
            }
             
            return (CGAL::Polynomial< ICoeff >( cv.begin(), cv.end() )).evaluate((*end));
        }
    };

// Now the version for the polynomials with all functors provided by all polynomials
template< class Coefficient_,
          class ICoeffAlgebraicCategory, class PolynomialAlgebraicCategory >
class Polynomial_traits_d_base< Polynomial< Coefficient_ >,
          ICoeffAlgebraicCategory, PolynomialAlgebraicCategory >
    : public Polynomial_traits_d_base_icoeff_algebraic_category< 
        Polynomial< Coefficient_ >, ICoeffAlgebraicCategory >,
      public Polynomial_traits_d_base_polynomial_algebraic_category<
        Polynomial< Coefficient_ >, PolynomialAlgebraicCategory > {

    CGAL_POLYNOMIAL_TRAITS_D_BASE_TYPEDEFS

public:

    //
    // Functors as defined in the reference manual (with sometimes slightly
    //   extended functionality)
    //


    //       Construct_polynomial;
    struct Construct_polynomial {
        
        typedef Polynomial_d  result_type;
        
        Polynomial_d operator()()  const {
            return Polynomial_d(0);
        }
        
        template <class T>
        Polynomial_d operator()( T a ) const {
            return Polynomial_d(a);
        }
        
        //! construct the constant polynomial a0
        Polynomial_d operator() (const Coefficient& a0) const
        {return Polynomial_d(a0);}
        
        //! construct the polynomial a0 + a1*x
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1) const
        {return Polynomial_d(a0,a1);}
        
        //! construct the polynomial a0 + a1*x + a2*x^2
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2) const
        {return Polynomial_d(a0,a1,a2);}
        
        //! construct the polynomial a0 + a1*x + ... + a3*x^3
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3) const
        {return Polynomial_d(a0,a1,a2,a3);}
        
        //! construct the polynomial a0 + a1*x + ... + a4*x^4
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3,
                const Coefficient& a4) const
        {return Polynomial_d(a0,a1,a2,a3,a4);}
        
        //! construct the polynomial a0 + a1*x + ... + a5*x^5
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3,
                const Coefficient& a4, const Coefficient& a5) const
        {return Polynomial_d(a0,a1,a2,a3,a4,a5);}
        
        //! construct the polynomial a0 + a1*x + ... + a6*x^6
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3,
                const Coefficient& a4, const Coefficient& a5, 
                const Coefficient& a6) const
        {return Polynomial_d(a0,a1,a2,a3,a4,a5,a6);}
        
        //! construct the polynomial a0 + a1*x + ... + a7*x^7
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3,
                const Coefficient& a4, const Coefficient& a5, 
                const Coefficient& a6, const Coefficient& a7) const
        {return Polynomial_d(a0,a1,a2,a3,a4,a5,a6,a7);}
        
        //! construct the polynomial a0 + a1*x + ... + a8*x^8
        Polynomial_d operator() (
                const Coefficient& a0, const Coefficient& a1,
                const Coefficient& a2, const Coefficient& a3,
                const Coefficient& a4, const Coefficient& a5, 
                const Coefficient& a6, const Coefficient& a7,
                const Coefficient& a8) const 
        {return Polynomial_d(a0,a1,a2,a3,a4,a5,a6,a7,a8);}
        
        template< class Input_iterator >
        inline
        Polynomial_d construct( 
                Input_iterator begin, 
                Input_iterator end , 
                Tag_true) const {
            return Polynomial_d(begin,end);
        }
        
        template< class Input_iterator >
        inline
        Polynomial_d construct( 
                Input_iterator begin, 
                Input_iterator end , 
                Tag_false) const {
            std::sort(begin,end); 
            return Create_polynomial_from_monom_rep< Coefficient >()
                ( begin, end ); 
        }
        

        template< class Input_iterator >
        Polynomial_d 
        operator()( Input_iterator begin, Input_iterator end ) const {
            if(begin == end ) return Polynomial_d(0);
            typedef typename Input_iterator::value_type value_type;
            typedef Boolean_tag<boost::is_same<value_type,Coefficient>::value> 
                Is_coeff;
            return construct(begin,end,Is_coeff());
        }
    
    private:
    public: 

        template< class T >
        class Create_polynomial_from_monom_rep {
        public:
            template <class Monom_rep_iterator>
            Polynomial_d operator()( 
                    Monom_rep_iterator begin,
                    Monom_rep_iterator end) const {
                
                std::vector< Innermost_coefficient > coefficients;
                for(Monom_rep_iterator it = begin; it != end; it++){
                    while( it->first[0] > (int) coefficients.size() ){
                        coefficients.push_back(Innermost_coefficient(0));
                    }
                    coefficients.push_back(it->second);
                }
                return Polynomial_d(coefficients.begin(),coefficients.end());
            }
        };
        template< class T >
        class Create_polynomial_from_monom_rep< Polynomial < T > > {
        public:
            template <class Monom_rep_iterator>
            Polynomial_d operator()( 
                    Monom_rep_iterator begin,
                    Monom_rep_iterator end) const {
                //std::cout << " ------\n "  << std::endl;
                
                typedef Polynomial_traits_d<Coefficient> PT;
                typename PT::Construct_polynomial construct;
                
                BOOST_STATIC_ASSERT(PT::d != 0); // Coefficient is a Polynomial
                std::vector<Coefficient> coefficients;
                
                Monom_rep_iterator it = begin; 
                while(it != end){
                    int current_exp = it->first[PT::d];
                    //std::cout <<"current_exp: " <<  current_exp << std::endl;
                    // fill up with zeros until current exp is reached
                    while( (int) coefficients.size() < current_exp){
                        coefficients.push_back(Coefficient(0));
                        //std::cout <<" insert "<< std::endl;
                    }
                    // collect all coeffs for this exp
                    Monom_rep monoms; 
                    while(  it != end && it->first[PT::d] == current_exp ){
                        Exponent_vector ev = it->first;
                        ev.pop_back();
                        monoms.push_back( Exponents_coeff_pair(ev,it->second));
                        it++;
                    }                    
                    coefficients.push_back(
                                construct(monoms.begin(), monoms.end()));
                }
                //std::cout << " ------\n "  << std::endl;
                return Polynomial_d(coefficients.begin(),coefficients.end());
            }
        };
    };

    //       Get_coefficient;
    struct Get_coefficient 
        : public Binary_function<Polynomial_d, int, Coefficient > {
        
        Coefficient operator()( const Polynomial_d& p, int i) const {
            CGAL_precondition( i >= 0 );
            typename PT::Degree degree;
            if( i >  degree(p) )
                return Coefficient(0);
            return p[i];
        }
            
    };
    
    //       Get_innermost_coefficient;
    struct Get_innermost_coefficient
        : public Binary_function< Polynomial_d, Exponent_vector, Innermost_coefficient > {
        
        Innermost_coefficient operator()( const Polynomial_d& p, Exponent_vector ev ) const {
            CGAL_precondition( !ev.empty() );
            typename PTC::Get_innermost_coefficient gic;
            typename PT::Get_coefficient gc;
            int exponent = ev.back();
            ev.pop_back();
            return gic( gc( p, exponent ), ev ); 
        }; 
    
    };
     
    //       Swap;
    // Swap variable x_i with x_j
    struct Swap {
        typedef Polynomial_d        result_type;  
        typedef Polynomial_d        first_argument_type;
        typedef int                 second_argument_type;
        typedef int                 third_argument_type;
        typedef Arity_tag< 3 >         Arity;
        
        // We use our own Strict Weak Ordering predicate in order to avoid
        // problems when calling sort for a Exponents_coeff_pair where the
        // coeff type has no comparison operators available.
      private:
        struct Compare_exponents_coeff_pair 
            : public Binary_function< std::pair< Exponent_vector, Innermost_coefficient >,
                                      std::pair< Exponent_vector, Innermost_coefficient >,
                                      bool > {
            bool operator()( const std::pair< Exponent_vector, Innermost_coefficient >& p1,
                             const std::pair< Exponent_vector, Innermost_coefficient >& p2 ) const {
                // TODO: Precondition leads to an error within test_translate in Polynomial_traits_d test
                //CGAL_precondition( p1.first != p2.first );
                return p1.first < p2.first;
            }            
        }; 
      
      public:
        
        Polynomial_d operator()(const Polynomial_d& p, int i, int j ) const {
            //std::cout << i <<" " << j << " : " ; 
            CGAL_precondition(0 <= i && i < d);
            CGAL_precondition(0 <= j && j < d);
            typedef std::pair< Exponent_vector, Innermost_coefficient >
                Exponents_coeff_pair;
            typedef std::vector< Exponents_coeff_pair > Monom_rep; 
            Get_monom_representation gmr;
            typename Construct_polynomial::template Create_polynomial_from_monom_rep< Coefficient > construct;
            Monom_rep mon_rep;
            gmr( p, std::back_inserter( mon_rep ) );
            for( typename Monom_rep::iterator it = mon_rep.begin(); 
                 it != mon_rep.end();
                 ++it ) {
                std::swap(it->first[i],it->first[j]);
                // it->first.swap( i, j );
            }
            std::sort( mon_rep.begin(), mon_rep.end(), Compare_exponents_coeff_pair() );
            return construct( mon_rep.begin(), mon_rep.end() );
        }
    };

    //       Move;    
    // move variable x_i to position of x_j
    // order of other variables remains 
    // default j = d makes x_i the othermost variable
    struct Move {
        typedef Polynomial_d        result_type;  
        typedef Polynomial_d        first_argument_type;
        typedef int                 second_argument_type;
        typedef int                 third_argument_type;
        typedef Arity_tag< 3 >         Arity;
        
        Polynomial_d operator()(const Polynomial_d& p, int i, int j = (d-1) ) const {
            //std::cout << x <<" " << y << " : " ; 
            CGAL_precondition(0 <= i && i < d);
            CGAL_precondition(0 <= j && j < d);
            typedef std::pair< Exponent_vector, Innermost_coefficient >
                Exponents_coeff_pair;
            typedef std::vector< Exponents_coeff_pair > Monom_rep; 
            Get_monom_representation gmr;
            Construct_polynomial construct;
            Monom_rep mon_rep;
            gmr( p, std::back_inserter( mon_rep ) );
            for( typename Monom_rep::iterator it = mon_rep.begin(); 
                 it != mon_rep.end();
                 ++it ) {
                // this is as good as std::rotate since it uses swap also
                if (i < j) 
                    for( int k = i; k < j; k++ )
                        std::swap(it->first[k],it->first[k+1]);
                else
                    for( int k = i; k > j; k-- )
                        std::swap(it->first[k],it->first[k-1]);
                
            }
            std::sort( mon_rep.begin(), mon_rep.end() );
            return construct( mon_rep.begin(), mon_rep.end() );
        }
    };

    //       Degree;    
    struct Degree : public Unary_function< Polynomial_d , int  >{
        int operator()(const Polynomial_d& p, int i = (d-1)) const {      
            if (i == (d-1)) return p.degree();                           
            else return Swap()(p,i,d-1).degree();
        }     
    };

    //       Total_degree;
    struct Total_degree : public Unary_function< Polynomial_d , int >{
        int operator()(const Polynomial_d& p) const {
            typedef Polynomial_traits_d<Coefficient> COEFF_POLY_TRAITS;
            typename COEFF_POLY_TRAITS::Total_degree total_degree;
            Degree degree;
            CGAL_precondition( degree(p) >= 0);

            int result = 0;
            for(int i = 0; i <= degree(p) ; i++){
                if( ! CGAL::is_zero( p[i]) )
                    result = std::max(result , total_degree(p[i]) + i );
            } 
            return result;
        }
    };

    //       Leading_coefficient;
    struct Leading_coefficient 
        : public Unary_function< Polynomial_d , Coefficient>{
        Coefficient operator()(const Polynomial_d& p) const {
            return p.lcoeff();
        }
        Coefficient operator()(Polynomial_d p, int i) const {
            return Swap()(p,i,PT::d-1).lcoeff();
        }
    };
    
    //       Innermost_leading_coefficient;
    struct Innermost_leading_coefficient 
        : public Unary_function< Polynomial_d , Innermost_coefficient>{
        Innermost_coefficient 
        operator()(const Polynomial_d& p) const {
            typename PTC::Innermost_leading_coefficient ilcoeff;
            typename PT::Leading_coefficient lcoeff;
            return ilcoeff(lcoeff(p));
        }
    };

    //       Canonicalize;
    struct Canonicalize
        : public Unary_function<Polynomial_d, Polynomial_d>{
        Polynomial_d
        operator()( const Polynomial_d& p ) const {
            return CGAL::POLYNOMIAL::canonicalize_polynomial(p);
        }  
     };

    //       Derivative;
     struct Derivative 
         : public Unary_function<Polynomial_d, Polynomial_d>{
        Polynomial_d
        operator()(Polynomial_d p, int i = (d-1)) const {
            if (i == (d-1) ){
                p.diff();
            }else{
                Swap swap;
                p = swap(p,i,d-1);
                p.diff();
                p = swap(p,i,d-1);
            }
            return p;
         }
    };

    //       Evaluate;
    struct Evaluate
        :public Unary_function<Polynomial_d,Innermost_coefficient>{
        Coefficient
        operator()(const Polynomial_d& p, Innermost_coefficient x, int i = (d-1)) 
            const {
            if(i == (d-1) )
                return p.evaluate(x);
            else{
                return Move()(p,i).evaluate(x);
            }
        }
        
        template< class Input_iterator >
        Innermost_coefficient operator()( const Polynomial_d& p, Input_iterator begin, Input_iterator end ) const {
            CGAL_precondition( begin != end );
            
            typename PT::Evaluate evaluatePoly;
            typename PTC::Evaluate evaluateCoeff;
            --end;
            return evaluateCoeff( evaluatePoly( p, (*end) ), begin, end );            
        }  
    };
    
    //       Evaluate_homogeneous;
    struct Evaluate_homogeneous{
        typedef Coefficient           result_type;  
        typedef Polynomial_d          first_argument_type;
        typedef Innermost_coefficient second_argument_type;
        typedef Innermost_coefficient third_argument_type;
        typedef Arity_tag< 3 >         Arity;
        
        Coefficient
        operator()(
                const Polynomial_d& p, 
                Innermost_coefficient a, 
                Innermost_coefficient b,
                int i = (PT::d-1) ) const {
            if (i == (d-1) )
                return p.evaluate_homogeneous(a,b);
            else
                return Move()(p,i,PT::d-1).evaluate_homogeneous(a,b);
        }
        
        template< class Input_iterator >
        Innermost_coefficient operator()( const Polynomial_d & p,
                                          Input_iterator begin,
                                          Input_iterator end ) const {
            typename PT::Total_degree total_degree;
            typename PT::Evaluate_homogeneous eval_hom;
            return eval_hom( p, begin, end, total_degree(p) );
        }
        
        template< class Input_iterator >
        Innermost_coefficient operator()( const Polynomial_d& p,
                                          Input_iterator begin,
                                          Input_iterator end,
                                          int total_degree ) const {
            --end;
            Evaluate_homogeneous_func< Polynomial_d > eval_hom;
            return eval_hom( p, begin, end, total_degree, (*end) );
        }
        
    };
    
    //       Is_zero_at;
    struct Is_zero_at {
        typedef bool    result_type;
        
        template< class Input_iterator >
        bool operator()( const Polynomial_d& p, Input_iterator begin, Input_iterator end ) const {
            typename PT::Evaluate evaluate;            
            return( CGAL::is_zero( evaluate( p, begin, end ) ) );
        } 
    };
    
    //       Is_zero_at_homogeneous;
    struct Is_zero_at_homogeneous {
        typedef bool    result_type;
        
        template< class Input_iterator >
        bool operator()( const Polynomial_d& p, Input_iterator begin, Input_iterator end ) const {
            typename PT::Evaluate_homogeneous evaluate_homogeneous;
            return( CGAL::is_zero( evaluate_homogeneous( p, begin, end ) ) );
        }
    };

    //       Sign_at;
    struct Sign_at {
        typedef Sign    result_type;
        
        template< class Input_iterator >
        Sign operator()( const Polynomial_d& p, Input_iterator begin, Input_iterator end ) const {
            typename PT::Evaluate evaluate;
            return CGAL::sign( evaluate( p, begin, end ) );
        }
    };
    
    //       Sign_at_homogeneous;
    struct Sign_at_homogeneous {
        typedef Sign    result_type;
        
        template< class Input_iterator >
        Sign operator()( const Polynomial_d& p, Input_iterator begin, Input_iterator end ) const {
            typename PT::Evaluate_homogeneous evaluate_homogeneous;
            return CGAL::sign( evaluate_homogeneous( p, begin, end ) );
        }
    };
    
    //      Compare;
    struct Compare
        : public Binary_function< Comparison_result, Polynomial_d, Polynomial_d > {
        
        Comparison_result operator()( const Polynomial_d& p1, const Polynomial_d& p2 ) const {
            return p1.compare( p2 );
        }
    };

    //
    // This is going to be in PolynomialToolBox
    // 
    struct Coefficient_begin                                                  
        : public Unary_function< Polynomial_d, Coefficient_iterator > {       
        Coefficient_iterator                                                  
        operator () (const Polynomial_d& p) { return p.begin(); }             
    };                                                                        
    struct Coefficient_end                                                    
        : public Unary_function< Polynomial_d, Coefficient_iterator > {       
        Coefficient_iterator                                                  
        operator () (const Polynomial_d& p) { return p.end(); }               
    };                                                                        
                                                                              
    struct Innermost_coefficient_begin                                        
        : public Unary_function< Polynomial_d, Innermost_coefficient_iterator > {  
        Innermost_coefficient_iterator                                             
        operator () (const Polynomial_d& p) {                                      
            return typename Coefficient_flattening::Flatten()(p.end(),p.begin());  
        }                                                                          
    };                                                                             
                                                                                   
    struct Innermost_coefficient_end                                               
      : public Unary_function< Polynomial_d, Innermost_coefficient_iterator > {    
        Innermost_coefficient_iterator                                             
        operator () (const Polynomial_d& p) {                                      
            return typename Coefficient_flattening::Flatten()(p.end(),p.end());    
        }                                                                          
    };                                                                             

    //       Is_square_free;
    struct Is_square_free 
        : public Unary_function< Polynomial_d, bool >{
        bool operator()( const Polynomial_d& p ) const {
            if( !POLYNOMIAL::may_have_multiple_factor( p ) )
                return true;
            
            Univariate_content_up_to_constant_factor ucontent_utcf;
            Integral_division_up_to_constant_factor  idiv_utcf;
            Derivative diff;
            
            Coefficient content = ucontent_utcf( p );
            typename PTC::Is_square_free isf;
            
            if( !isf( content ) )
                return false;
            
            Polynomial_d regular_part = idiv_utcf( p, Polynomial_d( content ) ); 
            
            Polynomial_d g = gcd_utcf(regular_part,diff(regular_part));
            return ( g.degree() == 0 );
        }
    };

                   
    //       Make_square_free;
    struct Make_square_free 
        : public Unary_function< Polynomial_d, Polynomial_d >{
        Polynomial_d
        operator()(const Polynomial_d& p) const {
            if (CGAL::is_zero(p)) return p;
            Univariate_content_up_to_constant_factor ucontent_utcf;
            Integral_division_up_to_constant_factor  idiv_utcf;
            Derivative diff;
            typename PTC::Make_square_free msf;
            
            Coefficient content = ucontent_utcf(p);
            Polynomial_d result = Polynomial_d(msf(content));
            
            Polynomial_d regular_part = idiv_utcf(p,Polynomial_d(content));
            Polynomial_d g = gcd_utcf(regular_part,diff(regular_part));
            
            
            result *= idiv_utcf(regular_part,g);
            return Canonicalize()(result);
                    
        }
    };

    //       Pseudo_division;
    struct Pseudo_division {
        typedef Polynomial_d        result_type;  
        void
        operator()(
                const Polynomial_d& f, const Polynomial_d& g,
                Polynomial_d& q, Polynomial_d& r, Coefficient& D) const {
            Polynomial_d::pseudo_division(f,g,q,r,D);
        }
    };
    
    //       Pseudo_division_quotient;
    struct Pseudo_division_quotient
        :public Binary_function<Polynomial_d, Polynomial_d, Polynomial_d> {
        
        Polynomial_d
        operator()(const Polynomial_d& f, const Polynomial_d& g) const {
            Polynomial_d q,r;
            Coefficient D;
            Polynomial_d::pseudo_division(f,g,q,r,D);
            return q;
        }
    };

    //       Pseudo_division_remainder;
    struct Pseudo_division_remainder
        :public Binary_function<Polynomial_d, Polynomial_d, Polynomial_d> {
        
        Polynomial_d
        operator()(const Polynomial_d& f, const Polynomial_d& g) const {
            Polynomial_d q,r;
            Coefficient D;
            Polynomial_d::pseudo_division(f,g,q,r,D);
            return r;
        }
    };
    
    //       Gcd_up_to_constant_factor;
    struct Gcd_up_to_constant_factor
        :public Binary_function<Polynomial_d, Polynomial_d, Polynomial_d> {
        Polynomial_d
        operator()(const Polynomial_d& p, const Polynomial_d& q) const {
            if (CGAL::is_zero(p) && CGAL::is_zero(q)) 
                return Polynomial_d(0);
            return gcd_utcf(p,q);
        }
    };
    
    //       Integral_division_up_to_constant_factor;
    struct Integral_division_up_to_constant_factor
        :public Binary_function<Polynomial_d, Polynomial_d, Polynomial_d> {
        Polynomial_d
        operator()(const Polynomial_d& p, const Polynomial_d& q) const {
            typedef Innermost_coefficient IC;

            typename PT::Construct_polynomial construct;
            typename PT::Innermost_leading_coefficient ilcoeff;
            typename PT::Innermost_coefficient_begin begin;
            typename PT::Innermost_coefficient_end end;
            typedef Algebraic_extension_traits<Innermost_coefficient> AET;
            typename AET::Denominator_for_algebraic_integers dfai;
            typename AET::Normalization_factor nfac;

            
            IC ilcoeff_q = ilcoeff(q);
            // this factor is needed in case IC is an Algebraic extension
            IC dfai_q = dfai(begin(q), end(q));
            // make dfai_q a 'scalar'
            ilcoeff_q *= dfai_q * nfac(dfai_q);
            
            Polynomial_d result = (p * construct(ilcoeff_q)) / q;

            return Canonicalize()(result);
        }
    };
    
    //       Univariate_content_up_to_constant_factor;
    struct Univariate_content_up_to_constant_factor
        :public Unary_function<Polynomial_d, Coefficient> {
        Coefficient
        operator()(const Polynomial_d& p) const {
            typename PTC::Gcd_up_to_constant_factor gcd_utcf;
            
            if(CGAL::is_zero(p)) return Coefficient(0);
            if(PT::d == 1) return Coefficient(1);

            Coefficient result(0);
            for(typename Polynomial_d::const_iterator it = p.begin();
                it != p.end();
                it++){
                result = gcd_utcf(*it,result);
            }
            return result;

        }
    };

    //       Square_free_factorization_up_to_constant_factor;
    struct Square_free_factorization_up_to_constant_factor {
        typedef int result_type;
    private:
        typedef Coefficient Coeff;
        typedef Innermost_coefficient ICoeff;
        
        // rsqff_utcf computes the sqff recursively for Coeff  
        // end of recursion: ICoeff
        template < class OutputIterator1, class OutputIterator2 >
        int rsqff_utcf  (ICoeff c, 
                OutputIterator1 factors, 
                OutputIterator2 mults) const{
            return 0;
        }        
        template < class OutputIterator1, class OutputIterator2 >
        int rsqff_utcf (
                typename First_if_different<Coeff,ICoeff>::Type c,
                OutputIterator1 fit, 
                OutputIterator2 mit) const {
            typename PTC::Square_free_factorization_up_to_constant_factor sqff;
            std::vector<Coefficient> factors;
            int n = sqff(c, std::back_inserter(factors), mit);
            for(int i = 0; i < (int)factors.size(); i++){
                *fit++=Polynomial_d(factors[i]);
            }
            return n; 
        }
    public:
        template < class OutputIterator1, class OutputIterator2 >
        int operator()(
                Polynomial_d p, 
                OutputIterator1 fit, 
                OutputIterator2 mit) const {
            
            if (CGAL::is_zero(p)) return 0;

            Univariate_content_up_to_constant_factor ucontent_utcf;
            Integral_division_up_to_constant_factor idiv_utcf;
            Coefficient c = ucontent_utcf(p);
            p = idiv_utcf( p , Polynomial_d(c));
            int n = square_free_factorization_utcf(p,fit,mit);
            if (Total_degree()(c) > 0) 
                    return rsqff_utcf(c,fit,mit)+n;
                else 
                    return n;
        }
    };

    //       Shift;
    struct Shift
        : public Unary_function< Polynomial_d, Polynomial_d >{
        Polynomial_d 
        operator()(const Polynomial_d& p, int e, int i = PT::d) const {
            Construct_polynomial construct; 
            Get_monom_representation gmr; 
            Monom_rep monom_rep;
            gmr(p,std::back_inserter(monom_rep));
            for(typename Monom_rep::iterator it = monom_rep.begin(); 
                it != monom_rep.end();
                it++){
                it->first[i-1]+=e;
            }
            return construct(monom_rep.begin(), monom_rep.end());
        }
    };
    
    //       Negate;
    struct Negate
        : public Unary_function< Polynomial_d, Polynomial_d >{
        
        Polynomial_d operator()(const Polynomial_d& p, int i = (d-1)) const {
            Construct_polynomial construct; 
            Get_monom_representation gmr; 
            Monom_rep monom_rep;
            gmr(p,std::back_inserter(monom_rep));
            for(typename Monom_rep::iterator it = monom_rep.begin(); 
                it != monom_rep.end();
                it++){
                if (it->first[i] % 2 != 0) 
                    it->second = - it->second; 
            }
            return construct(monom_rep.begin(), monom_rep.end());
        }
    };

    //       Invert;
    struct Invert
        : public Unary_function< Polynomial_d , Polynomial_d >{
        Polynomial_d operator()(Polynomial_d p, int i = (PT::d-1)) const {
            if (i == (d-1)){
                p.reversal(); 
            }else{
                p =  Swap()(p,i,PT::d-1);
                p.reversal();
                p = Swap()(p,i,PT::d-1);   
            }
            return p ;
        }
    };

    //       Translate;
    struct Translate
        : public Binary_function< Polynomial_d , Polynomial_d, 
                                  Innermost_coefficient >{
        Polynomial_d
        operator()(
                Polynomial_d p, 
                const Innermost_coefficient& c, 
                int i = (d-1)) 
            const {
            if (i == (d-1) ){
                p.translate(Coefficient(c)); 
            }else{
                Swap swap;
                p = swap(p,i,d-1);
                p.translate(Coefficient(c));
                p = swap(p,i,d-1); 
            }
            return p;
        }
    };

    //       Translate_homogeneous;
    struct Translate_homogeneous{
        typedef Polynomial_d result_type;
        typedef Polynomial_d first_argument_type;
        typedef Innermost_coefficient second_argument_type;
        typedef Innermost_coefficient third_argument_type;
        
        Polynomial_d
        operator()(Polynomial_d p, 
                const Innermost_coefficient& a, 
                const Innermost_coefficient& b,
                int i = (d-1) ) const {
            if (i == (d-1) ){
                p.translate(Coefficient(a), Coefficient(b) );  
            }else{
                Swap swap;
                p = swap(p,i,d-1);
                p.translate(Coefficient(a), Coefficient(b));
                p = swap(p,i,d-1);
            }
            return p;
         }
    };

    //       Scale;
    struct Scale 
        : public Binary_function< Polynomial_d, Innermost_coefficient, Polynomial_d > {
        
        Polynomial_d operator()( Polynomial_d p, const Innermost_coefficient& c,
                                 int i = (PT::d-1) ) {
            typename PT::Scale_homogeneous scale_homogeneous;
            
            return scale_homogeneous( p, c, Innermost_coefficient(1), i );
        }
        
    };
    
    //       Scale_homogeneous;
    struct Scale_homogeneous{
        typedef Polynomial_d result_type;
        typedef Polynomial_d first_argument_type;
        typedef Innermost_coefficient second_argument_type;
        typedef Innermost_coefficient third_argument_type;
        
        Polynomial_d
        operator()(
                Polynomial_d p, 
                const Innermost_coefficient& a, 
                const Innermost_coefficient& b,
                int i = (d-1) ) const {
            CGAL_precondition( ! CGAL::is_zero(b) );
            
            if (i == (d-1) ) p = Swap()(p,i,d-1);
          
            if(CGAL::is_one(b)) 
                p.scale_up(Coefficient(a));
            else 
                if(CGAL::is_one(a)) 
                    p.scale_down(Coefficient(b));
                else 
                    p.scale(Coefficient(a), Coefficient(b) );  
          
            if (i == (d-1) ) p = Swap()(p,i,d-1);
          
            return p;
         }
    };

    //       Resultant;
    struct Resultant
        : public Binary_function<Polynomial_d, Polynomial_d, Coefficient>{
        
        Coefficient
        operator()(
                const Polynomial_d& p, 
                const Polynomial_d& q,
                int i = (d-1) ) const {
            if(i == (d-1) )
                return resultant(p,q);
            else
                return resultant(Move()(p,i),Move()(q,i));
        }  
     };
    
    //
    // Functors not mentioned in the reference manual
    //

    struct Get_monom_representation {      
        typedef std::pair< Exponent_vector, Innermost_coefficient >
        Exponents_coeff_pair;
        typedef std::vector< Exponents_coeff_pair > Monom_rep; 
        
        template <class OutputIterator>
        void operator()( const Polynomial_d& p, OutputIterator oit ) const {
            typedef Boolean_tag< d == 1 > Is_univariat;
            create_monom_representation( p, oit , Is_univariat());
        }
      
      private:
        
        template <class OutputIterator>
        void
        create_monom_representation 
        ( const Polynomial_d& p, OutputIterator oit, Tag_true ) const{
            for( int exponent = 0; exponent <= p.degree(); ++exponent ) {
                // std::cout << "p[exponent]: "<<p[exponent];
                if ( p[exponent] != Coefficient(0) ){
                    Exponent_vector exp_vec;
                    exp_vec.push_back( exponent );
                    *oit = Exponents_coeff_pair( exp_vec, p[exponent] ); 
                }
            } 
        } 
        template <class OutputIterator>
        void 
        create_monom_representation 
        ( const Polynomial_d& p, OutputIterator oit, Tag_false ) const { 
            for( int exponent = 0; exponent <= p.degree(); ++exponent ) {
                Monom_rep monom_rep;
                typedef Polynomial_traits_d< Coefficient > PT;
                typename PT::Get_monom_representation gmr;
                gmr( p[exponent], std::back_inserter( monom_rep ) );
                for( typename Monom_rep::iterator it = monom_rep.begin();
                     it != monom_rep.end(); ++it ) {
                    it->first.push_back( exponent );
                }
                copy( monom_rep.begin(), monom_rep.end(), oit );               
            }
        }
    };

    // returns the Exponten_vector of the innermost leading coefficient 
    struct Degree_vector{
        typedef Exponent_vector           result_type;
        typedef Polynomial_d              argument_type;
        
        // returns the exponent vector of inner_most_lcoeff. 
        result_type operator()(const Polynomial_d& polynomial){
            
            typename PTC::Degree_vector degree_vector;
            
            Exponent_vector result = degree_vector(polynomial.lcoeff());
            result.insert(result.begin(),polynomial.degree());
            return result;
        }
    };        
};

} // namespace POLYNOMIAL

// Definition of Polynomial_traits_d
//
// In order to determine the algebraic category of the innermost coefficient,
//  the Polynomial_traits_d_base class with "Null_tag" is used.

template< class Polynomial >
class Polynomial_traits_d
    : public POLYNOMIAL::Polynomial_traits_d_base< Polynomial,  
    typename Algebraic_structure_traits<
        typename POLYNOMIAL::Innermost_coefficient<Polynomial>::Type >::Algebraic_category,
    typename Algebraic_structure_traits< Polynomial >::Algebraic_category > {

//------------ Rebind ----------- 
private:
template <class T, int d>
struct Gen_polynomial_type{
    typedef CGAL::Polynomial<typename Gen_polynomial_type<T,d-1>::Type> Type;
};
template <class T>
struct Gen_polynomial_type<T,0>{ typedef T Type; };

public:
template <class T, int d>
struct Rebind{
    typedef Polynomial_traits_d<typename Gen_polynomial_type<T,d>::Type> Other;
};
//------------ Rebind ----------- 

};


CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_TRAITS_D_H
