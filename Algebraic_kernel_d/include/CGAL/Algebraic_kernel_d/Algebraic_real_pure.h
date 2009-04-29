// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_rep.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/convert_to_bfi.h>

#include <iterator>
#include <list>
#include <vector>
#include <queue>

CGAL_BEGIN_NAMESPACE 
namespace CGALi {
template <class Coefficient_,  class Rational_, 
          class HandlePolicy , class RepClass  >
class Algebraic_real_pure;
}

template <class COEFF, class RAT, class POLICY, class REPCLASS >
typename Get_arithmetic_kernel<COEFF>::Arithmetic_kernel::Bigfloat_interval
inline
convert_to_bfi(
const CGALi::Algebraic_real_pure< COEFF, RAT, POLICY, REPCLASS >& x);

namespace CGALi { 
/*! \ingroup NiX_Algebraic_real
  \brief  An Algebraic_real_pure \a x is represented by a polynomial and an 
  isolating interval. It is guaranteed that the polynomial is square free. 
  The isolating interval is an open interval and it is guaranteed that the 
  polynomial is non-zero at the endpoints of the interval. 
  
  The algebraic real are reference counted by default.
  The template parameters are:
        - \b Coefficient_: a model of the \c IntegralDomain concept. 
        - \b Rational_: a model of the \c Field concept.
             Must be assignable to \c Field_with_sqrt_.
        - \b HandlePolicy: a model of the \c HandlePolicy concept or the
              \c Handle_policy_in_place class template that selects a specialized
              implementation without reference counting. Has the
              default \c Handle_policy_union.
  
  THIS CLASS IS CONSIDERED AS EXPERIMENTAL !
*/
template <class Coefficient_, 
          class Rational_,            
          class HandlePolicy = ::CGAL::Handle_policy_no_union,
          class RepClass = CGALi::Algebraic_real_rep< Coefficient_, Rational_ > >
class Algebraic_real_pure :
        public
::CGAL::Handle_with_policy< RepClass, HandlePolicy > {    
public :
    typedef Coefficient_                                         Coefficient;
    typedef Rational_                                            Rational;
    typedef typename
        CGAL::Polynomial_type_generator<Coefficient,1>::Type     Polynomial;
  
    typedef RepClass    Int_rep;

    typedef ::CGAL::Handle_with_policy<Int_rep,HandlePolicy>                       Base;
    typedef Algebraic_real_pure<Coefficient,Rational, HandlePolicy, RepClass> Self;
private:
    typedef CGAL::Fraction_traits<Rational> FT_rational;
    typedef typename FT_rational::Numerator_type Integer;
public:
 
    //! creates the algebraic real as \c Rational from \a i.
    explicit Algebraic_real_pure(int i = 0 ) : Base(i) { }
    
    //! creates the algebraic real as \c Rational from \a i.
    explicit Algebraic_real_pure(const Rational& i) : Base(i) { }

    /*! \brief creates the algebraic real as the unique root of \a P
     *   in the open interval <var>]low,high[</var>.
     *   \pre P is square free.
     *   \pre P is not zero at low
     *   \pre P is not zero at high
     */ 
    Algebraic_real_pure(const Polynomial& P, Rational low, Rational high) 
        : Base (Int_rep(P,low,high)) {}
 
    //! returns the polynomial defining \a x
    const Polynomial& polynomial() const { return this->ptr()->polynomial(); }

    /*! \brief returns the degree of the polynomial defining \a x
     *  This is not necessarily the algebraic degree of \a x, since 
     *  the polynomial may be reducible.
     */
    const int degree() const { return this->ptr()->polynomial().degree(); }
    
    //! returns the lower endpoint of the isolating interval
    Rational low() const { return this->ptr()->low(); }
  
    //! returns the upper endpoint of the isolating interval 
    Rational high() const { return this->ptr()->high(); }

    /*! \brief returns the sign of the defining polynomial 
     *  at the lower endpoint of the isolating interval  
     */
    int sign_at_low() const { return this->ptr()->sign_at_low(); }
    
    /*! \brief returns whether a rational representation 
     * (as Rational) is known. 
     */
    bool is_rational() const { return this->ptr()->is_rational(); }

    //! returns a rational representation of \a x 
    /*! \pre: type() == NiX::IS_RATIONAL */
    Rational rational() const { 
        CGAL_precondition(is_rational());
        return this->ptr()->rational(); 
    }

    //! compute a \c double approximation (without guarantees)
    double to_double() const { 
        return (to_interval().first +to_interval().second)/2;
    }

    /*! \brief returns a double Interval approximation, 
     *  it is guaranteed that \a x is contained in the Interval.
     *
     *  This is not the isolating interval.
     */
    std::pair<double, double> to_interval() const { 
        if( this->ptr()->interval_option ) {
            return *(this->ptr()->interval_option);
        } else {
            typedef typename Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Bigfloat_interval BFI;
            long old_precision = get_precision( BFI() );
            set_precision( BFI(), 53 );
            std::pair<double, double> interval = CGAL::to_interval( convert_to_bfi( (*this)));
            this->ptr()->interval_option = boost::optional< std::pair<double, double> >(interval);
            set_precision( BFI(), old_precision );
            return *(this->ptr()->interval_option);            
        }
    }
        
    /*! \brief Refines the isolating interval. */ 
    void refine() const{ this->ptr()->refine(); }  

    /*! \brief Bisects the isolating interval. */ 
    void bisect() const{ this->ptr()->bisect(); }  
  
    /*! \brief Refines the isolating interval until \a m is outside 
     *  the \c closed interval        
     */    
    template < class NTX >
    void strong_refine(const NTX& m) const{ compare(m); } 
        
    /*! \brief compares \a x with respect to \a y  
     *  It returns CGAL::SMALLER if \a x is smaller than <var>y</var>, 
     *  CGAL::EQUAL if they are equal and CGAL::LARGER otherwise. 
     */
    template < class NTX >
    CGAL::Comparison_result
    compare(const NTX& y) const { return intern_compare(y,false); } 
    

    /*! \brief compares \a x with respect to \a y
     *  \pre \a x != \a y 
     */  
    template < class NTX > 
    CGAL::Comparison_result
    compare_distinct(const NTX& y) const { return intern_compare(y,true); }
   
    template <class Algebraic_real_iterator>
    static void conjugate(Algebraic_real_iterator begin,
                          Algebraic_real_iterator end){
        if(begin == end) return;
        
        CGAL_precondition(begin->ptr()->next==begin->ptr());
        CGAL_precondition(begin->ptr()->prev==begin->ptr());
        
        for(Algebraic_real_iterator it = begin; it != end; ++it){
            CGAL_precondition(it->ptr()->next==it->ptr());
            CGAL_precondition(it->ptr()->prev==it->ptr());
            CGAL_precondition(it->polynomial()==begin->polynomial());
            it->ptr()->next         =begin->ptr();
            it->ptr()->prev         =begin->ptr()->prev;
            begin->ptr()->prev->next=it->ptr();
            begin->ptr()->prev      =it->ptr();  
        }
    }


private:
    CGAL::Comparison_result 
    intern_compare(const Self& y, bool distinct) const {
        if(distinct)            
            CGAL_precondition(this->ptr()!=y.ptr());
        else if(this->ptr()==y.ptr()) { 
            return CGAL::EQUAL;
        }
        CGAL::Comparison_result result = 
            this->ptr()->compare(*(y.ptr()), distinct);
        if (result == CGAL::EQUAL) {
            unify(y); 
        }
        return result;
    }
    
    CGAL::Comparison_result
    intern_compare(const Rational& y, bool distinct) const { 
        return this->ptr()->compare(y,distinct); 
    }
    
    template < class NTX >
    CGAL::Comparison_result
    intern_compare(const NTX& y, bool distinct) const {                 
        if(y           < NTX(low())) return CGAL::LARGER;
        if(NTX(high()) < y         ) return CGAL::SMALLER;

        if(!distinct){
            typename Real_embeddable_traits<NTX>::To_interval to_interval;
            if( CGAL::Interval_nt<true>( this->to_interval() ).do_overlap( 
                              CGAL::Interval_nt<true>(to_interval(y)) ) ) {
                if(polynomial().sign_at(y)==CGAL::ZERO)
                    return CGAL::EQUAL;
            }
        }

        while(NTX(low()) <= y && y <= NTX(high())) {
             refine();
        }
        if(y            <  NTX(low())) return CGAL::LARGER;
        CGAL_assertion(NTX(high())  <  y); return CGAL::SMALLER;
    }   

public:
    //! returns if y is contained in the \c closed isolating interval 
    template < class NTX >
    bool contains(const NTX& y) const { 
        return ((NTX(low()) <= y) && (y <= NTX(high())));
    }

    //! return if \a x is a root of Q
    bool is_root_of(const Polynomial& Q) const {return this->ptr()->is_root_of(Q); }

    /*! \brief returns a rational (Rational) between this number and \c y.
     *  \pre x != y 
     */
    Rational rational_between (const Self& y) const{
        CGAL::Comparison_result s = compare(y);
        CGAL_precondition(s != CGAL::EQUAL);
        if(s == CGAL::SMALLER){
            Rational r((high()+y.low())/Rational(2));
            CGAL::simplify(r);
            return r;
        }else{
            Rational r((y.high()+low())/Rational(2));
            CGAL::simplify(r);
            return r;
        }
    }

public:
    /*! \brief refines the isolating interval to ]<var>lo</var>,
     *  <var>hi</var>[. 
     *
     *  This function can be used to inform an Algebraic_real_pure \a x of an
     *  externally refined isolating interval. Its arguments must be
     *  the boundaries of an isolating interval that contains \a x in its
     *  interior. (Use other functions like \c .strong_refine() in case you
     *  want to communicate an explicit value.)
     */
    void refine_to(const Rational& lo, const Rational& hi) const {
        // test whether lo < x < hi 
        // and refines isolating interval until in ]lo,hi[
        CGAL::Comparison_result s;
        s  = compare_distinct(lo); CGAL_assertion(CGAL::LARGER == s);
        s = compare_distinct(hi) ; CGAL_assertion(CGAL::SMALLER == s);
    }


public:
    template <class NTX>
    bool operator==( const NTX& y) const {return compare(y)==CGAL::EQUAL;}
    template <class NTX>
    bool operator!=( const NTX& y) const {return compare(y)!=CGAL::EQUAL;}
    template <class NTX>
    bool operator< ( const NTX& y) const {return compare(y)==CGAL::SMALLER;}
    template <class NTX>
    bool operator> ( const NTX& y) const {return compare(y)==CGAL::LARGER;}
    template <class NTX>
    bool operator<=( const NTX& y) const {return compare(y)!=CGAL::LARGER;}
    template <class NTX>
    bool operator>=( const NTX& y) const {return compare(y)!=CGAL::SMALLER;}

    //! unary operator +
    const Self& operator+() const { return *this; }

    //! unary operator -
    Self operator-() const {
        Polynomial P(polynomial());
        P.scale_up(Coefficient(-1));
        Rational high_(-low());
        Rational low_ (-high());
        return Self(P,low_,high_);
    }
    
    //! Simplifies the algebraic number
    void simplify() const {
        this->ptr()->simplify();
    }

}; // class Algebraic_real_pure 

} // namespace CGALi 


//----------------------------------------------------------

/*! \ingroup NiX_Algebraic_real_pure
 *  \ingroup NiX_NT_traits_spec
 *  \brief NT_traits class for NiX::Algebraic_real_pure, which is a model of the
 *  RealComparable concept. 
 *
 *  NiX::Algebraic_real_pure does not support any arithmetic operations, thus they
 *  are not even a model of the IntegralDomainWithoutDiv concept. \see NiX_NT_Concepts 
 */
 template< class Coefficient, class Rational, class HandlePolicy, class RepClass >
 class Real_embeddable_traits< CGALi::Algebraic_real_pure< Coefficient, Rational, HandlePolicy, RepClass > >
   : public INTERN_RET::Real_embeddable_traits_base< CGALi::Algebraic_real_pure< Coefficient, Rational, HandlePolicy, RepClass > , CGAL::Tag_true > {

  public:
    
    typedef CGALi::Algebraic_real_pure< Coefficient, Rational, HandlePolicy, RepClass > Type;
        
    class Compare
        : public std::binary_function< Type, Type, CGAL::Comparison_result > {
        public:
            CGAL::Comparison_result operator()( const Type& a, const Type& b ) const {
                return a.compare( b );
            }        
    };
    
    class Sgn
        : public std::unary_function< Type, CGAL::Sign > {
        public:
            CGAL::Sign operator()( const Type& a ) const {
                return a.compare( Rational(0) );
            }
    };
    
    class To_double 
        : public std::unary_function< Type, double > {
        public:
            double operator()(const Type& a) const {
                return a.to_double();
            }
    };
    
    class To_interval 
        : public std::unary_function< Type, std::pair<double, double> > {
        public:
            typename std::pair<double, double> operator()(const Type& a) const {         
                return a.to_interval(); 
            }
    };        
};
  
/*! \relates NiX::Algebraic_real_pure
 *  \brief outputs \c x to \c os 
 */
/*template<class Coefficient, class Rational, class HandlePolicy>
std::ostream& 
operator << (std::ostream& os, 
             const Algebraic_real_pure<Coefficient,Rational,HandlePolicy>& x){
    os << 2 << " " << x.polynomial() << " ";
    os << "[ " << oformat(x.low()) << " , " << oformat(x.high()) << " ]";
    return os;
}*/
/*! \relates NiX::Algebraic_real_pure
 *  \brief read an NiX::Algebraic_real_pure from \c is into \c x.  
 */ 
/*template<class Coefficient, class Rational, class HandlePolicy>
std::istream& 
operator >> (std::istream& is, 
             Algebraic_real_pure<Coefficient, Rational,HandlePolicy>& x){  
    
    typedef Algebraic_real_pure<Coefficient,Rational,HandlePolicy> ALGNUM;    
    int i;    
    is >> i;
    
    // TODO: switch durch is_rational() und else case 
    switch (i){
    case CGAL::IS_RATIONAL:{
        Rational rational;
        is >> iformat(rational);
        x=ALGNUM(rational);
        break;
    }
    case CGAL::IS_GENERAL:{
        Rational low, high;
        typename CGAL::Polynomial_type_generator<Coefficient,1>::Type poly;
        is >> poly;
        swallow(is, '[');// read the "["
        is >> iformat(low);
        swallow(is, ',');// read the ","
        is >> iformat(high);
        swallow(is, ']');// read the "]"
        x = ALGNUM(poly, low, high);
        break;
    }
    default:{
        CGAL_error_msg("ALGEBRAIC_TYPE unexpected!" );
    }
    }
    return is;
}*/

template <class COEFF, class RAT, class POLICY, class REPCLASS >
typename Get_arithmetic_kernel<COEFF>::Arithmetic_kernel::Bigfloat_interval
inline
convert_to_bfi(const CGALi::Algebraic_real_pure< COEFF, RAT, POLICY, REPCLASS >& x){
    typedef typename Get_arithmetic_kernel<COEFF>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI; 
    typedef CGALi::Algebraic_real_pure< COEFF, RAT, POLICY, REPCLASS > ALG;

    if (x.is_rational()) return convert_to_bfi(x.rational());
    
    if(CGAL::sign(x) == CGAL::ZERO) return (BFI(0));
    
    CGAL_postcondition(CGAL::sign(x.low()) == CGAL::sign(x.high()));
    long final_prec = set_precision( BFI(),get_precision( BFI())+4);
  
    BFI bfi = CGAL::hull(convert_to_bfi(x.low()), convert_to_bfi(x.high()));
    
    while( !singleton(bfi) &&  get_significant_bits(bfi) < final_prec  ){
        x.refine();
        bfi = CGAL::hull(
                convert_to_bfi(x.low()), 
                convert_to_bfi(x.high()));
    }

    set_precision(BFI(),final_prec);
    return bfi; 
}

 
template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
struct Get_arithmetic_kernel<CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class> >{

    typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel 
        Arithmetic_kernel;

};

template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
inline CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class> 
min BOOST_PREVENT_MACRO_SUBSTITUTION(
const CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class>& x,
const CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class>& y){
  return (x<=y)?x:y; 
}
template <class Coefficient,class Rational,class Handle_policy,class Rep_class>
inline CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class> 
max BOOST_PREVENT_MACRO_SUBSTITUTION(
const CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class>& x,
const CGAL::CGALi::Algebraic_real_pure<Coefficient,Rational,Handle_policy,Rep_class>& y){
  return (x>=y)?x:y; 
}

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H

// EOF
