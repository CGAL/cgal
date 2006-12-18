// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Elmar Schoemer <schoemer@mpi-inf.mpg.de>
//                 Michael Hemmer <mhemmer@uni-mainz.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file CGAL/Sqrt_extension.h
    \brief Defines class CGAL::Sqrt_extension. 
    
    Provides the Number type Sqrt_extension that extends a given number type  
    \c NT by a square root. One can add several square roots recursively.

*/

#ifndef CGAL_SQRT_EXTENSION_H
#define CGAL_SQRT_EXTENSION_H

#include <CGAL/number_type_basic.h>

#include <numeric> // fro std::accumulate
#include <boost/numeric/interval.hpp> // Needed by To_interval

//#include <boost/iterator/transform_iterator.hpp>
//#include <boost/mpl/if.hpp>

//#include <CGAL/CGAL/number_type_basic.h>

// We have to define the macros befor including Polynomials, 
// since they cause a doxygen error otherwise.. (version 1.2.4)
// the error does not appear for doxygen version 1.2.6
/*! 
  \ingroup CGAL_Arithmetic_traits
  \brief locally define names from \c Arithmetic_traits_sqrt_extension
  
  This macro helps to keep your source code readable if used in the
  following way:
  
  <PRE>
  template \< class ArithmeticTraitsSqrtExtension \>
  class my_class {
  public:
      CGAL_SNAP_ARITHMETIC_TRAITS_SQRT_EXTENSION_TYPEDEFS(ArithmeticTraitsSqrtExtension);
      // ...
  };
  </PRE>
  
  It will declare all typedefs from the ArithmeticTraitsSqrtExtension template
  argument within the class. This makes them accessible for users
  of your class and saves typing the lengthy ArithmeticTraits::
  prefix.
*/
#define CGAL_SNAP_ARITHMETIC_TRAITS_SQRT_EXTENSION_TYPEDEFS(AT) \
  CGAL_SNAP_ARITHMETIC_TRAITS_TYPEDEFS(AT) \
  typedef typename AT::Extn Extn; \
  typedef typename AT::Nested_extn Nested_extn; \
  typedef typename AT::Poly_extn1 Poly_extn1; \
  typedef typename AT::Poly_extn2 Poly_extn2; \
  typedef typename AT::Poly_extn3 Poly_extn3; \
  typedef typename AT::Poly_nested_extn1 Poly_nested_extn1; \
  typedef typename AT::Poly_nested_extn2 Poly_nested_extn2; \
  typedef typename AT::Poly_nested_extn3 Poly_nested_extn3; \


/*#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Cofraction_traits.h>
#include <CGAL/number_type_utils.h>
#include <CGAL/Modular.h>
#include <functional> // only for conversion function at the end
#include <CGAL/Coercion_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Arithmetic_traits.h>
#include <CGAL/Algebraic_extension_traits.h>*/
//#include <CGAL/LiS/Mapping_iterator.h> // only for conversion function at the end
//#include <functional>
//#include <numeric> // only for Algebraic_extension_traits (std::accumulate)

CGAL_BEGIN_NAMESPACE

// From polynomial.h TODO: Where to put this?
inline static void swallow(std::istream &is, char d) {
    char c;
    do c = is.get(); while (isspace(c));
    if (c != d) CGAL_assertion_msg( false, 
                             "input error: unexpected character in polynomial");
}
//// END: From polynomial.h




template <class NT,class ROOT> class Sqrt_extension;
    
template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator + (const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator + (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator - (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> Sqrt_extension<NT,ROOT>
operator * (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template <class NT,class ROOT> inline Sqrt_extension<NT,ROOT>
operator / (const Sqrt_extension<NT,ROOT>&, const Sqrt_extension<NT,ROOT>&);
template<class NT,class ROOT>  
std::ostream& operator << (std::ostream& os, const Sqrt_extension<NT,ROOT>& p);
template <class NT,class ROOT>  
std::istream& operator >> (std::istream& is, Sqrt_extension<NT,ROOT>& p);

/*! \ingroup CGAL_Sqrt_extension
\brief represents an extension of a number type by one square root. 

 An instance of this class
represents  an extension of the type NT by a square root of the 
type ROOT. In case NT and ROOT do not coincide, 
NT must be constructible from ROOT.  The number type NT 
must be at least a model of the IntegralDomainWithoutDiv concept. 

An Sqrt_extension is a model of RealComparable if NT is RealComparable.\n  
The <B>algebraic type</B> of CGAL::Sqrt_extension depends on the algebraic type 
of NT: 
- IntegralDomainWithoutDiv -> IntegralDomainWithoutDiv
- IntegralDomain           -> IntegralDomain
- UFDomain                 -> IntegralDomain
- EuclideanRing            -> IntegralDomain
- Field                    -> Field
- FieldWithSqrt            -> Field (not recommended)


Note that NT and ROOT can themselves be an instance of
CGAL::Sqrt_extension, yielding a nested extension.\n
Note that the extension of an UFDomain or EuclideanRing is just an 
IntegralDomain, since the extension in general destroys the unique 
factorization property. 
*/
template <class NT_,class ROOT_> 
class Sqrt_extension {
public:
    typedef NT_ NT;
    typedef ROOT_ ROOT;
    typedef Sqrt_extension<NT,ROOT> Self;
private:
    NT a0_;
    NT a1_;
    ROOT root_;
    bool is_extended_;
    
    typedef CGAL::Algebraic_structure_traits<NT> Algebraic_structure_traits_nt;
    typedef CGAL::Real_embeddable_traits<NT> Real_embeddable_traits_nt;
public:
    //! Default constructor of \c Sqrt_extension 
    Sqrt_extension()
        : a0_( NT(0)), a1_( NT(0)), root_(ROOT(0)), is_extended_(false) {}
  
    /*!\brief Explicit constructor of Sqrt_extension, from any type NTX.
     * \pre NT must constructible from NTX */ 
    template <class NTX>
    explicit Sqrt_extension(const NTX& i) 
        : a0_(NT(i)), a1_(NT(0)), root_(ROOT(0)), is_extended_(false) {}

    /*! \brief copy constructor  */ 
    Sqrt_extension(const Sqrt_extension<NT,ROOT>& x)
        : a0_(x.a0()),
          a1_(x.a1()),
          root_(x.root()),
          is_extended_(x.is_extended()){}
    

    /*! \brief Expicite constructor of Sqrt_extension, from 
     *  \c Sqrt_extension<NTX,ROOTX>.
     *  \pre NT must constructible from NTX 
     *  \pre ROOT must constructible from ROOTX */ 
//     template <class NTX,class ROOTX>
//     explicit Sqrt_extension(const Sqrt_extension<NTX,ROOTX>& x)
//         : a0_(x.a0()),
//           a1_(x.a1()),
//           root_(x.root()),
//           is_extended_(x.is_extended()){}

    /*! \brief Constructor from some type NTX and ROOTX
     *  NT must be constructible from NTX\\
     *  ROOT must be construcible from ROOTX\\
     */ 
    template <class NTX,class ROOTX>
    explicit Sqrt_extension(const NTX& a0, const NTX& a1, const ROOTX& root)
        : a0_(a0),
          a1_(a1),
          root_(root),
          is_extended_(true){
        CGAL_precondition( !CGAL_NTS is_zero(root));
    }
    
    //! Access operator for a0_, \c const 
    const NT& a0() const { return a0_; }
    //! Access operator for a0_
    NT&        a0()       { return a0_; }
    //! Access operator for a1_, \c const 
    const NT& a1() const { return a1_; }
    //! Access operator for a1_
    NT&        a1()       { return a1_; }
    //! Access operator for root_, \c const 
    const ROOT& root() const { return root_; }
    //! Access operator for is_extended_, \c const 
    const bool& is_extended() const { return is_extended_; }
    //! Access operator for root_ 
    //ROOT& root() { return root_; }
  
public:
    
    //! propagates the simplify command to the members of xx   
    void simplify(){
        CGAL_NTS simplify(a0_);
        CGAL_NTS simplify(a1_);
        CGAL_NTS simplify(root_);
    }
    
    //! determines the sign of xx by repeated squaring.
    ::CGAL::Sign sign() const { 
        if (!is_extended()) 
            return CGAL_NTS sign(a0());
        
        ::CGAL::Sign s0,s1;

        s0 = CGAL_NTS sign(a0());
        s1 = CGAL_NTS sign(a1());

        if (s0 == s1) return s0;
        if (s0 == CGAL::ZERO) return s1;
        if (s1 == CGAL::ZERO) return s0;
        
        // s0*s1=-1
        NT r = a1()*a1()*NT(root())-a0()*a0();
        // if(r>0) return s1 else s0
        if (s1 == CGAL::POSITIVE)
            return CGAL_NTS sign(r);
        else
            return -CGAL_NTS sign(r); // TODO: Is this valid??? Was: -CGAL::sign(..)
    }
  

    
    template < class BOOL_TAG > 
    bool is_zero_(const BOOL_TAG&) const{
        // Is_real_comparable == CGAL::Tag_true 
        return sign()==CGAL::ZERO;
    }
    bool is_zero_(const CGAL::Tag_false) const {
        // Is_real_comparable == CGAL::Tag_false
        // i.e. CGAL::Modular
        if(is_extended()){
            if(a0() == (NT)0 && a1()== (NT)0) {
                return true;
            }else{
                return (a0()*a0()-a1()*a1()*(NT)root() == NT(0));
            }
        }else{
            return a0() == NT(0);
        }
    }

    //! returns \a true if xx is zero 
    bool is_zero() const { 
//        typedef typename Algebraic_structure_traits_nt::Is_real_comparable Is_real_comparable;
      typedef typename Real_embeddable_traits_nt::Is_real_embeddable 
                                                             Is_real_embeddable;
        return is_zero_(Is_real_embeddable());
    }

    //! returns the absolute value of xx
    Self abs() const { 
        if (sign() == CGAL::NEGATIVE) 
            return -*this; 
        return *this; 
    }
    
    friend Self operator + <>(const Self&);
    friend Self operator - <>(const Self&); 
    friend Self operator + <>(const Self&, const Self&);
    friend Self operator - <>(const Self&, const Self&);
    friend Self operator * <>(const Self&, const Self&);
    friend Self operator / <>(const Self& p1, const Self& p2);

    Self& operator += (const Self& p) {
        (*this)=(*this)+p; 
        return (*this); 
    }
    Self& operator -= (const Self& p){ 
        (*this)=(*this)-p; 
        return (*this); 
    }
    Self& operator *= (const Self& p){ 
        (*this)=(*this)*p; 
        return (*this); 
    }
    Self& operator /= (const Self& p) { 
        (*this)=(*this)/p; 
        return (*this); 
    }

//------------------------------------------------------------------
// SPECIALIZE_MEMBERS
    
    Self& operator += (const NT& num) { 
        a0() += NT(num); 
        return *this; 
    }
    Self& operator -= (const NT& num) { 
        a0() -= NT(num); 
        return *this; 
    }
    Self& operator *= (const NT& num) { 
        a0() *= NT(num); 
        a1() *= NT(num);
        return *this; 
    }
    Self& operator /= (const NT& num) { 
        typename Real_embeddable_traits_nt::Sign sign;       
        CGAL_assert(sign(num) != 0);
        a0() /= NT(num); 
        a1() /= NT(num);
        return *this; 
    }
public:    
    void output_maple(std::ostream& os) const;
public:    
    /*! \brief write Sqrt_extension to \c os in a format readable 
         by \c input_ascii()
    
         The output format is: 
         <TT><B>EXT[</B></TT><I>a0</I><TT><B>,</B></TT>
         <I>a1</I><TT><B>,</B></TT><I>root</I><TT><B>]</B></TT> 
     */
   void output_ascii(std::ostream& os) const{
       os<<"EXT["<<a0()<<","<<a1()<<","<<root()<<"]";
   }
    /*! \relates CGAL::Sqrt_extension
     *  \brief read a Sqrt_extension from \c is
     *
     * The expected format is:
     * <TT><B>EXT[</B></TT><I>a0</I><TT><B>,</B></TT>
     * <I>a1</I><TT><B>,</B></TT><I>root</I><TT><B>]</B></TT> 
     * 
     * The format of the coefficients must be understandable for
     * <TT> is >> iformat(ai) </TT>.
     *
     * Example: A \c CGAL::Sqrt_extension<int,Root<1,int> > with a value of
     * 4-2*sqrt(5) has to be written as
     * \c <TT>Self[4,-2,5]</TT> for this function.
     */
    static Sqrt_extension<NT,ROOT> input_ascii(std::istream& is);
};

template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator + (const Sqrt_extension<NT,ROOT>& p) { return p; }

template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator - (const Sqrt_extension<NT,ROOT>& p){
    if(p.is_extended())
        return  Sqrt_extension<NT,ROOT>(-p.a0(),-p.a1(),p.root()); 
    else 
        return  Sqrt_extension<NT,ROOT>(-p.a0()); 
}

template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator + (const Sqrt_extension<NT,ROOT>& p1, 
            const Sqrt_extension<NT,ROOT>& p2){
    typedef Sqrt_extension<NT,ROOT> EXT;
    if(p1.is_extended()){
        if (p2.is_extended()){CGAL_precondition(p1.root()==p2.root());}
        return EXT (p1.a0()+p2.a0(), p1.a1()+p2.a1(), p1.root());
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()+p2.a0(), p2.a1(), p2.root());
        else
            return EXT (p1.a0()+p2.a0()); 
    }
} 

template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator - (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) { 
    typedef Sqrt_extension<NT,ROOT> EXT;
    if(p1.is_extended()){
        if (p2.is_extended()){CGAL_precondition(p1.root()==p2.root());}
        return EXT (p1.a0()-p2.a0(), p1.a1()-p2.a1(), p1.root());
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()-p2.a0(), -p2.a1(), p2.root());
        else
            return EXT (p1.a0()-p2.a0()); 
    } 
}

template <class NT,class ROOT> Sqrt_extension<NT,ROOT> 
operator * (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) {
    typedef Sqrt_extension<NT,ROOT> EXT;
    
    if(p1.is_extended()){
        if (p2.is_extended()){
            CGAL_precondition(p1.root()==p2.root());
            return EXT (p1.a0()*p2.a0()+p1.a1()*p2.a1()*NT(p1.root()),
                        p1.a0()*p2.a1()+p1.a1()*p2.a0(),
                        p1.root());
        }else{
            return EXT (p1.a0()*p2.a0(), p1.a1()*p2.a0(), p1.root()); 
        }
    }else{
        if (p2.is_extended())
            return EXT (p1.a0()*p2.a0(), p1.a0()*p2.a1(), p2.root()); 
        else
            return EXT (p1.a0()*p2.a0()); 
    } 
}

template <class NT,class ROOT> inline Sqrt_extension<NT,ROOT> 
operator / (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2) { 
    
    typedef Sqrt_extension<NT,ROOT> EXT;
    CGAL_assertion(! p2.is_zero());
    typename CGAL::Algebraic_structure_traits<NT>::Integral_division Idiv;
   
    if(p1.is_extended()){
        if (p2.is_extended()){
            CGAL_precondition(p1.root()==p2.root());
            NT c = p2.a0()*p2.a0()-p2.a1()*p2.a1()*NT(p2.root());
            if (c == NT(0)) {                          //TR
                NT a0 = Idiv(p1.a0(), NT(2)*p2.a0()) 
                      + Idiv(p1.a1(), NT(2)*p2.a1());
                return EXT(a0);
            }
            NT a0 = Idiv(p1.a0()*p2.a0()-p1.a1()*p2.a1()*NT(p2.root()),c);    
            NT a1 = Idiv(p1.a1()*p2.a0()-p1.a0()*p2.a1(),c);    
            return EXT(a0,a1,p1.root());
        }else{
            NT a0 = Idiv(p1.a0(),p2.a0());    
            NT a1 = Idiv(p1.a1(),p2.a0());    
            return EXT(a0,a1,p1.root());
        }
    }else{
        if (p2.is_extended()){
            NT c = p2.a0()*p2.a0()-p2.a1()*p2.a1()*NT(p2.root());
            if(c == NT(0)){
                NT a0 = Idiv(p1.a0(), NT(2)*p2.a0());
                return EXT(a0);
            }
            NT a0 = Idiv(p1.a0()*p2.a0(),c);    
            NT a1 = Idiv(-p1.a0()*p2.a1(),c);  
            return EXT(a0,a1,p2.root()); 
        }else{
            return EXT (Idiv(p1.a0(),p2.a0())); 
        }
    } 
}

template <class NT,class ROOT> bool 
operator == (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return (p1-p2).is_zero() ; }   

template <class NT,class ROOT> bool 
operator != (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return !(p1-p2).is_zero() ; }    

template <class NT,class ROOT> bool 
operator <  (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() < 0 ); }    

template <class NT,class ROOT> bool 
operator <= (const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() <= 0 ); }    

template <class NT,class ROOT> bool operator >  
(const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
{ return ( (p1-p2).sign() > 0 ); }    

template <class NT,class ROOT> bool operator >= 
(const Sqrt_extension<NT,ROOT>& p1, const Sqrt_extension<NT,ROOT>& p2)
  { return ( (p1-p2).sign() >= 0 ); }    

// lefthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator + 
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) + p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator - 
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) - p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator * 
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num) * p2); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator / 
(const NT& num, const Sqrt_extension<NT,ROOT>& p2)
{ return (Sqrt_extension<NT,ROOT>(num)/p2); }

// righthand side
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator + 
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{ return (p1 + Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator - 
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{ return (p1 - Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator * 
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{ return (p1 * Sqrt_extension<NT,ROOT>(num)); }
template <class NT,class ROOT>    Sqrt_extension<NT,ROOT> operator / 
(const Sqrt_extension<NT,ROOT>& p1, const NT& num)
{ return (p1 / Sqrt_extension<NT,ROOT>(num)); }

// lefthand side
template <class NT,class ROOT>    bool operator ==  
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() == 0 );}
template <class NT,class ROOT>    bool operator != 
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() != 0 );}
template <class NT,class ROOT>    bool operator <  
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() < 0 );}
template <class NT,class ROOT>    bool operator <=  
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() <= 0 );}
template <class NT,class ROOT>    bool operator >  
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=  
(const NT& num, const Sqrt_extension<NT,ROOT>& p) 
{ return ( (Sqrt_extension<NT,ROOT>(num)-p).sign() >= 0 );}

// righthand side
template <class NT,class ROOT>    bool operator ==
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() == 0 );}
template <class NT,class ROOT>    bool operator !=
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() != 0 );}
template <class NT,class ROOT>    bool operator < 
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() < 0 );}
template <class NT,class ROOT>    bool operator <= 
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() <= 0 );}
template <class NT,class ROOT>    bool operator > 
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() > 0 );}
template <class NT,class ROOT>    bool operator >=
(const Sqrt_extension<NT,ROOT>& p, const NT& num) 
{ return ( (p-Sqrt_extension<NT,ROOT>(num)).sign() >= 0 );}

// Algebraic structure traits
template< class Type, class Algebraic_type >
class Sqrt_extension_algebraic_structure_traits_base;

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                        CGAL::Integral_domain_without_division_tag >
  : public Algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_without_division_tag Algebraic_category;
    
    class Simplify 
      : public Unary_function< Type&, void > {
      public:
        typedef void result_type;
        typedef Type& argument_type;

        void operator()( Type& x ) const {
          x.simplify();
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                    CGAL::Integral_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_without_division_tag > {
  public:
    typedef CGAL::Integral_domain_tag Algebraic_category;
    
    class Integral_division 
      : public Binary_function< Type, Type, 
                                Type > {
      public:
        Type operator()( const Type& x, 
                                        const Type& y ) const {
          return x/y;
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR( Type )
    };                                  
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                    CGAL::Unique_factorization_domain_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                    CGAL::Euclidean_ring_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                    CGAL::Field_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Integral_domain_tag > {
  public:
    typedef Field_tag Algebraic_category;
    
    class Unit_part 
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return( x == Type(0) ? Type(1) : x );
        }
    };
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                    CGAL::Field_with_sqrt_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      CGAL::Field_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                CGAL::Field_with_kth_root_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class Type >
class Sqrt_extension_algebraic_structure_traits_base< Type, 
                                                CGAL::Field_with_root_of_tag >
  : public Sqrt_extension_algebraic_structure_traits_base< Type,
                                      // TODO: Why not Fiel_tag?
                                      CGAL::Field_with_sqrt_tag > {
  // Nothing new
};

template< class COEFF, class ROOT>
class Algebraic_structure_traits< Sqrt_extension< COEFF, ROOT > >
  : public Sqrt_extension_algebraic_structure_traits_base< 
      Sqrt_extension< COEFF, ROOT >,
      typename Algebraic_structure_traits< COEFF >::Algebraic_category > {
  public:
    typedef Sqrt_extension< COEFF, ROOT > Type;
    
    // Tag_true if COEFF and ROOT are exact 
    typedef typename ::boost::mpl::if_c< 
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<ROOT >::Is_exact,::CGAL::Tag_true>::value )&&
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<COEFF>::Is_exact,::CGAL::Tag_true>::value )
           ,::CGAL::Tag_true,::CGAL::Tag_false>::type Is_exact;
    
    typedef typename Algebraic_structure_traits<COEFF>::Is_numerical_sensitive 
    Is_numerical_sensitive;
};

//
// Real embeddable traits
//

template< class COEFF, class ROOT > 
class Real_embeddable_traits< Sqrt_extension<COEFF, ROOT> > 
  : public INTERN_RET::Real_embeddable_traits_base_selector< 
                  Sqrt_extension<COEFF, ROOT>, 
                  typename Real_embeddable_traits<COEFF>::Is_real_embeddable > {
  public:
    typedef Sqrt_extension<COEFF, ROOT> Type;  
           
    class To_interval 
      : public Unary_function< Type, std::pair< double, double > > {
      private:
        typedef boost::numeric::interval<double> Double_interval;
      public:
        std::pair<double,double> operator()(const Type& x) const {
            if(x.is_extended()){
              std::pair<double, double> pair_a0 = CGAL_NTS to_interval( x.a0() );
              std::pair<double, double> pair_a1_root = CGAL_NTS to_interval(
                                         x.a1() * x.a1() * COEFF( x.root() ) );
                                         
              Double_interval result( pair_a0.first, pair_a0.second );
              result += ( Double_interval( (int) CGAL_NTS sign(x.a1())) * 
                      ::boost::numeric::sqrt( 
                              Double_interval( pair_a1_root.first,
                                               pair_a1_root.second ) ) ); 
              return std::pair<double, double>(result.lower(), result.upper() );
            } else {
                return CGAL_NTS to_interval( x.a0());
            }
        }
    };

    class To_double 
      : public Unary_function< Type, double > {
      public:
        // The main problem here is, that even tough the total 
        // expression fits into double, one of the coefficients 
        // or the root may not. ?? ! 
        double operator()(const Type& x) const {
            if(x.is_extended()){
                return CGAL_NTS to_double(x.a0()) 
                    +  int(CGAL_NTS sign(x.a1()))
                    * CGAL_NTS sqrt(CGAL_NTS to_double(x.a1()*x.a1() * 
                                                    Type(x.root()))); 
            }else{
                return CGAL_NTS to_double(x.a0());
            }
        }
    };
};



//##########################Algebraic_structure_traits<Sqrt_extension>######################


// ############## initializing for Algebraic_structure_traits_comparable_base 
/*template <class EXT, class TypeTAG>
class Sqrt_ext_NTtr_comp_base;

template <class NT, class ROOT >
class Sqrt_ext_NTtr_comp_base<Sqrt_extension<NT,ROOT> , CGAL::Tag_true >
    :public Algebraic_structure_traits_comparable_base <Sqrt_extension<NT,ROOT> > {
    typedef Sqrt_extension<NT,ROOT> EXT;
public:
    class To_Interval {
    public:
        //! type for the \c AdaptableUnaryFunction concept.
        typedef EXT argument_type;
        //! type for the \c AdaptableUnaryFunction concept.
        typedef Interval result_type;
        //! the function call.
        Interval operator()(const EXT& a) const {
            if(a.is_extended()){               
                return CGAL::to_Interval(a.a0())
                    + Interval(int(CGAL::sign(a.a1())))
                    * CGAL::sqrt(CGAL::to_Interval(a.a1()*a.a1()*NT(a.root()))); 
            }else{
                return CGAL::to_Interval(a.a0());
            }
        }
    };
    class To_double {
    public:
        // the result type.
        typedef double result_type;
        // the argument type.
        typedef EXT argument_type;
        // The main problem here is, that even tough the total 
        // expression fits into double, one of the coefficients 
        // or the root may not. ?? ! 
        double operator()(const EXT& a) const {
            if(a.is_extended()){
                return CGAL::to_double(a.a0()) 
                    +  int(CGAL::sign(a.a1()))
                    *sqrt(CGAL::to_double(a.a1()*a.a1()*NT(a.root()))); 
            }else{
                return CGAL::to_double(a.a0());
            }
        }
    };
    
};
template <class EXT>
class Sqrt_ext_NTtr_comp_base<EXT,CGAL::Tag_false>
    :public Algebraic_structure_traits_comparable_base < LiS::Null_tag > {
// nothing new
};*/


// ############### actual Algebraic_structure_traits for Sqrt_extension

/*! \ingroup CGAL_Sqrt_extension
    \ingroup CGAL_Algebraic_structure_traits_spec
  
    \brief Specialization of CGAL::Algebraic_structure_traits for CGAL::Sqrt_extension. 
*/
/*template < class COEFF, class ROOT >
class Algebraic_structure_traits< Sqrt_extension< COEFF, ROOT > > 
    : public Sqrt_extension_NTtr_base< Sqrt_extension< COEFF, ROOT >, 
                             typename Algebraic_structure_traits< COEFF >::Algebra_type>,
      public Sqrt_ext_NTtr_comp_base< Sqrt_extension< COEFF, ROOT >, 
                             typename Algebraic_structure_traits<COEFF>::Is_real_comparable>,
      public Algebraic_structure_traits_log2_null_base
{
public:
    typedef Sqrt_extension< COEFF, ROOT > NT;  

    // Tag_true if COEFF and ROOT are exact 
    typedef typename ::boost::mpl::if_c< 
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<ROOT >::Is_exact,::CGAL::Tag_true>::value )&&
       bool( ::boost::is_same<typename CGAL::Algebraic_structure_traits<COEFF>::Is_exact,::CGAL::Tag_true>::value )
       ,::CGAL::Tag_true,::CGAL::Tag_false>::type Is_exact;
};*/


// #################### IO BEGIN ############################################
// #################### INPUT


template<class NT, class ROOT> 
Sqrt_extension<NT,ROOT> 
Sqrt_extension<NT,ROOT>::input_ascii(std::istream& is){
    // expected input format: EXT[ext.a0(),ext.a1(),root()]
    typedef Sqrt_extension<NT,ROOT> EXT;
    char c;
    NT a0;
    NT a1;
    ROOT root;

    swallow(is, 'E');
    swallow(is, 'X');
    swallow(is, 'T');
    swallow(is, '[');
    is >> iformat(a0);
    do c = is.get(); while (isspace(c));
    // TODO: Replace CGAL_assertion_msg( false, ... ) with CGAL_error
    if (c != ',') CGAL_error( "input error: , expected" );

    is >> iformat(a1);
    do c = is.get(); while (isspace(c));
    if (c != ',') CGAL_error( "input error: , expected" );

    is >> iformat(root);
    do c = is.get(); while (isspace(c));
    if (c != ']') CGAL_error( "input error: ] expected" );
    
    if ( root  < ROOT(0)) CGAL_error("input error: non-negative root expected");
    if ( root == ROOT(0)) return EXT(a0);
    return EXT(a0,a1,root);
}

// ##################### OUTPUT
//! write Sqrt_extension to \c os in \c LiS::IO::PRETTY format
/*! The output is intended to be Maple-readable; see module
 *  \link CGAL_io CGAL I/O Support \endlink.
 */
template<class NT, class ROOT> 
void 
Sqrt_extension<NT,ROOT>::output_maple(std::ostream& os) const{
    CGAL::IO::Mode o_mode=::CGAL::get_mode(os);        
    ::CGAL::set_mode(os,CGAL::IO::PRETTY); 

    if ( a0() != NT(0)){
        if ( a1() != NT(0)){               
            os << a0() 
               << "+" << CGAL::oformat(a1(),CGAL::Parens_as_product_tag()) 
               << "*sqrt(" << root() << ")";
        }else{
            os << a0(); 
        }
    }
    else{
        if (a1() != NT(0)){
            os << CGAL::oformat(a1(),CGAL::Parens_as_product_tag()) 
               << "*sqrt(" << root() << ")"; 
        }else{
            os << 0; 
        }
    }
    ::CGAL::set_mode(os,o_mode);   
    return;
}
 
/*! \relates CGAL::Sqrt_extension
 *  \brief output \c ext to \c os
 *
 *  Output \c ext in a format as specified by
 *  \c LiS::get_mode(os), see \link LiS_io LiS I/O Support \endlink.
 *  Currently, the output for \c LiS::IO::BINARY happens to be
 *  identical to \c LiS::IO::ASCII.
 */
template <class NT,class ROOT>
std::ostream& operator << (std::ostream& os, 
                           const Sqrt_extension<NT,ROOT>& ext){
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        ext.output_maple(os); break;
    default:
        ext.output_ascii(os); break;
    }
    return os;
}

/*! \relates CGAL::Sqrt_extension
 *  \brief try to read a CGAL::Sqrt_extension from \c is into \c ext
 *
 *  \c is must be in a mode that supports input of CGAL::Sqrt_extension
 *  (\c LiS::IO::ASCII or \c LiS::IO::BINARY) and the input from
 *  \c is must have the format of output to a stream of the same mode.
 */
template <class NT,class ROOT>
std::istream& operator >> (std::istream& is, Sqrt_extension<NT,ROOT>& ext) {
    CGAL_precondition(!CGAL::is_pretty(is));
    ext = Sqrt_extension<NT,ROOT>::input_ascii(is);
    return is;
}
// ################################ IO END ###############################

/*!
  \ingroup CGAL_Sqrt_extension
  \ingroup CGAL_Modular_traits_spec
  \brief Specialization of CGAL::Modular_traits for \c CGAL::Sqrt_extension.
  
  A model of the concept ModularTraits.
 
  CGAL::Modular_traits::Modular_image maps the coefficients of a  
  CGAL::Sqrt_extension to their Modular_image and returns the resulting 
  compound type.  
*/
/*template< class COEFF, class ROOT>
class Modular_traits< Sqrt_extension<COEFF,ROOT> > {
    
private:
    typedef Sqrt_extension<COEFF,ROOT> EXT; 
    typedef Modular_traits<COEFF> MT_COEFF;
    typedef Modular_traits<ROOT>  MT_ROOT;
    typedef typename MT_COEFF::Modular_NT Modular_NT_coeff;
    typedef typename MT_ROOT::Modular_NT  Modular_NT_root;
public:
    typedef Sqrt_extension<COEFF, ROOT > NT;
    typedef typename MT_COEFF::Is_convertible Is_convertible;
    typedef Sqrt_extension<Modular_NT_coeff, Modular_NT_root> Modular_NT;
    
    struct Modular_image{
        Modular_NT operator()(const EXT& a){
            typename MT_ROOT::Modular_image mod_image_root;
            typename MT_COEFF::Modular_image mod_image_coeff;
            Modular_NT_root  root_mod = mod_image_root(a.root());
            if(root_mod != Modular_NT_root(0)){
                return Modular_NT(mod_image_coeff(a.a0()),
                                  mod_image_coeff(a.a1()),
                                  root_mod);
            }else{
                return Modular_NT(mod_image_coeff(a.a0()));
            }
        }
    };
};*/

/*! \relates CGAL::Sqrt_extension
 *  \brief convert Sqrt_extension to root expression in a FieldWithSqrt
 *
 *  If \c NT is a type from which numbers of type \c FieldWithSqrt
 *  can be constructed, then this function converts a one-root expression
 *  from \c Sqrt_extension<NT,ROOT> to \c FieldWithSqrt in the canonical
 *  way. In its present form, this function does not work for nested
 *  Sqrt_extensions.
 */
/*template <class FieldWithSqrt, class NT, class ROOT>
FieldWithSqrt to_rootexp(Sqrt_extension<NT, ROOT> x) {
    if (!x.is_extended()) {
        return FieldWithSqrt(x.a0());
    } else {
        return FieldWithSqrt(x.a0())
             + FieldWithSqrt(x.a1())*sqrt(FieldWithSqrt(x.root()));
    }
};*/

/*! \relates CGAL::Sqrt_extension
 *  \brief convert coefficients of univariate Polynomial
 *
 *  Convert a polynomial from \c Sqrt_extension to \c FieldWithSqrt by
 *  applying \c to_rootexp() to each coefficient.
 *  \pre Types must be such that \c to_rootexp() works for the
 *  \c Sqrt_extension type.
 */
/*template <class FieldWithSqrt, class NT, class ROOT>
Polynomial<FieldWithSqrt> to_rootexp(Polynomial< Sqrt_extension<NT, ROOT> > p) {
    typedef Polynomial<FieldWithSqrt> POLY;
    FieldWithSqrt (*conv)(Sqrt_extension<NT, ROOT> ext)
        = &to_rootexp<FieldWithSqrt, NT, ROOT>;
    return POLY(::boost::make_transform_iterator(p.begin(), std::ptr_fun(conv)),
                ::boost::make_transform_iterator(p.end(),   std::ptr_fun(conv))
    );
};*/

/*! \relates CGAL::Sqrt_extension
 *  \brief convert \c CGAL::Algebraic_real into \c Sqrt_extension
 *
 *  Let \c x be an CGAL::Algebraic_real with defining \c x.polynomial()
 *  of degree at most 2. Then \c x can be written as a one-root
 *  expression. This function converts it to a one-root expression,
 *  represented by an object of class \c Sqrt_extension (which is expected
 *  to be a specialization of \c CGAL::Sqrt_extension<> ).
 *
 *  It is the caller's responsibility to take care of the problem
 *  that an algebraic real with defining polynomial of degree 2
 *  may actually be rational.
 *
 *  This function is designed for cases where \c Sqrt_extension::NT models
 *  the rationals, \c Sqrt_extension::ROOT models is integer 
 *  and \c CGAL::Algebraic_real<>::Coefficient is integer
 *  or rational. (You may get away in other cases if you know
 *  what you are doing). If parameter \c simplify_radicant is set to \c true
 *  the function tries to extract square parts out of the root. Note that the
 *  current method is really slow at the moment.
 *
 *  \pre Sqrt_extension::NT must have a constructor from \c Coefficient.
 *  \pre It must be possible to call \c CGAL::Algebraic_real<>.compare()
 *  with an argument of type \c NT .
 */
 // TODO: No Coercion_traits available yet.
/*template <class Extension, class AlgebraicReal >
Extension to_Sqrt_extension(AlgebraicReal x, bool simplify_radicant = false) { 
    
    typedef typename AlgebraicReal::Coefficient Coefficient;
    typedef typename AlgebraicReal::Rational Rational;

    typedef typename Extension::NT NT;
    typedef typename Extension::ROOT ROOT;
    
    Polynomial< Coefficient > p = CGAL::canonicalize_polynomial(x.polynomial());
    switch (p.degree()) {
    case 2: {
        typedef CGAL::Coercion_traits< NT, Coefficient > CT;
        typedef typename CT::Type Type;
        typename CGAL::Algebraic_structure_traits<Type>::integral_division idiv;
        typename CT::Cast cast;
    
        Type den(Type(2) * cast(p[2]));
        CGAL_assert(CGAL::sign(den) == CGAL::POSITIVE);
        Type mid = idiv(cast(-p[1]), den);
        CGAL::simplify(mid);
        ROOT rad = p[1]*p[1] - ROOT(4)*p[0]*p[2];
        CGAL_assert(rad != ROOT(0)); // double root of algebraic real polynomial
        ROOT f = 1;
        Extension e;
        // is radicant a square number?
        typename CGAL::Algebraic_structure_traits< ROOT >::Sqrt isqrt;
        ROOT root = isqrt(rad);
        if (root * root == rad) {
            Type c = idiv(f, den);
            switch (x.compare(mid)) {
            case CGAL::LARGER:             break;
            case CGAL::SMALLER:    c = -c; break;
            default:               CGAL_error("bogus comparison result");
            }
            e = Extension(mid + c * NT(root));
        } else {
            if (simplify_radicant) {
                // root out square parts of radicant -> very slow!!!
                ROOT s = 4;
                ROOT i = 2;
                while (CGAL::compare(s,rad) == CGAL::SMALLER) {
                    typename CGAL::Algebraic_structure_traits< ROOT >::Mod mod;
                    while (mod(rad, s) == ROOT(0)) {
                        typename CGAL::Algebraic_structure_traits< ROOT >::Div div;
                        rad = div(rad,s);
                        f *= i;
                    }
                    s += (2 * i) + 1;
                    i += 1;
                }
            }
            Type c = idiv(f, den);
            switch (x.compare(mid)) {
            case CGAL::LARGER:             break;
            case CGAL::SMALLER:    c = -c; break;
            default:             CGAL_error("bogus comparison result");
            }
            CGAL::simplify(c);
            e = Extension(mid, c, rad);
        }
        CGAL_postcond(p.evaluate(e) == Extension(0));
        CGAL_postcond_code(
                typename AlgebraicReal::Field_with_sqrt x0;
                CGAL::convert_to(e,x0);
        );
        CGAL_postcond(x.real() == x0);
        return e; 
    }
    case 1:
        CGAL_assert(x.is_rational());
        // warning: Extension::ROOT::root is not set!
        return Extension(x.rational());
    default:
        CGAL_error("CGAL::to_Extension(x) called with x.poly.degree != 1, 2");
        // NOT REACHED 
        return Extension();
    }
};*/


//################################# CGAL::Scalar_factor_traits ##################
/*! \ingroup CGAL_Sqrt_extension
    \ingroup CGAL_Scalar_factor_traits_spec
    \brief specialization of CGAL::Scalar_factor_traits for CGAL::Sqrt_extension
 */
template <class COEFF, class INTERNAL>
class Scalar_factor_traits< Sqrt_extension<COEFF, INTERNAL> > {
public:
    
    //! the number type for which this instance has been instantiated
    typedef Sqrt_extension<COEFF, INTERNAL> NT;
      //! the number type of scalars that can be extracted from NT
    typedef typename Scalar_factor_traits<COEFF>::Scalar Scalar;
    
    class Scalar_factor 
    {
    public:
        //! argument type
        typedef NT argument_type;
        //! first argument type
        typedef NT first_argument_type;
        //! second argument type
        typedef Scalar second_argument_type;
        //! result type
        typedef Scalar result_type;
        
        Scalar 
        operator () (const NT& x, const Scalar& d_ = Scalar(0) ) {
            typename Scalar_factor_traits<COEFF>::Scalar_factor sfac;  

            Scalar d(d_);
            Scalar unity(1);
            if(d==unity) return d;
            d=sfac(x.a0(),d);
            if(d==unity) return d;
            if(x.is_extended())
                d=sfac(x.a1(),d);
            return d;
        }
    };
    
    class Scalar_div 
    {
    public:
        //! first_argument_type
        typedef NT first_argument_type;
        //! second_argument_type
        typedef Scalar second_argument_type;
        //! divides an extension \c a by a scalar factor \c b
        void operator () (NT& a, const Scalar& b) {
            CGAL_precondition(b != Scalar(0));
            typename Scalar_factor_traits<COEFF>::Scalar_div sdiv;
            sdiv(a.a0(), b); sdiv(a.a1(), b); // perform division in place
        }
    };
};

//################################# CGAL::Fraction_traits ##################
// Select the right alternative as Fraction_traits
// The actual Type traits is Intern::Sqrt_ext_Ftr_base_2
// The selction is done in two steps:
// 1. Inter::Sqrt_ext_Ftr_base_1 selects by the BOOL_TAG whether the COEFF type
//    Is_fraction 
// 2. Intern::Sqrt_ext_Ftr_base_2 checks whether the internal type of the ROOT
//    is still implicite convertible to the new COEFF type. 
//    since the ROOT type it self can not be converted.
namespace Intern{ 
    template <class EXT, bool> class Sqrt_ext_Ftr_base_2;
    template <class EXT, class BOOL_TAG> class Sqrt_ext_Ftr_base_1;
} 

/*! \ingroup CGAL_Sqrt_extension
    \ingroup CGAL_Fraction_traits_spec
    \brief Specialisation of CGAL::Fraction_traits for CGAL::Sqrt_extension.  
 *
 *  Extensions provide suitable specializations of \c CGAL::Fraction_traits.
 *  They are decomposable iff their coefficient type is.
 *  The denominator \e d of a Extension \e ext is a low common multiple
 *  (see \c CGAL::Fraction_traits::Common_factor for details) of the
 *  denominators of its coefficients.  The numerator is the Extenion
 *  \e d*ext with a fraction-free coefficient type.
 *
 *  This works for nested Sqrt_extensions, too.
 */

template <class COEFF, class ROOT_NT >
class Fraction_traits< Sqrt_extension<COEFF,ROOT_NT > >
    : public Intern::Sqrt_ext_Ftr_base_1< 
    Sqrt_extension<COEFF,ROOT_NT >, 
    typename CGAL::Fraction_traits<COEFF>::Is_fraction >              
{
    // nothing new
};

namespace Intern {

// Use this if the coefficients cannot be decomposed 
// into numerator and denominator
template <class NT_ > 
class Sqrt_ext_Ftr_base_2< NT_, false > {
public:
    typedef NT_ NT;
    typedef ::CGAL::Tag_false Is_fraction;
    typedef ::CGAL::Null_tag Numerator_type;
    typedef ::CGAL::Null_tag Denominator_type;
    typedef ::CGAL::Null_tag Common_factor;
    typedef ::CGAL::Null_tag Decompose;
    typedef ::CGAL::Null_tag Compose;
};

template <class COEFF, class ROOT_NT>
class Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT_NT>, true > {
private:
    typedef Fraction_traits<COEFF> CFT;
public:
    typedef Sqrt_extension<COEFF,ROOT_NT> NT;
    typedef CGAL::Tag_true Is_fraction;
    typedef Sqrt_extension<typename CFT::Numerator_type,ROOT_NT> Numerator_type;
    typedef typename CFT::Denominator_type Denominator_type;
    typedef typename Algebraic_structure_traits<Denominator_type>::Gcd Common_factor;

    class Decompose {
    public:
        typedef NT first_argument_type;
        typedef Numerator_type second_argument_type;
        typedef Denominator_type& third_argument_type;
        void operator () (const NT& ext,
                          Numerator_type&   num,
                          Denominator_type& den){
            typename CFT::Decompose decompose;
            typename CFT::Common_factor common_factor;
            typedef typename CFT::Numerator_type NUM;
            typedef typename CFT::Denominator_type DEN;
            
            if(ext.is_extended()){
                NUM a0_num, a1_num;
                DEN a0_den, a1_den;
                DEN common_den;
                decompose(ext.a0(),a0_num,a0_den);
                decompose(ext.a1(),a1_num,a1_den);
                common_den=common_factor(a0_den,a1_den);
                
                a0_num = a0_num*CGAL::integral_division(a1_den,common_den);
                a1_num = a1_num*CGAL::integral_division(a0_den,common_den);
                den = CGAL::integral_division(a0_den,common_den)*a1_den;
                num = Numerator_type(a0_num,a1_num,ext.root());
            }else{
                NUM a0_num;
                decompose(ext.a0(),a0_num,den);
                num = Numerator_type(a0_num);
            }
        }
    };
    class Compose {
    public:
        typedef Numerator_type first_argument_type;
        typedef Denominator_type second_argument_type;
        typedef NT result_type;
        NT operator () (const Numerator_type&   num,
                        const Denominator_type& den){
            if(num.is_extended()){
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                COEFF a1=compose(num.a1(),den);
                return NT(a0,a1,num.root());
            }else{
                typename CFT::Compose compose;
                COEFF a0=compose(num.a0(),den);
                return NT(a0);
            }
        }
    };
};

template <class EXT, class BOOL_TAG>
class Sqrt_ext_Ftr_base_1;

template <class COEFF, class ROOT_NT>
class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT_NT >, CGAL::Tag_true >
    : public Sqrt_ext_Ftr_base_2< 
    Sqrt_extension<COEFF,ROOT_NT >, 
    ::boost::is_same< typename CGAL::Coercion_traits<ROOT_NT,typename CGAL::Fraction_traits<COEFF>::Numerator_type>::Type,
                        typename CGAL::Fraction_traits<COEFF>::Numerator_type>::value >
{
    //nothing new
};

    template <class COEFF, class ROOT_NT>
    class Sqrt_ext_Ftr_base_1< Sqrt_extension<COEFF,ROOT_NT>, CGAL::Tag_false >
        : public Sqrt_ext_Ftr_base_2< Sqrt_extension<COEFF,ROOT_NT >, false>
    {
        //nothing new
    };
} // namespace Intern


/*
namespace Intern{
    template <class SqrtExt,class BoolTag> class Sqrt_ext_Coftr_base_1;
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_false >{
    public:
        typedef SqrtExt          Numerator_type;
        typedef ::CGAL::Tag_false Is_composable;
        typedef ::CGAL::Null_tag Denominator_type;
        typedef ::CGAL::Null_tag Type;
        typedef ::CGAL::Null_tag Compose;  
    }; 
    template <class SqrtExt>
    class Sqrt_ext_Coftr_base_1< SqrtExt, CGAL::Tag_true >{
        typedef typename SqrtExt::NT Coeff;
        typedef typename SqrtExt::ROOT Root;
        typedef typename CGAL::Cofraction_traits<Coeff> CFT;
        typedef typename CFT::Type Type_coeff;
       
    public:
        typedef SqrtExt                                       Numerator_type;
        typedef ::CGAL::Tag_true                               Is_composable;
        typedef typename CFT::Denominator_type                Denominator;
        typedef CGAL::Sqrt_extension<Type_coeff,Root> Type;
        
        class Compose {
    public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type 
            typedef Denominator_type second_argument_type; 
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type   num, 
                                      Denominator_type den){
                if(num.is_extended()){
                    typename CFT::Compose compose_coeff;
                    Type_coeff a0_new(compose_coeff(num.a0(),den));
                    Type_coeff a1_new(compose_coeff(num.a1(),den));
                    return result_type(a0_new, a1_new, num.root());
                }else{
                    typename CFT::Compose compose_coeff;
                    return result_type(compose_coeff(num.a0(),den));                    
                }
            };
        };     
    };
}

template <class Coeff, class Root>
class Cofraction_traits<Sqrt_extension<Coeff,Root> >
    :public Intern::Sqrt_ext_Coftr_base_1< 
    Sqrt_extension<Coeff,Root>, 
    typename CGAL::Cofraction_traits<Coeff>::Is_composable>{
    //nothing new;
};
*/


template <class COEFF, class ROOT>
class Needs_parens_as_product< Sqrt_extension<COEFF,ROOT> >{
public:
    typedef Sqrt_extension<COEFF,ROOT> NT;
    bool operator()(const NT& t){         
        if( t.a0() != NT(0) && t.a1() != NT(0)){
            return true;
        } 
        if( t.a1() == NT(0) ){
            Needs_parens_as_product<COEFF> npap;
            return npap(t.a0());
        }
        return false;
    }
};




        
/////////// COERCION_TRAITS BEGIN 

// <EXT,EXT>
template <class A_coeff, class B_coeff, class Root>
struct Coercion_traits_for_level<Sqrt_extension<A_coeff, Root>, 
                           Sqrt_extension<B_coeff, Root>,
                           CTL_SQRT_EXT>{
private:
    typedef Coercion_traits<A_coeff, B_coeff> CT;
    typedef Sqrt_extension<A_coeff,Root> A;
    typedef Sqrt_extension<B_coeff,Root> B;
            
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type, Root> Type;   

    struct Cast{  
    private:
        inline Type cast(const Type& x) const{ return x; }
        
        template <class T>
        inline Type cast(const T& x) const{ 
            typename CT::Cast cast;
            if (x.is_extended()) {
                return result_type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return result_type(cast(x.a0()));
            }
        }
    public:
        typedef Type result_type;       
        // this is in order to allow A and B only 
        Type operator()(const A& x) const { return cast(x);}
        Type operator()(const B& x) const { return cast(x);}      
    };
};

template <class Coeff, class Root_1, class Root_2>
struct Coercion_traits_for_level<Sqrt_extension<Sqrt_extension<Coeff,Root_1>, 
Root_2>, 
                           Sqrt_extension<Coeff,Root_1>,
                           CTL_SQRT_EXT>{
private:
    typedef Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2> A;
    typedef Sqrt_extension<Coeff,Root_1> B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    // Type = A 
    typedef Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const A& x) const { return x;}
        Type operator()(const B& x) const { return Type(x);} 
    };  
};
        
template <class Coeff, class Root_1, class Root_2>
struct Coercion_traits_for_level
<
            Sqrt_extension<Coeff,Root_1>,
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2> 
            ,CTL_SQRT_EXT> 
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_2>, 
            Sqrt_extension<Coeff,Root_1> 
            ,CTL_SQRT_EXT>
{}; 

template <class Coeff, class Root_1>
struct Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1>, 
            Sqrt_extension<Coeff,Root_1> 
            ,CTL_SQRT_EXT>{
private:
    typedef  Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1> A;
    typedef  Sqrt_extension<Coeff,Root_1> B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef  Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const A& x) const { return x;} 
        Type operator()(const B& x) const { return Type(x);} 
    }; 
};
        
template <class Coeff, class Root_1>
struct Coercion_traits_for_level
<
            Sqrt_extension<Coeff,Root_1>,
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1> 
            ,CTL_SQRT_EXT> 
    :public Coercion_traits_for_level
<
            Sqrt_extension<Sqrt_extension<Coeff,Root_1>, Root_1>, 
            Sqrt_extension<Coeff,Root_1> 
            ,CTL_SQRT_EXT>
{};
        
        
namespace INTERN_CT{
// Coercion_traits for Sqrt_extenison to FieldWithSqrt
template <class A, class B> class CT_ext_to_fwsqrt;
// Coercion_traits for Sqrt_extenison not with FieldWithSqrt
template <class A, class B> class CT_ext_not_to_fwsqrt;
} // namespace INTERN_CT

        
//<EXT,ANY>
template <class Coeff, class Root, class B>
struct Coercion_traits_for_level<Sqrt_extension<Coeff, Root>, B , CTL_SQRT_EXT> 
:public ::boost::mpl::if_c< 
             // if B is fwsqrt
              ::boost::is_base_and_derived< 
                  Field_with_sqrt_tag, 
typename Algebraic_structure_traits<B>::Algebraic_category >::value || 
              ::boost::is_same< 
                  Field_with_sqrt_tag, 
typename Algebraic_structure_traits<B>::Algebraic_category >::value
            ,
            //then take Intern::Coercion_traits for fwsqrt
            INTERN_CT::CT_ext_to_fwsqrt<Sqrt_extension<Coeff,Root>, B>
            ,
            //else take Intern::Coercion_traits not for fwsqrt
            INTERN_CT::CT_ext_not_to_fwsqrt< Sqrt_extension<Coeff,Root> ,B>
              >::type
{};
        
// <ANY,EXT>
template <class Coeff, class Root, class B>
struct Coercion_traits_for_level
<B,Sqrt_extension<Coeff, Root>,CTL_SQRT_EXT >
    :public Coercion_traits_for_level<Sqrt_extension<Coeff,Root>,B,CTL_SQRT_EXT>
{};
        
namespace INTERN_CT{
// EXT coercion with FieldWithSqrt
template <class Coeff, class Root, class FieldWithSqrt>
struct CT_ext_to_fwsqrt<Sqrt_extension<Coeff,Root>,
                                         FieldWithSqrt>{
private:
    typedef Sqrt_extension<Coeff,Root> A;
    typedef FieldWithSqrt B;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef FieldWithSqrt Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const A& x) const { 
            typedef Coercion_traits<Coeff,FieldWithSqrt> CT_coeff;
            typedef Coercion_traits<Root ,FieldWithSqrt> CT_root;
            typename CT_coeff::Cast coeff_cast;
            typename CT_root::Cast root_cast;
            if (x.is_extended()) {   
                typename CGAL::Algebraic_structure_traits<
                typename CT_root::Type>::Sqrt sqrt;
                return // a0+a1*sqrt(root) 
                    coeff_cast(x.a0())+
                    coeff_cast(x.a1())*
                    sqrt(root_cast(x.root()));
            } else {
                return coeff_cast(x.a0());
            }
        } 
        Type operator()(const B& x) const { return x;} 
    }; 
};
    
// EXT coercion not with FieldWithSqrt
template <class Coeff, class Root, class B_>
struct CT_ext_not_to_fwsqrt<Sqrt_extension<Coeff,Root>, B_>{
private:
    typedef Sqrt_extension<Coeff,Root> A;
    typedef B_ B;
    typedef Coercion_traits<Coeff,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Sqrt_extension<typename CT::Type,Root> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const A& x) const { 
            typename CT::Cast cast;
            if (x.is_extended()) {
                return Type(cast(x.a0()),cast(x.a1()),x.root());
            } else {
                return Type(cast(x.a0()));
            }
        } 
        Type operator()(const B& x) const { 
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };   
};
} // namespace INTERN_CT

/////////// COERCION_TRAITS END 

#if 1 // Algebraic_extension_traits not supported yet 
/////////// ALGEBRAIC_NUMER_TRAITS BEGIN
template <class COEFF, class ROOT >
class Algebraic_extension_traits<Sqrt_extension<COEFF,ROOT> > {
// needed to 'add up' sqrt_extensions in iterator range such that all roots are 
//   collected in order to keep operation time minimal all scalar coeffs are set 
//   to 1 by standardise. 
//   TODO .. find a better name, if you want to.
//
template <class NT_>
class Standardise {
public:
    typedef NT_ argument_type;
    typedef NT_ result_type;
    NT_ operator () (const NT_& a) const {
        return NT_(1);
    }
};
    
template <class COEFF_, class ROOT_>
class Standardise<Sqrt_extension<COEFF_,ROOT_> > {
    Standardise<COEFF_> standardise;
public:
    typedef Sqrt_extension<COEFF_,ROOT_> NT_;
    typedef NT_ argument_type;
    typedef NT_ result_type;
    NT_ operator () (const NT_& a) const {       
        if (a.a1() != COEFF_(0)){
            return NT_(standardise(a.a0()),standardise(a.a1()),a.root());
        }else{
            return NT_(standardise(a.a0()));
        }
    }
};

public:
    //! \name Typedefs 
    //! the number type for which this instance has been instantiated
    typedef Sqrt_extension<COEFF,ROOT> NT;
    //! Sqrt_extension as a number type is extended
    typedef ::CGAL::Tag_true Is_extended;
    
    //! computes the factor which normalizes a number to be integral 
    //!  after multiplication
    //!
    class Normalization_factor{
    public:
        //! argument type
        typedef NT argument_type;
        //! result type
        typedef NT result_type;
        //! determine normalization factor
        NT operator () (const NT& a) const {
            typedef Algebraic_extension_traits<COEFF> SET;
            typename  SET::Normalization_factor normalization_factor;
            CGAL_precondition(a != NT(0));

            NT result;
            if(a.is_extended() && CGAL::sign(a.a1())!= CGAL::ZERO){
                NT tmp1(a.a0(),-a.a1(),a.root());
                
                NT tmp2= a*tmp1;
                CGAL_postcondition(tmp2.a1()==COEFF(0));
                result = tmp1*NT(normalization_factor(tmp2.a0()));
                CGAL_postcondition(CGAL::sign(result.a1()) != CGAL::ZERO);
            }else{
                result = NT(normalization_factor(a.a0())); 
            }       
            return result;
        }
    };

    //! returns the extension factor needed for the gcd_utcf computation 
    //! for more details see ... TODO!!
    //!
    class Denominator_for_algebraic_integers {
    public:
        //! argument type
        typedef NT argument_type;
        //! result type
        typedef NT result_type;
        //! determine denominator for algebraic integers
    public:
        NT operator () (const NT& a) const {
            typedef Algebraic_extension_traits<COEFF> ANT;
            typename ANT::Denominator_for_algebraic_integers dfai;

            Standardise<COEFF> standardise;
            if (a.a1() != COEFF(0)) {
                COEFF tmp = 
                    standardise(a.a0())
                    + standardise(a.a1())
                    + standardise(COEFF(a.root()));
                return  NT(COEFF(4) * COEFF(a.root()))* NT(dfai(tmp));
            } else {
                return NT(dfai(a.a0()));
            };
        }

        //! overloaded operator for computing the denominator out of an iterator
        //!  range accumulates all root expressions which appear in the range to 
        //!  the most complex term and uses this term to determine the denominator 
        //!  for algebraic integers
        //!
        template <class InputIterator>
        NT operator () (InputIterator begin, InputIterator end) const {
            NT a = std::accumulate(::boost::make_transform_iterator(begin,Standardise<NT>()), 
                                   ::boost::make_transform_iterator(end  ,Standardise<NT>()), NT(0));
            return (*this)(a);
        }
    };
};
#endif
/////////// ALGEBRAIC_NUMER_TRAITS BEGIN

CGAL_END_NAMESPACE

#endif  // CGAL_SQRT_EXTENSION_H

// EOF
