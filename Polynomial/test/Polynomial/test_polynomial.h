
#include <CGAL/basic.h>
#include <CGAL/use.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

#include <iostream>
#include <cassert>

#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Coercion_traits.h>


#define CGAL_DEFINE_TYPES_FROM_AK(AK)                           \
    typedef typename AK::Integer Integer;                       \
    typedef typename AK::Rational Rational;                     \
    typedef typename AK::Field_with_sqrt Field_with_sqrt;  CGAL_USE_TYPE(Field_with_sqrt); \
    typedef CGAL::Polynomial<Integer> Poly_int1;                \
    typedef CGAL::Polynomial<Poly_int1> Poly_int2;              \
    typedef CGAL::Polynomial<Poly_int2> Poly_int3; CGAL_USE_TYPE(Poly_int3); \
    typedef CGAL::Polynomial<Rational> Poly_rat1;               \
    typedef CGAL::Polynomial<Poly_rat1> Poly_rat2;              \
    typedef CGAL::Polynomial<Poly_rat2> Poly_rat3;  CGAL_USE_TYPE(Poly_rat3);
   

// TODO: copied from number_type_utils.h
namespace CGAL {
template <class NT , class RT>
inline 
void convert_to(const NT& x, RT& r){ 
    typedef CGAL::Coercion_traits<NT,RT> CT;
    CGAL_assertion_code(typedef typename CT::Coercion_type RET;)
    CGAL_static_assertion((::boost::is_same<RET,RT>::value));
    r = typename CT::Cast()(x);
}
} //namespace CGAL

template <class NT>
void datastru() {
    typedef CGAL::Polynomial<NT> POLY;
    int i;

    // uninitialized construction
    POLY p;
    assert( p.degree() == 0 );

    // construction from iterator range, assignment
    NT array[] = { NT(3), NT(1), NT(4), NT(1), NT(5), NT(9) };
    NT* const array_end = array + sizeof array / sizeof *array;
    POLY q(array, array_end);
    p = q;
    assert( array_end - array == p.degree() + 1 );
    assert( std::equal(p.begin(), p.end(), array) );

    // immediate construction, reduce(), copy construction, (in)equality
    POLY r(NT(3), NT(1), NT(4), NT(1), NT(5), NT(9), NT(0), NT(0));
    assert( q.degree() == r.degree() );
    assert( q == r );
    assert( r == r );
    q = POLY(NT(2), NT(7), NT(1), NT(8), NT(1));
    assert( q != r );
    POLY s(p);
    assert( r == s );

    // selection of coefficients
    for (i = 0; i <= p.degree(); i++) {
        assert( p[i] == array[i] );
    }
    assert( i == array_end - array );
    assert( p.lcoeff() == array_end[-1] );

    // reversal()
    assert( reversal(q) == POLY(NT(1), NT(8), NT(1), NT(7), NT(2)) );

    // divide_by_x()
    r.divide_by_x();
    assert( r.degree() == 4 );
    assert( r == POLY(NT(1), NT(4), NT(1), NT(5), NT(9)) );

    // zero test
    assert( !p.is_zero() );
    q = POLY(NT(0));
    assert( q.is_zero() );
}

template <class NT> void signs(::CGAL::Tag_false);
template <class NT> void signs(::CGAL::Tag_true);
template <class NT> void division(CGAL::Integral_domain_without_division_tag);
template <class NT> void division(CGAL::Integral_domain_tag);

template <class NT>
void arithmetic() {
    typedef CGAL::Polynomial<NT> POLY;
    // evaluation
    POLY p(NT(6), NT(-5), NT(1));
    assert( p.evaluate(NT(1)) == NT(2) );
    assert( p.evaluate(NT(2)) == NT(0) );
    assert( p.evaluate(NT(3)) == NT(0) );
    // TODO: No evaluate_absolute available for all NTs.
    //assert( p.evaluate_absolute(NT(1)) == NT(12) );

    // diff
    POLY q(NT(3), NT(2), NT(-2));
    POLY z(NT(0));
    assert( CGAL::differentiate(z) == z );
    assert( p.degree() == 2 );
    assert( CGAL::differentiate(p) == POLY(NT(-5), NT(2)) );
    assert( p.degree() == 2 );
    assert( q.degree() == 2 );
    q.diff();
    assert( q == POLY(NT(2), NT(-4)) );

    // scale and translate
    p.scale_up(NT(4));
    assert( p == POLY(NT(6), NT(-20), NT(16)) );
    assert( scale_down(p, NT(2)) == POLY(NT(24), NT(-40), NT(16)) );
    assert( p == POLY(NT(6), NT(-20), NT(16)) );
    p = POLY(NT(3), NT(2), NT(1));
    // CGAL::translate takes an Innermost_coefficient 
    assert( translate_by_one(p) == CGAL::translate(p, 1) );
    p.translate(NT(-2));
    assert( p == POLY(NT(3), NT(-2), NT(1)) );

    // addition
    POLY r;
    p = POLY(NT(3), NT(2), NT(1), NT(1));
    q = POLY(NT(1), NT(2), NT(-1), NT(-1));
    r = p + q;
    assert( r.degree() == 1 );
    assert( r[0] == NT(4) );
    assert( r[1] == NT(4) );
    q += p;
    assert( q == r );
    q += NT(1);
    assert( q == POLY(NT(5), NT(4)) );
    assert( q == r + NT(1) );
    assert( q == NT(1) + r );

    // subtraction
    p = POLY(NT(3), NT(2), NT(1), NT(1));
    q = POLY(NT(1), NT(2), NT(-1), NT(1));
    r = p - q;
    assert( r.degree() == 2 );
    assert( r == POLY(NT(2), NT(0), NT(2)) );
    r -= -q;
    assert( r == p );
    r -= NT(2);
    assert(  r == p - NT(2) );
    assert( -r == NT(2) - p );

    // multiplication
    p = POLY(NT(1), NT(2), NT(3));
    q = POLY(NT(-2), NT(5), NT(2));
    r = p * q;
    assert( r == POLY(NT(-2), NT(1), NT(6), NT(19), NT(6)) );
    q *= p;
    assert( q == r );
    q *= NT(2);
    assert( q == POLY(NT(-4), NT(2), NT(12), NT(38), NT(12)) );
    assert( q == r * NT(2) );
    assert( q == NT(2) * r );

    // sign etc., if applicable
    // TODO: Replaced Is_real_comparable with RET::Is_real_embeddable, OK?
    typedef typename CGAL::Real_embeddable_traits<NT>::Is_real_embeddable 
      Embeddable;
    signs<NT>(Embeddable());

    // division, if applicable
    typedef typename CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
      Algebra_type;
    division<NT>(Algebra_type());
}

template <class NT> void signs(::CGAL::Tag_false) {
    // nothing to be done
}

template <class POLY>
void test_greater(POLY a, POLY b) {
    // check .compare()
    assert( a.compare(b) == CGAL::LARGER );
    assert( a.compare(a) == CGAL::EQUAL );
    assert( b.compare(a) == CGAL::SMALLER );

    // check that all comparison operators reflect a > b
    assert(   a >  b  );
    assert( !(b >  a) );
    assert( !(a >  a) );
    assert(   b <  a  );
    assert( !(a <  b) );
    assert( !(a <  a) );
    assert(   a >= b  );
    assert( !(b >= a) );
    assert(   a >= a  );
    assert(   b <= a  );
    assert( !(a <= b) );
    assert(   a <= a  );
}

template <class NT>
void signs(::CGAL::Tag_true) {
    typedef CGAL::Polynomial<NT> POLY;

    POLY p(NT(6), NT(-5), NT(1));
    POLY pp(NT(6), NT(-5), NT(1));
    POLY q(NT(3), NT(2), NT(-2));
    POLY z(NT(0));

    // sign, abs
    assert( z.sign() == CGAL::ZERO );
    assert( p.sign() == CGAL::POSITIVE );
    assert( p.abs() == p );
    assert( q.sign() == CGAL::NEGATIVE );
    assert( q.abs() == POLY(NT(-3), NT(-2), NT(2)) );

    // comparison operators
    test_greater(p + NT(1), p);
    test_greater(POLY(NT(1), NT(2), NT(3), NT(4)), p);
    
    // sign_at 
    q=POLY(NT(3),NT(2),NT(-1));
    assert(q.sign_at(-2)==CGAL::NEGATIVE);
    assert(q.sign_at(-1)==CGAL::ZERO);
    assert(q.sign_at( 0)==CGAL::POSITIVE);
    assert(q.sign_at( 3)==CGAL::ZERO);
    assert(q.sign_at( 4)==CGAL::NEGATIVE);

    assert(z.sign_at( 3)==CGAL::ZERO);
}

template <class NT>
void division(CGAL::Integral_domain_without_division_tag) {
    // nothing to be done
}

template <class NT>
void division(CGAL::Integral_domain_tag) {
    typedef CGAL::Polynomial<NT> POLY;

    // integral division (remainder zero)
    POLY p(NT(2), NT(4), NT(6)), q(NT(-2), NT(5), NT(2));
    POLY r(NT(-4), NT(2), NT(12), NT(38), NT(12)), s;
    assert( q == r / p );
    r /= q;
    assert( r == p );
    r /= NT(2);
    assert( r == POLY(NT(1), NT(2), NT(3)) );
    assert( NT(6) / POLY(NT(3)) == POLY(NT(2)) );

    // division with remainder
    p = POLY(NT(2), NT(-3), NT(1));    //  (x-1)(x-2)
    q = POLY(NT(24), NT(-14), NT(2));  // 2(x-3)(x-4)
    r = POLY(NT(-5), NT(1));           //       (x-5)
    s = p*q + r;

    assert(!CGAL::divides(q,s));

    POLY Q, R;
    POLY::euclidean_division(s, p, Q, R);  // works, since p.lcoeff() == 1
    assert( Q == q );
    assert( R == r );

    NT D;
    POLY::pseudo_division(s, q, Q, R, D);
    assert( Q == D*p );
    assert( R == D*r );

    // codecover 
    p = POLY(NT(2),NT(3));
    s = POLY(NT(2),NT(3),NT(4));
    POLY::euclidean_division(p, s, Q, R);  // works, since p.lcoeff() == 1
    assert( Q == POLY(NT(0)) );
    assert( R == p );

    POLY::pseudo_division(p, s, Q, R, D);  // works, since p.lcoeff() == 1
    assert( Q == POLY(NT(0)) );
    assert( R == D*p );
        
}

template <class NT>
void io() {
    typedef CGAL::Polynomial<NT> POLY;

    {
        // successful re-reading of output
        POLY p(NT(-3), NT(5), NT(0), NT(0), NT(-7), NT(9)), q;
        std::ostringstream os;
        os << p;
        std::istringstream is(os.str());
        is >> q;
        assert( p == q );
    }{       
        std::ostringstream os;
        CGAL::set_pretty_mode(os);
        os << oformat(POLY(NT(3)));
        //std::cout <<os.str()<<std::endl;
        assert( os.str() == "3" );
    }{       
        std::ostringstream os;
        CGAL::set_pretty_mode(os);
        os << oformat(POLY(NT(-3)));
        assert( os.str() == "(-3)" );
    }{       
        std::ostringstream os;
        CGAL::set_pretty_mode(os);
        os << oformat(POLY(NT(-3)),CGAL::Parens_as_product_tag());
        assert( os.str() == "(-3)" );
    }{       
        std::ostringstream os;
        CGAL::set_pretty_mode(os);
        os << oformat(POLY(NT(-3),NT(4)));
        if( CGAL::Polynomial_traits_d<POLY>::d == 1)
            assert( os.str() == "4*x + (-3)" );
        else
            assert( os.str() == "4*y + (-3)" );
    }{       
        std::ostringstream os;
        CGAL::set_pretty_mode(os);
        os << oformat(POLY(NT(-3),NT(4)), CGAL::Parens_as_product_tag());
        
        if( CGAL::Polynomial_traits_d<POLY>::d == 1)
            assert( os.str() == "(4*x + (-3))" );
        else
            assert( os.str() == "(4*y + (-3))" ); 
        
    }
    
}

template <class NT>
void canon(CGAL::Integral_domain_without_division_tag) {
    // dummy for cases where CGAL::canonicalize() does not apply
}

template <class NT>
void canon(CGAL::Unique_factorization_domain_tag) {
    typedef CGAL::Polynomial<NT>     UNPOLY;
    typedef CGAL::Polynomial<UNPOLY> BIPOLY;

    BIPOLY p(0), c;
    assert(CGAL::canonicalize(p) == BIPOLY(0));

    p = BIPOLY(UNPOLY(NT(15)), UNPOLY(NT(-12), NT( 9)), UNPOLY(NT(0), NT(-6)));
    c = BIPOLY(UNPOLY(NT(-5)), UNPOLY(NT(  4), NT(-3)), UNPOLY(NT(0), NT( 2)));
    assert(CGAL::canonicalize(p) == c);
}

template <class NT>
void canon(CGAL::Field_tag) {
    typedef CGAL::Polynomial<NT>     UNPOLY;
    typedef CGAL::Polynomial<UNPOLY> BIPOLY;

    BIPOLY p(0), c;
    assert(CGAL::canonicalize(p) == BIPOLY(0));

    p = BIPOLY(UNPOLY(NT(15)), UNPOLY(NT(-12), NT( 9)), UNPOLY(NT(0), NT(-6)));
    c = BIPOLY(UNPOLY(NT(-5)/NT(2)), UNPOLY(NT(2), NT(-3)/NT(2)),
            UNPOLY(NT(0), NT(1)));
    assert(CGAL::canonicalize(p) == c);
}

#ifdef NiX_POLY_USE_NT_TESTS
template <class NT>
void nt_tests(::CGAL::Tag_false) {
    typedef CGAL::Polynomial<NT> POLY;
    if ((volatile bool)0) { // test compilation only
        CGAL::test_number_type<POLY,
                              typename CGAL::Algebraic_structure_traits<POLY>::Algebric_structure_tag>();
    }
}

template <class NT>
void nt_tests(::CGAL::Tag_true) {
    typedef CGAL::Polynomial<NT> POLY;
    CGAL::test_real_comparable<POLY>();
    nt_tests<NT>(::CGAL::Tag_false());
}
#endif // NiX_POLY_USE_NT_TESTS

template <class NT>
void unigcdres(CGAL::Field_tag) {
    typedef CGAL::Polynomial<NT> POLY;

    // test univariate polynomial gcd and resultant over a field
    POLY f(NT(2), NT(7), NT(1), NT(8), NT(1), NT(8));
    POLY g(NT(3), NT(1), NT(4), NT(1), NT(5), NT(9));
    POLY h(NT(3), NT(4), NT(7), NT(7));
    f /= NT(3);
    g /= NT(5);
    h /= h.unit_part();
    POLY d;

    d = gcd(f, g);
    assert( d == POLY(NT(1)) );
    assert( prs_resultant(f, g) == NT(230664271L)/NT(759375L) ); // Maple

    POLY fh(f*h), gh(g*h);
    d = gcd(fh, gh);
    assert( d == h );
    assert( prs_resultant(fh, gh) == NT(0) );

    POLY a, b;
    d = gcdex(fh, gh, a, b);
    assert( d == h );
    assert( d == a*fh + b*gh );
}

template <class NT>
void unigcdres(CGAL::Integral_domain_tag) {
    typedef CGAL::Polynomial<NT> POLY;

    // test univariate polynomial gcd and resultant over a non-field
    POLY f(NT(2), NT(7), NT(1), NT(8), NT(1), NT(8));
    POLY g(NT(3), NT(1), NT(4), NT(1), NT(5), NT(9));
    POLY h(NT(3), NT(4), NT(7), NT(7));
    POLY d;
    NT c(42);

    d = gcd((-c)*f, c*g);
    assert( d == POLY(c) );
    assert( prs_resultant(f, g) == NT(230664271L) ); // as computed by Maple

    POLY fh(f*h), gh(g*h);
    d = gcd((-c)*fh, c*gh);
    assert( d == c*h );
    assert( prs_resultant(fh, gh) == NT(0) );

    POLY a, b; NT v;
    d = pseudo_gcdex((-c)*fh, c*gh, a, b, v);
    assert( d == c*h );
    assert( v*d == (-c)*a*fh + c*b*gh );

    // Michael Kerber's example for the hgdelta_update() bug:
    // These polynomials cretate a situation where h does not divide g,
    // but h^(delta-1) divides g^delta (as predicted by subresultant theory).
    // TODO: Uncomment following code
    /*CGAL::Creator_1<int, NT> int2nt;
    static const int cp[] = { 196, 0, -140, 0, 69, 0, -8, 0, -4 },
                     cq[] = { -15156, 0, -58572, 0, -82309, 0, -27032, 0,
                         20190, 0, 1696, 0, -6697, 0, 1644, 0, 1664, 0, 256 };
    static const size_t np = sizeof cp / sizeof *cp,
                        nq = sizeof cq / sizeof *cq;
    POLY p(LiS::mapper(cp, int2nt), LiS::mapper(cp+np, int2nt)),
         q(LiS::mapper(cq, int2nt), LiS::mapper(cq+nq, int2nt));
    d = gcd(p,q);
    assert ( d == POLY(1) );*/
}

template <class NT>
void bigcdres(CGAL::Field_tag) {
    typedef CGAL::Polynomial<NT>    POLY1;
    typedef CGAL::Polynomial<POLY1> POLY2;

    // test bivariate polynomial gcd and resultant over a field
    POLY2 f = POLY2(POLY1(NT(2), NT(7), NT(1)),
                    POLY1(NT(8), NT(1)), POLY1(NT(8)));
    POLY2 g = POLY2(POLY1(NT(3), NT(1), NT(4)),
                    POLY1(NT(1), NT(5)), POLY1(NT(9)));
    POLY2 h = POLY2(POLY1(NT(3), NT(4), NT(1)),
                    POLY1(NT(9), NT(7)), POLY1(NT(7)));
    f /= NT(3);
    g /= NT(5);
    h /= h.unit_part();
    POLY2 d;

    d = gcd(f, g);
    assert( d == POLY2(1) );
    POLY1 r(NT(1444), NT(-1726), NT(3295), NT(-2501), NT(560));
    r /= NT(225);
    assert( prs_resultant(f, g) == r ); // says Maple

    POLY2 fh(f*h), gh(g*h);
    d = gcd(fh, gh);
    assert( d == h );
    assert( prs_resultant(fh, gh) == POLY1(0) );
}

template <class NT>
void bigcdres(CGAL::Integral_domain_tag) {
    typedef CGAL::Polynomial<NT>    POLY1;
    typedef CGAL::Polynomial<POLY1> POLY2;

    // test bivariate polynomial gcd and resultant over a non-field
    POLY2 f = POLY2(POLY1(NT(2), NT(7), NT(1)),
                    POLY1(NT(8), NT(1)), POLY1(NT(8)));
    POLY2 g = POLY2(POLY1(NT(3), NT(1), NT(4)),
                    POLY1(NT(1), NT(5)), POLY1(NT(9)));
    POLY2 h = POLY2(POLY1(NT(3), NT(4), NT(1)),
                    POLY1(NT(9), NT(7)), POLY1(NT(7)));
    POLY2 c(42), d;

    d = gcd(-c*f, c*g);
    assert( d == c );
    POLY1 r(NT(1444), NT(-1726), NT(3295), NT(-2501), NT(560));
    assert( prs_resultant(f, g) == r ); // says Maple

    POLY2 fh(f*h), gh(g*h);
    d = gcd(-c*fh, c*gh);
    assert( d == c*h );
    assert( prs_resultant(fh, gh) == POLY1(0) );
}

template <class POLY>
void test_sqff_utcf_(const POLY& poly, int n){
    typedef CGAL::Polynomial_traits_d<POLY> PT;

    std::vector<POLY> fac;
    std::vector<int>  mul;
    std::back_insert_iterator<std::vector<POLY> > fac_bi(fac);
    std::back_insert_iterator<std::vector<int>  > mul_bi(mul);
    
    int tmp = CGAL::internal::square_free_factorize_utcf(poly, fac_bi, mul_bi);

    assert(n == tmp);

    assert((int) mul.size() == n);
    assert((int) fac.size() == n);
    if (typename PT::Total_degree()(poly) == 0 ){
        assert(n == 0); 
    }else{
        POLY check(1);
        for(int i = 0; i < n ; i++){
            assert(typename PT::Total_degree()(fac[i]) > 0);
            //std::cout << fac[i] << std::endl;
            //std::cout << typename PT::Canonicalize()(fac[i]) << std::endl;
            //std::cout << std::endl;
            assert(fac[i] == typename PT::Canonicalize()(fac[i]));
            check *= CGAL::ipower(fac[i],mul[i]);
        }
        assert(check == typename PT::Canonicalize()(poly));
    }
}

template <class POLY_1>
void test_sqff_utcf(const POLY_1& a1,const POLY_1& b1,const POLY_1& c1){
    typedef CGAL::Polynomial<POLY_1>  POLY_2; 
 
    test_sqff_utcf_(POLY_1(0), 0);
    test_sqff_utcf_(POLY_1(1), 0);
    test_sqff_utcf_(POLY_1(2), 0);
    test_sqff_utcf_(a1, 1);
    test_sqff_utcf_(a1*a1, 1);
    test_sqff_utcf_(a1*a1*b1,2);
    test_sqff_utcf_(a1*POLY_1(5),1);

    POLY_2 a2(a1,b1,c1);
    POLY_2 b2(a1,c1);
    POLY_2 c2(b1,a1,b1);
 
    test_sqff_utcf_(POLY_2(0), 0);
    test_sqff_utcf_(POLY_2(1), 0);
    test_sqff_utcf_(POLY_2(2), 0);

    test_sqff_utcf_(a2,1);
    test_sqff_utcf_(a2*a2,1);
    test_sqff_utcf_(a2*a2*b2,2);
    test_sqff_utcf_(a2*POLY_2(5),1);
    test_sqff_utcf_(a2*a2*POLY_2(5),1);
    
    // non regular polynomials
    test_sqff_utcf_(a1*a2,2);
    test_sqff_utcf_(a1*a2*a2,2);
    test_sqff_utcf_(a1*a2*b2*b2,3);
    test_sqff_utcf_(a1*a2*POLY_2(5),2);
    test_sqff_utcf_(a1*a2*a2*POLY_2(5),2);
    test_sqff_utcf_(a1*a2*a2*b2*POLY_2(5),3);
    test_sqff_utcf_(a1*a1*a2,2);
    test_sqff_utcf_(a1*a1*a2*a2,2);
    test_sqff_utcf_(a1*b1*b1*a2*a2,3);

}



template <class AT>
void psqff(){
    {
        typedef typename AT::Integer NT;
        typedef CGAL::Polynomial<NT> POLY_1; 
        
        // monic factors 
        POLY_1 a1(NT(std::string("4352435")), NT(std::string("4325245")));
        POLY_1 b1(NT(std::string("123")), NT(std::string("432235")),NT(std::string("43324252")));
        POLY_1 c1(NT(std::string("12324")), NT(std::string("25332235")),NT(std::string("24657252")));
 
        test_sqff_utcf(a1,b1,c1);
    }{
        typedef typename AT::Rational NT;
        typedef CGAL::Polynomial<NT> POLY_1; 
        
        // monic factors 
        POLY_1 a1(NT(std::string("4352435")), NT(std::string("4325245")));
        POLY_1 b1(NT(std::string("123")), NT(std::string("432235")),NT(std::string("43324252")));
        POLY_1 c1(NT(std::string("12324")), NT(std::string("25332235")),NT(std::string("24657252")));
 
        test_sqff_utcf(a1,b1,c1);
    }
    {
       typedef typename AT::Integer Integer; //UFD domain
       typedef CGAL::Sqrt_extension<Integer,Integer> NT;
       typedef CGAL::Polynomial<NT> POLY;
       // square-free factorization (i.e. factorization by multiplicities)
       POLY a1(NT(3,4,7), NT(2,3,7), NT(3));
       POLY b1(NT(9,2,7), NT(5,11,7), NT(2));
       POLY c1(NT(2,8,7), NT(6,3,7), NT(5));
       test_sqff_utcf(a1,b1,c1);
    }{
       typedef typename AT::Integer Integer; //UFD domain
       typedef typename AT::Rational Rational; //UFD domain
       typedef CGAL::Sqrt_extension<Rational,Integer> NT;
       typedef CGAL::Polynomial<NT> POLY;
       // square-free factorization (i.e. factorization by multiplicities)
       POLY a1(NT(3,4,7), NT(2,3,7), NT(3));
       POLY b1(NT(9,2,7), NT(5,11,7), NT(2));
       POLY c1(NT(2,8,7), NT(6,3,7), NT(5));
       test_sqff_utcf(a1,b1,c1);
    }
    {
      typedef typename AT::Integer Integer; CGAL_USE_TYPE(Integer);
       typedef typename AT::Rational Rational;
       typedef CGAL::Sqrt_extension<Rational,Rational> NT;
       typedef CGAL::Polynomial<NT> POLY;
       // square-free factorization (i.e. factorization by multiplicities)
       POLY a1(NT(3,4,7), NT(2,3,7), NT(3));
       POLY b1(NT(9,2,7), NT(5,11,7), NT(2));
       POLY c1(NT(2,8,7), NT(6,3,7), NT(5));
       test_sqff_utcf(a1,b1,c1);
    }
}
template <class NT>
void sqff() {
    typedef CGAL::Polynomial<NT> POLY;

    // square-free factorization (i.e. factorization by multiplicities)
    POLY p1(NT(3), NT(4), NT(1));
    POLY p3(NT(9), NT(5), NT(1));
    POLY p4(NT(2), NT(6), NT(1));
    POLY p = NT(3) * p1 * p3*p3*p3 * p4*p4*p4*p4;

    typedef std::vector<POLY> PVEC;
    typedef std::vector<int>  IVEC;
    PVEC fac;
    IVEC mul;
    std::back_insert_iterator<PVEC> fac_bi(fac);
    std::back_insert_iterator<IVEC> mul_bi(mul);
    unsigned n;
    NT alpha;
    n = CGAL::internal::square_free_factorize(p, fac_bi, mul_bi );

//    assert(alpha == 3);
    assert(n == 3);
    assert(mul.size() == n);
    assert(fac.size() == n);
    assert(mul[0] == 1 && fac[0] == p1);
    assert(mul[1] == 3 && fac[1] == p3);
    assert(mul[2] == 4 && fac[2] == p4);
    
    /*p = POLY( NT(1), NT(-2), NT(1) );
    std::cerr << p << std::endl;
    fac.clear();
    mul.clear();
    fac_bi = std::back_insert_iterator<PVEC>(fac);
    mul_bi = std::back_insert_iterator<IVEC>(mul);
    std::cerr << CGAL::internal::square_free_factorize( p, fac_bi, mul_bi ) << std::endl;
    std::cerr << fac[0] << std::endl;*/
    //std::cerr << fac[1] << std::endl;
    
    typedef CGAL::Polynomial< CGAL::Polynomial< NT > > BPOLY;
    typedef std::vector< BPOLY > BPVEC;
    BPOLY bp( POLY( NT(1), NT(-2), NT(1) ) );
//    std::cerr << bp << std::endl;
//    typename CGAL::Polynomial_traits_d< BPOLY >::Make_square_free make_square_free;
//    std::cerr << make_square_free( bp ) << std::endl;
//    std::cerr << "^^^^^^^^^^" << std::endl;
    BPVEC bfac;
    mul.clear();
    std::back_insert_iterator< BPVEC > bfac_bi( bfac );
    mul_bi = std::back_insert_iterator<IVEC>(mul);
//    std::cerr << CGAL::POLYNOMIAL::square_free_factorize( bp, bfac_bi, mul_bi ) << std::endl;
//    std::cerr << bfac[0] << std::endl;
    
}

template <class FNT>
void integr() {
    typedef CGAL::Polynomial<FNT> FPOLY;
    typedef typename CGAL::Fraction_traits<FPOLY>::Numerator_type IPOLY;
    typedef typename CGAL::Fraction_traits<FPOLY>::Denominator_type DENOM;
    typename CGAL::Fraction_traits<FPOLY>::Decompose decompose;
    typename CGAL::Fraction_traits<FPOLY>::Compose compose; (void) compose;
    typedef typename IPOLY::NT INT;

    // making polynomials integral and fractional
    FPOLY fp(FNT(1)/FNT(2), FNT(2)/FNT(3), FNT(3)/FNT(5));
    IPOLY ip; DENOM d;
    decompose(fp, ip,d);
    assert( d == DENOM(30) );
    assert( ip == IPOLY(INT(15), INT(20), INT(18)) );
    assert( fp == compose(ip, d) );
}

template <class NT>
void basic_tests() {
    // tests that should work for all number types
    typedef CGAL::Polynomial<NT> POLY;
    datastru<NT>();
    datastru<POLY>();
    arithmetic<NT>();
    arithmetic<POLY>();
    io<NT>();
    io<POLY>();
    canon<NT>(typename CGAL::Algebraic_structure_traits<NT  >::Algebraic_category());
#ifdef NiX_POLY_USE_NT_TESTS
    // TODO: Replaced Is_real_comparable with RET::Is_real_embeddable, OK?
    nt_tests<NT  >(typename CGAL::Real_embeddable_traits<NT  >::Is_real_embeddable());
    // TODO: Replaced Is_real_comparable with RET::Is_real_embeddable, OK?
    nt_tests<POLY>(typename CGAL::Real_embeddable_traits<POLY>::Is_real_embeddable());
#endif // NiX_POLY_USE_NT_TESTS
}

template <class NT>
void exact_tests() {
    // tests that probably don't work for inexact or bounded number types
    unigcdres<NT>(typename CGAL::Algebraic_structure_traits<NT>::Algebraic_category());
    bigcdres<NT>(typename CGAL::Algebraic_structure_traits<NT>::Algebraic_category());
    sqff<NT>();
}

template < class Rational >
void test_coefficients_to() {
    typedef CGAL::Polynomial< int >       Poly_i1;
    typedef CGAL::Polynomial< Poly_i1 >   Poly_i2;
    typedef CGAL::Polynomial< Poly_i2 >   Poly_i3;
    
    typedef CGAL::Polynomial< Rational >  Poly_rat1;
    typedef CGAL::Polynomial< Poly_rat1 > Poly_rat2;
    typedef CGAL::Polynomial< Poly_rat2 > Poly_rat3;
    
    // univariate
    {
        Poly_i1 p = Poly_i1(-3,4,5,-9,0,2);
        
        Poly_rat1 r = Poly_rat1(Rational(-3),Rational(4),Rational(5),
                                Rational(-9),Rational(0),Rational(2));
        Poly_rat1 tmp;
        CGAL::convert_to(p,tmp); 
        assert(tmp == r);                
    }
    
    // bivariate
    {
        Poly_i1 p1 = Poly_i1( 3,-1, 2);
        Poly_i1 p2 = Poly_i1(-5, 2,-4);
        Poly_i2 p = Poly_i2(p1, p2);
        
        Poly_rat1 r1 = Poly_rat1(Rational( 3),Rational(-1),Rational( 2));
        Poly_rat1 r2 = Poly_rat1(Rational(-5),Rational( 2),Rational(-4));
        Poly_rat2 r = Poly_rat2(r1,r2);
        
        Poly_rat2 tmp;
        CGAL::convert_to(p,tmp); 
        assert(tmp == r);     
    }

    // trivariate
    {
        Poly_i1 p11 = Poly_i1( 7,-6, 9);
        Poly_i1 p12 = Poly_i1(-1, 3);
        Poly_i2 p1  = Poly_i2(p11, p12);
        
        Poly_i1 p21 = Poly_i1(-17,  16, -19);
        Poly_i1 p22 = Poly_i1( 11, -13);
        Poly_i2 p2  = Poly_i2(p21, p22);
        
        Poly_i3 p   = Poly_i3(p1,p2);
 
        Poly_rat1 r11 = Poly_rat1(Rational( 7),Rational(-6),Rational( 9));
        Poly_rat1 r12 = Poly_rat1(Rational(-1),Rational( 3));
        Poly_rat2 r1  = Poly_rat2(r11, r12);
        
        Poly_rat1 r21 = Poly_rat1(Rational(-17),Rational( 16),Rational(-19));
        Poly_rat1 r22 = Poly_rat1(Rational( 11),Rational(-13));
        Poly_rat2 r2  = Poly_rat2(r21, r22);
        
        Poly_rat3 r   = Poly_rat3(r1,r2);
        
        Poly_rat3 tmp;
        CGAL::convert_to(p,tmp); 
        assert(tmp == r);
    }
}

template <class AT>
void test_evaluate(){
    CGAL_DEFINE_TYPES_FROM_AK(AT);
    { 
        Poly_int1 P(3,0,2);
        Integer x(2);        
        assert(P.evaluate(x)==Integer(11));
    }{    
        Poly_rat1 P(3,0,2);
        Integer x(2);
        assert(P.evaluate(x)==Rational(11));
    }{    
        Poly_int1 P(3,0,2);
        Rational x(2);
        assert(P.evaluate(x)==Rational(11));
    }{ 
        Poly_int1 P(3,0,2);
        CGAL::Interval_nt<true> x(2);    
        assert(P.evaluate(x).inf()<=11);
        assert(P.evaluate(x).sup()>=11);
    }{
        CGAL::Polynomial< CGAL::Interval_nt<true> > P(3,0,2);
        CGAL::Interval_nt<true> x(2);
        assert(P.evaluate(x).inf()<=11);
        assert(P.evaluate(x).sup()>=11);
    }{
        CGAL::Polynomial< CGAL::Interval_nt<true> > P(3,0,2);
        Integer x(2);
        assert(P.evaluate(x).inf()<=11);
        assert(P.evaluate(x).sup()>=11);
        }
}

template <class AT>
void test_evaluate_homogeneous(){
    CGAL_DEFINE_TYPES_FROM_AK(AT);
    { 
        Poly_int1 P(3,0,2);
        Integer u(2);
        Integer v(1);
        assert(P.evaluate_homogeneous(u,v)==Integer(11));
    }{ 
        Poly_rat1 P(2,-1,3);
        Rational u(3);
        Rational v(5);
        assert(P.evaluate_homogeneous(u,v)==Rational(62));
        assert(P.evaluate(u/v)==Rational(62)/Rational(25));
    }
}
void test_total_degree(){
    typedef CGAL::Polynomial<int>    Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    
    assert(CGAL::total_degree(5) == 0);
    
    assert(CGAL::total_degree(Poly_1(0))      == Poly_1(0).degree());
    assert(CGAL::total_degree(Poly_1(1))      == Poly_1(1).degree());      
    assert(CGAL::total_degree(Poly_1(1,1))    == Poly_1(1,1).degree());
    assert(CGAL::total_degree(Poly_1(0,0,0,1))== Poly_1(0,0,0,1).degree());
    
    assert(CGAL::total_degree(Poly_2(0))      ==  0);
    assert(CGAL::total_degree(Poly_2(1))      ==  0);
    assert(CGAL::total_degree(Poly_2(Poly_1(1),Poly_1(1))) ==  1);
    assert(CGAL::total_degree(Poly_2(Poly_1(0),Poly_1(0),Poly_1(0),Poly_1(1)))==  3);
    
    assert(CGAL::total_degree(Poly_2(Poly_1(1,1),Poly_1(1)))  == 1);
    assert(CGAL::total_degree(Poly_2(Poly_1(0),Poly_1(0),Poly_1(1)))==  2);
    assert(CGAL::total_degree(Poly_2(Poly_1(0),Poly_1(0),Poly_1(1)))== 2);
    assert(CGAL::total_degree(Poly_2(Poly_1(0),Poly_1(0),Poly_1(1,1)))== 3);
    assert(CGAL::total_degree(Poly_2(Poly_1(1),Poly_1(0),Poly_1(1,1)))== 3);
    assert(CGAL::total_degree(Poly_2(Poly_1(1),Poly_1(1),Poly_1(1,1)))== 3);
    assert(CGAL::total_degree(Poly_2(Poly_1(1),Poly_1(1,1,1,1),Poly_1(1,1)))== 4);
}

template<class AT>
void test_scalar_factor_traits(){
    
    {   
        typedef typename AT::Integer Integer;
        typedef CGAL::Polynomial<Integer> Polynomial; 
        typedef CGAL::Scalar_factor_traits<Polynomial> SFT;
        typedef typename AT::Integer Scalar;
        CGAL_assertion_code(typedef typename SFT::Scalar Scalar_;)
        CGAL_static_assertion((::boost::is_same<Scalar_, Scalar>::value));
        
        typename SFT::Scalar_factor sfac;
        
        assert(sfac(Polynomial(0))==Scalar(0));
        assert(sfac(Polynomial(9))==Scalar(9));
        assert(sfac(Polynomial(9,15,30))==Scalar(3));
        
        assert(sfac(Polynomial(0), Integer(0))==Scalar(0));
        assert(sfac(Polynomial(0), Integer(1))==Scalar(1));
        assert(sfac(Polynomial(0), Integer(2))==Scalar(2));
        
        assert(sfac(Polynomial(9), Integer(0))==Scalar(9));
        assert(sfac(Polynomial(9), Integer(6))==Scalar(3));
        
        assert(sfac(Polynomial(15,0 ,30) , Integer(9))==Scalar(3));
        assert(sfac(Polynomial(0 ,15,30) , Integer(9))==Scalar(3));
        assert(sfac(Polynomial(18,15,30) , Integer(0))==Scalar(3));
    }{
        typedef typename AT::Integer Integer;
        typedef CGAL::Sqrt_extension<Integer,Integer> EXT_1;
        typedef CGAL::Polynomial<EXT_1       > Poly_1_ext_1;
        typedef CGAL::Polynomial<Poly_1_ext_1> Poly_2_ext_1;
        typedef CGAL::Scalar_factor_traits<Poly_2_ext_1> SFT;
        CGAL_assertion_code(typedef typename AT::Integer Scalar;)
        CGAL_assertion_code(typedef typename SFT::Scalar Scalar_;)
        CGAL_static_assertion((::boost::is_same<Scalar_, Scalar>::value));
            
        typename SFT::Scalar_factor sfac;

        assert(sfac(Poly_2_ext_1( ))==Integer(0));
        assert(sfac(Poly_2_ext_1(1))==Integer(1));
        assert(sfac(Poly_2_ext_1(2))==Integer(2));
        
        EXT_1 a(18,27,456); 
        EXT_1 b( 0,15,456);
        Poly_2_ext_1 p(Poly_1_ext_1(a,b),Poly_1_ext_1(b,a));
        
        assert(sfac(p)==Integer(3));
  }
}


template<class Polynomial_d, class T>
void test_interoperable_with(){
    // construction
  assert(Polynomial_d(-2) == T(-2));

  // operators 
  Polynomial_d f(4); 
  assert( f + T(2) == Polynomial_d( 6));
  assert( f - T(2) == Polynomial_d( 2));
  assert( f * T(2) == Polynomial_d( 8));
  assert( f / T(2) == Polynomial_d( 2));
  
  assert( T(2) + f == Polynomial_d( 6));
  assert( T(2) - f == Polynomial_d(-2));
  assert( T(2) * f == Polynomial_d( 8));
  assert( T(8) / f == Polynomial_d( 2));

  f+=T(2);
  f-=T(4);
  f*=T(3);
  f/=T(2);
  assert( f == Polynomial_d(3));
  

  assert(  T(2) == Polynomial_d(2));
  assert(  T(2) != Polynomial_d(4));
    
  assert( (T(2) <  Polynomial_d(4)));
  assert(!(T(4) <  Polynomial_d(4)));
  assert(!(T(6) <  Polynomial_d(4)));

  assert(!(T(2) >  Polynomial_d(4)));
  assert(!(T(4) >  Polynomial_d(4)));
  assert( (T(6) >  Polynomial_d(4)));

  assert(  T(2) <= Polynomial_d(4));
  assert(  T(4) <= Polynomial_d(4));
  assert(!(T(6) <= Polynomial_d(4)));

  assert(!(T(2) >= Polynomial_d(4)));
  assert( (T(4) >= Polynomial_d(4)));
  assert( (T(6) >= Polynomial_d(4)));
}

template<class Polynomial_d>
void test_interoperable_poly(){
  typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
  typedef typename PT::Coefficient_type                Coefficient_type;
  typedef typename PT::Innermost_coefficient_type Innermost_coefficient_type;

  test_interoperable_with<Polynomial_d,int>();
  test_interoperable_with<Polynomial_d,Innermost_coefficient_type>();
  test_interoperable_with<Polynomial_d,Coefficient_type>();
  test_interoperable_with<Polynomial_d,Polynomial_d>();
  
  
}

template<class AT>
void test_interoperable_at(){
  typedef typename AT::Integer Integer;
  typedef CGAL::Sqrt_extension<Integer,Integer> EXT; 
  {
    typedef int Coefficient_type;
    typedef CGAL::Polynomial<Coefficient_type>    Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    test_interoperable_poly<Poly_1>();
    test_interoperable_poly<Poly_2>();
  }{
    typedef Integer Coefficient_type;
    typedef CGAL::Polynomial<Coefficient_type>    Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    test_interoperable_poly<Poly_1>();
    test_interoperable_poly<Poly_2>();
  }{
    typedef EXT Coefficient_type;
    typedef CGAL::Polynomial<Coefficient_type>    Poly_1;
    typedef CGAL::Polynomial<Poly_1> Poly_2;
    test_interoperable_poly<Poly_1>();
    test_interoperable_poly<Poly_2>();
  }
}


template <class AT>
void test_AT(){
  
  test_interoperable_at<AT>();

    basic_tests<int>();
    basic_tests<double>();
    {
        typedef typename AT::Integer Integer;
        typedef typename CGAL::Polynomial<Integer> Polynomial;
        typedef CGAL::Unique_factorization_domain_tag Tag;
        typedef CGAL::Tag_true Is_exact;
        CGAL::test_algebraic_structure<Polynomial,Tag, Is_exact>();
        basic_tests<Integer>();
        exact_tests<Integer>();
    }{
        typedef typename AT::Rational Rational;
//        typedef typename CGAL::Polynomial<Rational> Polynomial;
//        typedef CGAL::Euclidean_ring_tag Tag;
//        typedef CGAL::Tag_true Is_exact;
        //can't use this test for Polynomials 
        //CGAL::test_algebraic_structure<Polynomial,Tag, Is_exact>();
        basic_tests<Rational>();
        exact_tests<Rational>();
        integr<Rational>(); 
    }

    // TODO: This also leads to trouble with flattening. Needs to be fixed
    psqff<AT>();
    test_evaluate<AT>();
    test_evaluate_homogeneous<AT>();
    
    test_total_degree();
    //test_scalar_factor_traits<AT>();
}
