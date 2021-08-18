#include <cassert>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Test/_test_coercion_traits.h>

template<class A, class B> void test_coercion_from_to(A, B){
  CGAL::test_explicit_interoperable<A,B,B>();
}

template <typename AK>
void test_coercion_traits(){
  //CGAL::IO::set_pretty_mode(std::cout);

  typedef typename AK::Integer  Integer;
  typedef typename AK::Rational Rational;

  typedef typename CGAL::Polynomial_type_generator<Integer, 1>::Type POLY_INT_1;
  typedef typename CGAL::Polynomial_type_generator<Integer, 2>::Type POLY_INT_2;
  typedef typename CGAL::Polynomial_type_generator<Integer, 3>::Type POLY_INT_3;
  typedef typename CGAL::Polynomial_type_generator<Rational, 1>::Type POLY_RAT_1;
  typedef typename CGAL::Polynomial_type_generator<Rational, 2>::Type POLY_RAT_2;
  typedef typename CGAL::Polynomial_type_generator<Rational, 3>::Type POLY_RAT_3;

  test_coercion_from_to(POLY_INT_1(), POLY_INT_1());
  test_coercion_from_to(POLY_INT_1(), POLY_INT_2());
  test_coercion_from_to(POLY_INT_1(), POLY_INT_3());
  test_coercion_from_to(POLY_INT_2(), POLY_INT_2());
  test_coercion_from_to(POLY_INT_2(), POLY_INT_3());
  test_coercion_from_to(POLY_INT_3(), POLY_INT_3());

  test_coercion_from_to(POLY_RAT_1(), POLY_RAT_1());
  test_coercion_from_to(POLY_RAT_1(), POLY_RAT_2());
  test_coercion_from_to(POLY_RAT_1(), POLY_RAT_3());
  test_coercion_from_to(POLY_RAT_2(), POLY_RAT_2());
  test_coercion_from_to(POLY_RAT_2(), POLY_RAT_3());
  test_coercion_from_to(POLY_RAT_3(), POLY_RAT_3());

  // Though the coercion type is clear, the problem is how to match the
  // variables. The recursive definition of Polynomial<Coeff> suggest that
  // the coercion type of two polynomial types Polynomial<A> and Polynomial<B>
  // is defined as Polynomial<C>, where C is the coercion type.
  // However, this is not in line with the fact that a Polynomial<A>
  // is interoperable with its coefficient type A, that is, if A is a polynomial
  // the variables of A should not be moved outward while casting A to
  // Polynomial<A>. This is tested in the sequel.
  {
    POLY_INT_1 x_1 = CGAL::shift(POLY_INT_1(1),1,0);

    POLY_INT_3 x_3 = CGAL::shift(POLY_INT_3(1),1,0);
    POLY_INT_3 y_3 = CGAL::shift(POLY_INT_3(1),1,1);
    POLY_INT_3 z_3 = CGAL::shift(POLY_INT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_INT_1,POLY_INT_3> CT;
    assert(typename CT::Cast()(x_1) == x_3);
    assert(typename CT::Cast()(x_1) != y_3);
    assert(typename CT::Cast()(x_1) != z_3);
    assert(typename CT::Cast()(x_3) == x_3);
    assert(typename CT::Cast()(y_3) == y_3);
    assert(typename CT::Cast()(z_3) == z_3);
  }{
    POLY_INT_1 x_1 = CGAL::shift(POLY_INT_1(1),1,0);

    POLY_RAT_3 x_3 = CGAL::shift(POLY_RAT_3(1),1,0);
    POLY_RAT_3 y_3 = CGAL::shift(POLY_RAT_3(1),1,1);
    POLY_RAT_3 z_3 = CGAL::shift(POLY_RAT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_INT_1,POLY_RAT_3> CT;
    assert(typename CT::Cast()(x_1) == x_3);
    assert(typename CT::Cast()(x_1) != y_3);
    assert(typename CT::Cast()(x_1) != z_3);
    assert(typename CT::Cast()(x_3) == x_3);
    assert(typename CT::Cast()(y_3) == y_3);
    assert(typename CT::Cast()(z_3) == z_3);
  }{
    POLY_RAT_1 x_1 = CGAL::shift(POLY_RAT_1(1),1,0);

    POLY_INT_3 x_3 = CGAL::shift(POLY_INT_3(1),1,0);
    POLY_INT_3 y_3 = CGAL::shift(POLY_INT_3(1),1,1);
    POLY_INT_3 z_3 = CGAL::shift(POLY_INT_3(1),1,2);

    POLY_RAT_3 x_3r = CGAL::shift(POLY_RAT_3(1),1,0);
    POLY_RAT_3 y_3r = CGAL::shift(POLY_RAT_3(1),1,1);
    POLY_RAT_3 z_3r = CGAL::shift(POLY_RAT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_RAT_1,POLY_INT_3> CT;
    assert(typename CT::Cast()(x_1) == x_3r);
    assert(typename CT::Cast()(x_1) != y_3r);
    assert(typename CT::Cast()(x_1) != z_3r);
    assert(typename CT::Cast()(x_3) == x_3r);
    assert(typename CT::Cast()(y_3) == y_3r);
    assert(typename CT::Cast()(z_3) == z_3r);
  }

  {
    POLY_INT_1 x_1 = CGAL::shift(POLY_INT_1(1),1,0);

    POLY_INT_3 x_3 = CGAL::shift(POLY_INT_3(1),1,0);
    POLY_INT_3 y_3 = CGAL::shift(POLY_INT_3(1),1,1);
    POLY_INT_3 z_3 = CGAL::shift(POLY_INT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_INT_3,POLY_INT_1> CT;
    assert(typename CT::Cast()(x_1) == x_3);
    assert(typename CT::Cast()(x_1) != y_3);
    assert(typename CT::Cast()(x_1) != z_3);
    assert(typename CT::Cast()(x_3) == x_3);
    assert(typename CT::Cast()(y_3) == y_3);
    assert(typename CT::Cast()(z_3) == z_3);
  }{
    POLY_INT_1 x_1 = CGAL::shift(POLY_INT_1(1),1,0);

    POLY_RAT_3 x_3 = CGAL::shift(POLY_RAT_3(1),1,0);
    POLY_RAT_3 y_3 = CGAL::shift(POLY_RAT_3(1),1,1);
    POLY_RAT_3 z_3 = CGAL::shift(POLY_RAT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_RAT_3,POLY_INT_1> CT;
    assert(typename CT::Cast()(x_1) == x_3);
    assert(typename CT::Cast()(x_1) != y_3);
    assert(typename CT::Cast()(x_1) != z_3);
    assert(typename CT::Cast()(x_3) == x_3);
    assert(typename CT::Cast()(y_3) == y_3);
    assert(typename CT::Cast()(z_3) == z_3);
  }{
    POLY_RAT_1 x_1 = CGAL::shift(POLY_RAT_1(1),1,0);

    POLY_INT_3 x_3 = CGAL::shift(POLY_INT_3(1),1,0);
    POLY_INT_3 y_3 = CGAL::shift(POLY_INT_3(1),1,1);
    POLY_INT_3 z_3 = CGAL::shift(POLY_INT_3(1),1,2);

    POLY_RAT_3 x_3r = CGAL::shift(POLY_RAT_3(1),1,0);
    POLY_RAT_3 y_3r = CGAL::shift(POLY_RAT_3(1),1,1);
    POLY_RAT_3 z_3r = CGAL::shift(POLY_RAT_3(1),1,2);

    typedef CGAL::Coercion_traits<POLY_INT_3,POLY_RAT_1> CT;
    assert(typename CT::Cast()(x_1) == x_3r);
    assert(typename CT::Cast()(x_1) != y_3r);
    assert(typename CT::Cast()(x_1) != z_3r);
    assert(typename CT::Cast()(x_3) == x_3r);
    assert(typename CT::Cast()(y_3) == y_3r);
    assert(typename CT::Cast()(z_3) == z_3r);
  }
  /*
    {
    typedef CGAL::Coercion_traits<POLY_RAT_1,CGAL::Null_functor> CT;
    CGAL_static_assertion((
    ::boost::is_same< typename CT::Are_implicit_interoperable,
    CGAL::Tag_false>::value));
    }
  */

}

int main(){
#if CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
  typedef CGAL::Arithmetic_kernel AK;
  test_coercion_traits<AK>();
#endif
return 0;
}
