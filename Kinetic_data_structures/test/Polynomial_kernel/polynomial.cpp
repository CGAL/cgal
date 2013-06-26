#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/basic.h>
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Fixed_polynomial.h>
#include <CGAL/Polynomial/internal/Rational/Rational_traits_base.h>

//#include <CGAL/Polynomial/internal/Filtered_rational/Filtered_rational_traits.h>
//#include <CGAL/Polynomial/Default_filtering_traits.h>
#ifdef CGAL_USE_CORE
#include <CGAL/Polynomial/CORE_kernel.h>
#include <CGAL/CORE_Expr.h>
#endif

#include <vector>

//#include "write_maple_functions.h"

bool for_maple=false;

template <class Polynomial>
void write(const char *expr, const Polynomial &v)
{
  if (for_maple) {
    std::cout << "evalb(simplify(expand(" << expr << ")) = " << v << ");\n";
  }
  else {
    std::cout << expr << ": " << v << std::endl;
  }
}


template <class Polynomial>
void write_variable(const char *name, const Polynomial &v)
{
  if (for_maple) {
    std::cout << name << ":=" << v << ":\n";
  }
  else {
    std::cout << name << ": " << v << std::endl;
  }
}

template <class Traits>
void check_equal(const typename Traits::Function&a, const typename Traits::Function &b) {
  if (a != b) {
    std::cerr << a <<  " != " << b << std::endl;
    assert(a==b);
  }
}


template <class Traits>
void test_polynomial(const Traits &tr)
{
  typedef typename Traits::Construct_function CF;
  typedef typename Traits::Function Polynomial;
  typedef typename Polynomial::NT NT;

  CF cf= tr.construct_function_object();

  NT v[] = {-1, 2, 27, -17, 0, 0};
  //  NT v[] = {0, 0, 0, 0, 0, 0};

  Polynomial p(v, v+6);
  Polynomial q(v, v+3);

  NT a(2);

  write_variable( "p", p);
  write_variable("q", q );

  write("-p", (-p));
  check_equal<Traits>(-p , cf(1,-2,-27,17));

  write("p-p",(p-p));
  check_equal<Traits>(p-p , cf(0));

  write("p+q" , (p+q) );
  check_equal<Traits>(p+q , cf( -2, 4, 54, -17));

  write("p-q" , (p-q));
  check_equal<Traits>(p-q , cf(0,0,0,-17));

  write("q*(p-q)" , q*(p-q) );
  check_equal<Traits>(q*(p-q) , cf(0,0,0,17,-34,-459));

  write_variable( "a", a);

  Polynomial diff(p-q);
  Polynomial dp(diff +Polynomial(a));

  write("(p-q)+a" , dp );
  check_equal<Traits>((p-q)+Polynomial(a) , cf(2, 0, 0, -17));

  write("(p-q)-a" , ((p-q)-a) );
  check_equal<Traits>((p-q)-Polynomial(a) , cf(-2, 0, 0, -17));

  write("a*(p-q)" , (a*(p-q)) );
  check_equal<Traits>((Polynomial(a)*(p-q)) , cf(0,0,0,-34));

  write("(p-q)*a" , ((p-q)*a) );
  check_equal<Traits>(((p-q)*Polynomial(a)) , cf(0,0,0,-34));

  write("(p-q)/a" , ((p-q)/a) );
  check_equal<Traits>((Polynomial(p-q)/a) , cf(0,0,0,-NT(.5)*NT(17)));

  write("subs(t=-t, p)", tr.negate_variable_object()(p) );
  check_equal<Traits>(tr.negate_variable_object()(p) , cf( -1, -2, 27, 17));

  /*write("t^degree(p) * subs(t=(1/t), p)", tr.invert_variable_object()(p) );
    check_equal<Traits>( tr.invert_variable_object()(p) , cf(-17, 27, 2, -1));*/


  NT v1[] = {-1, 1};
  NT v2[] = {-2, 1};

  Polynomial r = Polynomial(v1, v1+2);
  Polynomial s = Polynomial(v2, v2+2);

  p = r * r * s + Polynomial(NT(1));
  write_variable( "p", p);
  check_equal<Traits>(p , cf(-1, 5, -4, 1));

  q = r * s;
  write_variable("q", q );
  check_equal<Traits>(q , cf( 2, -3, 1));

  write("rem(p,q,t)",
        tr.remainder_object()(p,q) );
  check_equal<Traits>(tr.remainder_object()(p,q) , cf(1));

  write("prem(p,q,t)",
        tr.pseudo_remainder_object()(p,q) );
  check_equal<Traits>(tr.pseudo_remainder_object()(p,q) , cf(1));

  write("quo(p,q,t)",
        tr.quotient_object()(p,q) );
  check_equal<Traits>(tr.quotient_object()(p,q) , cf(-1,1));

  write("pquo(p,q,t)",
        tr.pseudo_quotient_object()(p,q) );
  check_equal<Traits>(tr.pseudo_quotient_object()(p,q) , cf(-1,1));

  p = r * r * s * s * s;
  write_variable( "p", p);
  check_equal<Traits>(p , cf(-8, 28, -38, 25, -8, 1));

  q = r * s;
  write_variable("q", q );
  check_equal<Traits>(q , cf(2,-3,1));

  write("rem(p,q,t)",
        tr.remainder_object()(p,q) );
  check_equal<Traits>(tr.remainder_object()(p,q) , cf(0));

  write("prem(p,q,t)",
        tr.pseudo_remainder_object()(p,q) );
  check_equal<Traits>(tr.pseudo_remainder_object()(p,q) , cf(0));

  write("quo(p,q,t)",
        tr.quotient_object()(p,q) );
  check_equal<Traits>(tr.quotient_object()(p,q) , cf(-4, 8, -5, 1));

  write("pquo(p,q,t)",
        tr.pseudo_quotient_object()(p,q) );
  check_equal<Traits>(tr.pseudo_quotient_object()(p,q) , cf(-4, 8, -5, 1));

  /*  int shift = 6;

  write("p * t^6", tr.shift_power_object(shift)(p) );
  check_equal<Traits>(tr.shift_power_object(shift)(p)
  , cf(0,0,0,0,0,0,-8, 28, -38, 25, -8, 1));*/

  NT v3[] = {0, 1};

  Polynomial t = Polynomial(v3, v3+2);
  p = t * t;

  write_variable( "p", p);

  /*NT new_zero = NT(-1);
  write("subs(t=t-1,p)", tr.rational_translate_zero_object(new_zero)(p));
  check_equal<Traits>(tr.rational_translate_zero_object(new_zero)(p)
  , cf(1, -2, 1));*/
}


int main(int argc, char* argv[])
{

  //CORE::extLong pi=CORE_posInfty;
  //  CORE::Expr ep(CORE_posInfty);

  //std::cout << /*pi << " " <<*/ ep << std::endl;

  if ( argc > 1 ) {
    for_maple = (atoi(argv[1]) != 0);
  }

  /*
    typedef CORE::BigRat NT;
    typedef CORE::Polynomial<NT> P;
    NT pc[4];
    pc[0]=NT("-1/1");
    pc[1]=NT("2/1");
    pc[2]=NT("27/1");
    pc[3]=NT("-17/1");
    P p(3, pc);
    NT qc[4];
    qc[0]=NT("-1/1");
    qc[1]=NT("2/1");
    qc[2]=NT("27/1");
    P q(2, qc);
    std::cout << p << " " << q << std::endl;
    NT a("2");
    P r= (p+q) + a;
    }*/

  /*if (for_maple){
    write_maple_functions(std::cout);
    }*/
  {
    std::cout << "Testing regular poly.\n";
    typedef CGAL::POLYNOMIAL::Default_field_nt                        NT;
    typedef CGAL_POLYNOMIAL_NS::Polynomial<NT>     Polynomial;

    typedef CGAL_POLYNOMIAL_NS::internal::Rational_traits_base<Polynomial>
      Rational_traits;
    Rational_traits tr;
    test_polynomial(tr);
  }

  std::cout <<"\n\n\n\n\n";

  /*{
    std::cout << "Testing filtered poly.\n";
    typedef CGAL::POLYNOMIAL::Default_field_nt NT;
    typedef CGAL_POLYNOMIAL_NS::Default_filtering_traits<NT> FT;
    typedef CGAL_POLYNOMIAL_NS::internal::Filtered_rational_traits<FT> Tr;
    Tr tr;
    test_polynomial(tr);
    }*/
  std::cout <<"\n\n\n\n\n";
#ifdef CGAL_USE_CORE
  {
    std::cout << "Testing core poly.\n";
    typedef CGAL_POLYNOMIAL_NS::CORE_kernel CORE_kernel;
    CORE_kernel tr;
    test_polynomial(tr);
  }
  std::cout <<"\n\n\n\n\n";
#endif
  {
    std::cout << "Testing fixed poly.\n";
    typedef CGAL::POLYNOMIAL::Default_field_nt NT;
    typedef CGAL_POLYNOMIAL_NS::Fixed_polynomial<NT, 30> Poly;
    typedef CGAL_POLYNOMIAL_NS::internal::Rational_traits_base<Poly>
      Rational_traits;
    Rational_traits tr;
    test_polynomial(tr);
  }
  return 0;
}
