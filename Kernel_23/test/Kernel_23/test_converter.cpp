#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Interval_nt.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Homogeneous_converter.h>

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_Expr.h>
#endif

#include <iostream>

int main()
{
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
  using SCLD = CGAL::Simple_cartesian<long double>;
  using SCI = CGAL::Simple_cartesian<CGAL::Interval_nt<> >;
  using SHI = CGAL::Simple_homogeneous<CGAL::Interval_nt<> >;
  using NT_exact = CGAL::internal::Exact_field_selector<double>::Type;
  using SHD = CGAL::Simple_homogeneous<double>;
  using SHE = CGAL::Simple_homogeneous<NT_exact>;

  CGAL::Cartesian_converter<SCI, EPICK> sci_to_epick;
  CGAL::Cartesian_converter<SCLD, EPICK> scld_to_epick;
//  CGAL::Cartesian_converter<SCI, EPECK> sci_to_epeck;
  CGAL::Cartesian_converter<EPECK, EPICK> epeck_to_epick;
  CGAL::Homogeneous_converter<SHI, SHD> shi_to_shd;
  CGAL::Homogeneous_converter<SHE, SHD> she_to_shd;
  CGAL::Homogeneous_converter<SHE, SHE> she_to_she;

  assert(sci_to_epick(SCI::FT(2)) == EPICK::FT(2));
//  assert(sci_to_epick((long int)(2)) == EPICK::FT(2)); // static assertion expected
//  assert(sci_to_epick(2.) == EPICK::FT(2)); // static assertion expected
  assert(sci_to_epick(bool(2.)) == true);

  // fundamental type, but it's K1::FT so OK
  assert(scld_to_epick((long double)(2.)) == EPICK::FT(2));

#ifdef CGAL_USE_CORE
  using SSCE = CGAL::Simple_cartesian<CORE::Expr>;
  CGAL::Cartesian_converter<SSCE, EPICK> scce_to_epick;

  // fundamental types, so static assertion failures
//  assert(scce_to_epick(2.) == EPICK::FT(2));
//  assert(scce_to_epick((signed long long int)(1)) == true);
  scce_to_epick(CORE::Expr(2.) == EPICK::FT(2));
#endif

  // int* is not a fundamental type, so this goes through Enum_converter(bool)
  assert(sci_to_epick(true) == true);
  assert(sci_to_epick(false) == false);
  assert(sci_to_epick(CGAL::ON_POSITIVE_SIDE) == CGAL::ON_POSITIVE_SIDE);

  int a = 123;
  int* a_ptr = &a;
  CGAL_USE(a_ptr);
//  assert(sci_to_epick(a_ptr) == a_ptr); // unauthorized since 'int*' isn't a fundamental type

  // fundamental type, static assertion failure
//  assert(sci_to_epeck(2.) == EPECK::FT(2));

  assert(epeck_to_epick(EPECK::FT(2)) == EPICK::FT(2));
//  assert(epeck_to_epick(2.) == EPICK::FT(2));

  // Homogeneous
  assert(shi_to_shd(2.) == 2.);
  assert(shi_to_shd(SHI::RT(2.)) == SHD::RT(2.));

  assert(she_to_shd(SHE::RT(2)) == SHD::RT(2));
  assert(she_to_shd(SHE::FT(2.)) == SHD::FT(2.));
  assert(she_to_shd(2.) == 2.);
  assert(she_to_she(SHE::RT(2)) == SHE::RT(2));

  std::cout << "NT_exact is " << typeid(NT_exact).name() << std::endl;
  CGAL::NT_converter<CGAL::Quotient<NT_exact>, CGAL::Quotient<NT_exact> > qnte_to_qnte;
  assert(qnte_to_qnte(2.) == CGAL::Quotient<NT_exact>(2.));

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
