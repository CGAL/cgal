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
  using SCI = CGAL::Simple_cartesian<CGAL::Interval_nt<> >;
  using SHI = CGAL::Simple_homogeneous<CGAL::Interval_nt<> >;
  using NT_exact = CGAL::internal::Exact_field_selector<double>::Type;
  using SHD = CGAL::Simple_homogeneous<double>;
  using SHE = CGAL::Simple_homogeneous<NT_exact>;

  CGAL::Cartesian_converter<SCI, EPICK> sci_to_epick;
//  CGAL::Cartesian_converter<SCI, EPECK> sci_to_epeck;
  CGAL::Cartesian_converter<EPECK, EPICK> epeck_to_epick;
  CGAL::Homogeneous_converter<SHI, SHD> shi_to_shd;
  CGAL::Homogeneous_converter<SHE, SHD> she_to_shd;
  CGAL::Homogeneous_converter<SHE, SHE> she_to_she;

  assert(sci_to_epick(SCI::FT(2)) == EPICK::FT(2));
  assert(sci_to_epick((long int)(2)) == EPICK::FT(2)); // long int --> SCI::FT --> EPICK::FT
  assert(sci_to_epick(2.) == EPICK::FT(2)); // double --> SCI::FT --> EPICK::FT
  assert(sci_to_epick(bool(2.)) == true);

#ifdef CGAL_USE_CORE
  using SSCE = CGAL::Simple_cartesian<CORE::Expr>;
  CGAL::Cartesian_converter<SSCE, EPICK> scce_to_epick;

  // double is a fundamental type, but CORE::Expr has no implicit double --> CORE::Expr conversion
  assert(scce_to_epick(2) == EPICK::FT(2));

  // This does not compile because 'signed long long int' is a fundamental type,
  // but there is no CORE::Expr(signed long long int) constructor
//  assert(scce_to_epick((signed long long int)(1)) == true);
#endif

  // int* is not a fundamental type, so this goes through Enum_converter(bool)
  int a = 123;
  auto a_ptr = &a;
  assert(sci_to_epick(a_ptr) == true);

  // this does not compile because there is no conversion between Interval_nt and EPECK::FT
//  assert(sci_to_epeck(2.) == EPECK::FT(2));

  assert(epeck_to_epick(EPECK::FT(2)) == EPICK::FT(2));
  assert(epeck_to_epick(2.) == EPICK::FT(2)); // double --> EPECK::FT --> EPICK::FT

  assert(shi_to_shd(2.) == SHD::FT(2)); // double --> SHI::FT --> SHD::FT

  assert(she_to_shd(SHE::RT(2)) == SHD::RT(2));
  assert(she_to_shd(SHE::FT(2.)) == SHD::FT(2.));
  assert(she_to_shd(2.) == SHD::FT(2)); // double --> SHE::RT --> SHD::FT
  assert(she_to_she(2.) == 2.);

  std::cout << "NT_exact is " << typeid(NT_exact).name() << std::endl;
  CGAL::NT_converter<CGAL::Quotient<NT_exact>, CGAL::Quotient<NT_exact> > qnte_to_qnte;
  assert(qnte_to_qnte(2.) == CGAL::Quotient<NT_exact>(2.));

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
