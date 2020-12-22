#include <CGAL/config.h>
#include <iostream>

#if (!CGAL_USE_CORE) && (!CGAL_USE_LEDA) && (!(CGAL_USE_GMP && CGAL_USE_MPFI))
int main ()
{
  std::cout << "Sorry, this example needs CORE, LEDA, or GMP+MPFI ..."
            << std::endl;
  return 0;
}
#else

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#if CGAL_USE_GMP && CGAL_USE_MPFI
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz Integer;
#elif CGAL_USE_CORE
#include <CGAL/CORE_BigInt.h>
typedef CORE::BigInt Integer;
#else
#include <CGAL/leda_integer.h>
typedef LEDA::integer Integer;
#endif

typedef CGAL::Arr_algebraic_segment_traits_2<Integer> Arr_traits_2;
typedef CGAL::Arrangement_2<Arr_traits_2> Arrangement_2;
typedef Arr_traits_2::Curve_2 Curve_2;
typedef Arr_traits_2::Polynomial_2 Polynomial_2;

int main() {

    // For nice printouts
    CGAL::set_pretty_mode(std::cout);

    Arr_traits_2 arr_traits;

    // Functor to create a curve from a Polynomial_2
    Arr_traits_2::Construct_curve_2 construct_curve
        = arr_traits.construct_curve_2_object();

    Polynomial_2 x = CGAL::shift(Polynomial_2(1),1,0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1),1,1);

    Arrangement_2 arr(&arr_traits);

    // Construct an (unbounded line) with equation 3x-5y+2=0
    Polynomial_2 f1 = 3*x-5*y+2;
    Curve_2 cv1 = construct_curve(f1);
    std::cout << "Adding curve " << f1 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv1);

    // Construct the ellipse x^2+3*y^2-10=0
    Polynomial_2 f2 = CGAL::ipower(x,2)+3*CGAL::ipower(y,2)-10;
    Curve_2 cv2 = construct_curve(f2);
    std::cout << "Adding curve " << f2 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv2);

    // Construct a cubic curve with isoated point, and vertical asymptote
    // x^2+y^2+xy^2
    Polynomial_2 f3 = CGAL::ipower(x,2)+CGAL::ipower(y,2)+x*CGAL::ipower(y,2);
    Curve_2 cv3 = construct_curve(f3);
    std::cout << "Adding curve " << f3 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv3);

    // Construct a curve of degree 6 with equation x^6+y^6-x^3y^3-12
    Polynomial_2 f4 = CGAL::ipower(x,6)+CGAL::ipower(y,6)-
                      CGAL::ipower(x,3)*CGAL::ipower(y,3)-12;
    Curve_2 cv4 = construct_curve(f4);
    std::cout << "Adding curve " << f4 << " to the arrangement" << std::endl;
    CGAL::insert(arr,cv4);

    // Print the arrangement size.
    std::cout << "The arrangement size:" << std::endl
              << "   V = " << arr.number_of_vertices()
              << ",  E = " << arr.number_of_edges()
              << ",  F = " << arr.number_of_faces() << std::endl;

    return 0;
}

#endif
