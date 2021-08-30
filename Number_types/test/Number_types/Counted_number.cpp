#include <iostream>
#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/Counted_number.h>
#include <CGAL/Algebraic_structure_traits.h>

#include<CGAL/int.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>

template<class T, class Tag >
void test_counted_number(T,Tag){
    typedef CGAL::Counted_number<T> NT;
    typedef typename CGAL::Algebraic_structure_traits<T>::Is_exact Is_exact;

     CGAL::test_algebraic_structure<NT,Tag, Is_exact>();

     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6),NT(15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(6),NT(-15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(6), NT(15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(4),NT(-6),NT(-15));
     CGAL::test_algebraic_structure<NT,Tag, Is_exact>(NT(-4),NT(-6),NT(-15));


     CGAL::test_real_embeddable<NT>();
}

template< class AK >
void call_tests_with_types_from_ak() {
  test_counted_number( typename AK::Integer(), CGAL::Euclidean_ring_tag() );
  test_counted_number( typename AK::Rational(), CGAL::Field_tag() );
  test_counted_number( typename AK::Field_with_sqrt(),
                       typename CGAL::Algebraic_structure_traits<
                                                    typename AK::Field_with_sqrt
                                                 >::Algebraic_category() );
}

int main() {
    test_counted_number(int(), CGAL::Euclidean_ring_tag());
    test_counted_number(double(),CGAL::Field_with_kth_root_tag()); // works
#ifdef CGAL_USE_LEDA
    call_tests_with_types_from_ak< CGAL::LEDA_arithmetic_kernel >();
#endif // CGAL_USE_LEDA
#ifdef CGAL_USE_CORE
    call_tests_with_types_from_ak< CGAL::CORE_arithmetic_kernel >();
#endif //CGAL_USE_CORE
    return 0;
}
