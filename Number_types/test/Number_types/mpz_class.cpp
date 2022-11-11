#include <CGAL/config.h>
#include <CGAL/tags.h>
#ifdef CGAL_USE_GMPXX

#include <iostream>
#include <CGAL/mpz_class.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>



int main() {
  typedef mpz_class NT;
  typedef CGAL::Euclidean_ring_tag Tag;
  typedef CGAL::Tag_true Is_exact;
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

  return 0;
}
#else
int main()
{
  return 0;
}
#endif //CGAL_USE_GMPXX
