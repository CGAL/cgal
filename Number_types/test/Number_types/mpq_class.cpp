#include <CGAL/config.h>

#ifdef CGAL_USE_GMPXX

#include <iostream>
#include <CGAL/mpq_class.h>
#include <CGAL/Lazy_exact_nt.h>

#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_fraction_traits.h>
#include <CGAL/Test/_test_rational_traits.h>

int main() {
    {   typedef mpq_class NT;
        typedef CGAL::Field_tag Tag;
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
        CGAL::test_fraction_traits<NT>();
        // backward compatiblity
        CGAL::test_rational_traits<NT>();
    }
    {
      mpq_class q;
      std::istringstream in("12.34");
      in >> CGAL::IO::iformat(q);
      assert(in);
      assert(q.get_num() == 617);
      assert(q.get_den() == 50);
    }
    {
      CGAL::Lazy_exact_nt<mpq_class> x;
      std::istringstream in("12.34");
      in >> x;
      mpq_class q = x.exact();
      assert(in);
      assert(q.get_num() == 617);
      assert(q.get_den() == 50);
    }
    return 0;
}


#else
int main()
{
  return 0;
}
#endif //CGAL_USE_GMPXX
