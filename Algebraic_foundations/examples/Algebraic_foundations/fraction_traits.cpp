#include <CGAL/Fraction_traits.h>
#include <CGAL/IO/io.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
int main(){
    typedef CGAL::Fraction_traits<CGAL::Gmpq> FT;
    typedef FT::Numerator_type Numerator_type;
    typedef FT::Denominator_type Denominator_type;

    CGAL_static_assertion((boost::is_same<Numerator_type,CGAL::Gmpz>::value));
    CGAL_static_assertion((boost::is_same<Denominator_type,CGAL::Gmpz>::value));

    Numerator_type numerator;
    Denominator_type denominator;
    CGAL::Gmpq fraction(4,5);
    FT::Decompose()(fraction,numerator,denominator);

    CGAL::set_pretty_mode(std::cout);
    std::cout << "decompose fraction: "<< std::endl;
    std::cout << "fraction   : " << fraction << std::endl;
    std::cout << "numerator  : " << numerator<< std::endl;
    std::cout << "denominator: " << denominator << std::endl;

    std::cout << "re-compose fraction: "<< std::endl;
    fraction = FT::Compose()(numerator,denominator);
    std::cout << "fraction   : " << fraction << std::endl;
}
#else
int main(){ std::cout << "This examples needs GMP" << std::endl; }
#endif
