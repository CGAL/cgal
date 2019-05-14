#include <CGAL/basic.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/IO/io.h>
#include <vector>

template <class Fraction>
std::vector<typename CGAL::Fraction_traits<Fraction>::Numerator_type >
integralize(
        const std::vector<Fraction>& vec,
        typename CGAL::Fraction_traits<Fraction>::Denominator_type& d
) {
    typedef CGAL::Fraction_traits<Fraction> FT;
    typedef typename FT::Numerator_type Numerator_type;
    typedef typename FT::Denominator_type Denominator_type;
    typename FT::Decompose decompose;

    std::vector<Numerator_type>   num(vec.size());
    std::vector<Denominator_type> den(vec.size());

    // decompose each coefficient into integral part and denominator
    for (unsigned int i = 0; i < vec.size(); i++) {
        decompose(vec[i], num[i], den[i]);
    }

    // compute 'least' common multiple of all denominator
    // We would like to use gcd, so let's think of Common_factor as gcd.
    typename FT::Common_factor        gcd;
    d = 1;
    for (unsigned int i = 0; i < vec.size(); i++) {
        d *= CGAL::integral_division(den[i], gcd(d, den[i]));
    }

    // expand each (numerator, denominator) pair to common denominator
    for (unsigned int i = 0; i < vec.size(); i++) {
        // For simplicity ImplicitInteroperability is expected in this example
        num[i] *= CGAL::integral_division(d, den[i]);
    }
    return num;
}

#ifdef CGAL_USE_GMP

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

int main(){
    std::vector<CGAL::Gmpq> vec(3);
    vec[0]=CGAL::Gmpq(1,4);
    vec[1]=CGAL::Gmpq(1,6);
    vec[2]=CGAL::Gmpq(1,10);
    std::cout<< "compute an integralized vector" << std::endl;
    std::cout<<"input vector:  ["
             << vec[0] << "," << vec[1] << "," << vec[2] << "]" << std::endl;
    CGAL::Gmpz d;
    std::vector<CGAL::Gmpz> integral_vec = integralize(vec,d);
    std::cout<<"output vector: ["
             << integral_vec[0] << ","
             << integral_vec[1] << ","
             << integral_vec[2] << "]" << std::endl;
    std::cout<<"denominator  : "<< d <<std::endl;
}
#else
int main(){ std::cout << "This examples needs GMP" << std::endl; }
#endif
