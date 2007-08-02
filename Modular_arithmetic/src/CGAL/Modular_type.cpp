//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

#include <CGAL/Modular_arithmetic/Modular_type.h>

namespace CGAL{
    int Modular::prime_int = 67111067;
    double Modular::prime =67111067.0;
    double Modular::prime_inv =1/67111067.0;
    
    const double Modular::CST_CUT = std::ldexp( 3., 51 );

}
