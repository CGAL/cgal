
#include "test_polynomial.h"
#include <CGAL/CORE_arithmetic_kernel.h>

int main() {

    // Set wrong rounding mode to test modular arithmetic
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

    CGAL::IO::set_pretty_mode(std::cout);

#ifdef CGAL_USE_CORE
    {
        typedef CGAL::CORE_arithmetic_kernel AT;
        test_AT<AT>();
    }
#endif

    return 0;
}


