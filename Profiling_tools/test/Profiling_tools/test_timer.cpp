#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>
#include <CGAL/FPU.h>

#include <cassert>

template <class T>
double test_timer() {
    T t;
    assert( ! t.is_running());
    t.start();
    assert( t.is_running());
    t.reset();
    assert( t.is_running());
    t.stop();
    assert( ! t.is_running());
    std::cout.precision(17);
    std::cout << "time()                 : " << t.time() << "\n";
    assert( t.time() >= 0.0);
    assert( t.intervals() == 1);
    assert( t.precision() >= 0.0);
    std::cout << "precision()            : " << t.precision() << "\n"; 
    assert( (t.max)() > 0.0);

    T s;
    s.start();
    double p = 0.0;     
    for ( int i = 0; i < 5; i++) {
        for ( int j = 0; j < 1000000; j++)
            p = p + 1.0;
        std::cout << "time() in " << (i+1) << ". iteration : " << s.time()
                  << "\n";
    }
    // make use of the computed result to stop optimizers from removing
    // the complete loop.
    return p;
}

int main(){
    // The following ensures that the conversion from internal timer
    // to double is done in a consistent way.  Otherwise we may get
    // double-rounding effects (show up on x86 at -O2 typically), which
    // makes the timer non-monotone.
    CGAL::Set_ieee_double_precision pfr;
    double p = 0.0;
    std::cout << "CGAL::Timer:\n";
    p += test_timer<CGAL::Timer>();
    std::cout << std::endl;
    std::cout << "CGAL::Real_timer:\n";
    p += test_timer<CGAL::Real_timer>();
    std::cout << std::endl;
    return (p > 0) ? 0 : 1;
}
// EOF //
