

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

void test_timer() {
    CGAL::Timer t;
    CGAL_assertion( ! t.is_running());
    t.start();
    CGAL_assertion( t.is_running());
    t.reset();
    CGAL_assertion( t.is_running());
    t.stop();
    CGAL_assertion( ! t.is_running());
    CGAL_assertion( t.time() >= 0.0);
    std::cout << t.time() << "\n";
    CGAL_assertion( t.intervals() == 1);
    CGAL_assertion( t.precision() >= 0.0);
    std::cout << t.precision() << "\n"; 
    CGAL_assertion( t.max() > 0.0);
}

void test_real_timer() {
    CGAL::Real_timer t;
    CGAL_assertion( ! t.is_running());
    t.start();
    CGAL_assertion( t.is_running());
    t.reset();
    CGAL_assertion( t.is_running());
    t.stop();
    CGAL_assertion( ! t.is_running());
    CGAL_assertion( t.time() >= 0.0);
    std::cout << t.time() << "\n";
    CGAL_assertion( t.intervals() == 1);
    CGAL_assertion( t.max() > 0.0);
}

int main(){
    test_timer();
    test_real_timer();
    CGAL::Real_timer Tm;
    CGAL::Timer Tm2;
    Tm.start();
    Tm2.start();
    int i1,i2;
    double p = 0;
     
    std::cout << "\n";
    for (i1=0;i1<5;i1++) {
    for (i2=0;i2<800000;i2++) p=p+1.0 ;
     std::cout << Tm.time() << "\n";
    }
    std::cout << "\n";

    for (i1=0;i1<5;i1++) {
    for (i2=0;i2<800000;i2++) p=p+1.0 ;
     std::cout << Tm2.time() << "\n";
    }
    std::cout << "\n";
  
    return 0;
}
// EOF //
