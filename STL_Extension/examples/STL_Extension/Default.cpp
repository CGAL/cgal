#include <CGAL/Default.h>

struct A {};

template < typename A1_ = A, // we could also write CGAL::Default here instead
           typename A2 = int >
struct B
{
    typedef typename CGAL::Default::Get<A1_, A>::type  A1;
    A1 a1;
};

int main ()
{
    B<CGAL::Default, double> b;

    A a = b.a1; //  It is really of type A.
}
