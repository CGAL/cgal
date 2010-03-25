#include <CGAL/Default.h>

// A is a concrete type
struct A {};

// B is the template class which has 2 template parameters
// with default arguments : A and int.
template < typename A1_ = A, typename A2 = int >
struct B
{
    B()
      : a1()
    {}

    // Note that it is also possible to use CGAL::Default
    // instead of A as the default argument for A1_ above.

    // Extract the desired type for A1 :
    typedef typename CGAL::Default::Get<A1_, A>::type  A1;

    A1 a1;
};

int main ()
{
    B<CGAL::Default, double> b;

    A a = b.a1; //  It is really of type A.
}
