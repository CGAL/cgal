#ifndef CGAL_USE_CLN
#  error You need to have CLN installed in order to run this example
#endif

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/CLN/cl_integer.h>

#include <cl_io.h>
#include <cl_integer_io.h>

int
main ()
{
    using std::cout;
    using std::endl;

    cl_I b = "31";
    cl_I p = b * 75;
    bool ok = (p==2325) && (CGAL::compare(p,cl_I(31*75))==CGAL::EQUAL);

    cout << (ok ? "ok" : "ERROR") << endl;
    cout << p << " == " << 31*75 << endl;
    return (ok) ? 0 : -1;
}
