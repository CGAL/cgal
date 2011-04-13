// This program tests some basic IO functions for some CLN types.

#include <CGAL/basic.h>
#include <iostream>

#ifndef CGAL_USE_CLN
int main()
{
  std::cout << "CGAL was not installed with CLN support." << std::endl;
  return 0;
}
#else // CGAL_USE_CLN

#include <fstream>
#include <cassert>

#include <CGAL/CLN/cl_integer.h>
#include <CGAL/CLN/cl_rational.h>
#include <CGAL/Quotient.h>

#include <cl_output.h> // for cl_default_print_flags

int main()
{
    using std::cout;
    using std::cin;
    using std::endl;

    bool ok;
    const char * filename = "test_io.tmp";
    cl_I pp, p = "1234567890123456";
    cl_I qq, q = "9876543210";
    cl_RA rr, r = cl_RA(p)/cl_RA(q);
    CGAL::Quotient<cl_I> qqII, qI (p, q);
    ofstream pout;
    ifstream pin;

    // Check that this is the correct base (10) for printing.
    cout << "Print in base " << cl_default_print_flags.rational_base << endl;
    assert(cl_default_print_flags.rational_base == 10);

    // Output to a file.
    pout.open(filename);
    pout << p << endl << q << endl << r << endl << qI << endl;
    cout << p << endl << q << endl << r << endl << qI << endl;
    pout.close();

    // Read from that file.
    pin.open(filename);
    pin >> pp >> qq >> rr >> qqII;
    pin.close();

    // Check read(print()) = id().
    ok = (p == pp) && (q == qq) && (r == rr) && (qI == qqII);

    ok = ok && CGAL::is_valid(p) && CGAL::is_valid(r);

    return (ok) ? 0 : -1;
}
#endif // CGAL_USE_CLN
