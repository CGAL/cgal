#ifndef CGAL_USE_CLN
#  error "You need to have CLN installed in order to run this example"
#endif

// This program tests some basic IO functions for some CLN types.

#include <CGAL/basic.h>
#include <iostream>
#include <fstream.h>
#include <cassert>

#include <CGAL/CLN/cl_integer.h>
#include <CGAL/CLN/cl_rational.h>
#include <CGAL/Quotient.h>

#include <cl_io.h>
#include <cl_integer_io.h>
#include <cl_rational_io.h>
#include <cl_output.h>

// #include <cl_input.h>
// extern cl_read_flags cl_I_read_flags;

int
main ()
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

    // cout << "Enter integer  : ";
    // cin >> p;
    // cout << "fprintdecimal: "; fprintdecimal(cout, p); cout << endl;
    // cout << "fprintbinary : "; fprintbinary(cout, p);  cout << endl;
    // cout << "fprint       : "; fprint(cout, p);        cout << endl;
    // cout << "print_integer: ";
    // print_integer(cout, cl_default_print_flags, p);
    // cout << endl;
    // cout << "print_integer: "; print_integer(cout, 10, p); cout << endl;
    // cout << "read value   : " << p << endl;
    // cout << "Calling read_integer() : ";
    // p = read_integer(cin, cl_I_read_flags);
    // cout << "Enter rational : ";
    // cin >> r;
    // cout << "read value: " << r << endl;
    // cout << "Enter Quotient<cl_I> : ";
    // cin >> qI;
    // cout << "read value: " << qI << endl;

    return (ok) ? 0 : -1;
}
