#include <CGAL/basic.h>
#include <iostream>

#ifndef CGAL_USE_CLN
int main ()
{
  std::cout << "CGAL was installed with no CLN support." << std::endl;
  return 0;
}
#else // CGAL_USE_CLN

#include <CGAL/CLN/cl_integer.h>

int
main ()
{
    using std::cout;
    using std::endl;

    cl_I b = "31";
    cl_I p = b * 75;
    bool ok = (p==2325) && (CGAL_NTS compare(p,cl_I(31*75))==CGAL::EQUAL);

    cout << (ok ? "ok" : "ERROR") << endl;
    cout << p << " == " << 31*75 << endl;
    return (ok) ? 0 : -1;
}

#endif // CGAL_USE_CLN
