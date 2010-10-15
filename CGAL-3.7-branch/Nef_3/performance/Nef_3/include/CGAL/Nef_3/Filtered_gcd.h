#ifndef CGAL_FILTERED_GCD_H
#define CGAL_FILTERED_GCD_H

#include <CGAL/basic.h>
#include <CGAL/Lazy_exact_nt.h>
#undef _DEBUG
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

/*
template <typename NT> 
Lazy_exact_nt<NT> operator%=(const Lazy_exact_nt<NT>& a, const Lazy_exact_nt<NT>& b) {

  NT res = a.exact();
  res %= b.exact();
  return Lazy_exact_nt<NT>(res);
}
*/
template <typename NT> 
Lazy_exact_nt<NT> gcd(const Lazy_exact_nt<NT>& a, const Lazy_exact_nt<NT>& b) {

  NT nta = a.exact();
  NT ntb = b.exact();
  return Lazy_exact_nt<NT>(CGAL_NTS gcd(nta,ntb));
}

} //namespace CGAL
#endif // CGAL_FILTERED_GCD_H
