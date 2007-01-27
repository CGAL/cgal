// Tests that all mixed operators on Quotient are defined.
// Sylvain Pion

#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

int main()
{
  typedef CGAL::MP_Float     RT;
  typedef CGAL::Quotient<RT> QT;

  RT r(1);
  QT q(1);

  q+q; q+r; r+q; q+1; 1+q;
  q-q; q-r; r-q; q-1; 1-q;
  q*q; q*r; r*q; q*1; 1*q;
  q/q; q/r; r/q; q/1; 1/q;
  -q;
  q<q; q<r; r<q; q<1; 1<q;
  q>q; q>r; r>q; q>1; 1>q;
  q<=q; q<=r; r<=q; q<=1; 1<=q;
  q>=q; q>=r; r>=q; q>=1; 1>=q;
  q==q; q==r; r==q; q==1; 1==q;
  q!=q; q!=r; r!=q; q!=1; 1!=q;

  return 0;
}
