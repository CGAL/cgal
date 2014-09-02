#include <CGAL/Epick_d.h>
#include <CGAL/Tangential_complex.h>

int main()
{
  typedef CGAL::Epick_d<CGAL::Dimension_tag<3> > Kernel;
  const int INTRINSIC_DIMENSION = 2;

  CGAL::Tangential_complex<Kernel, INTRINSIC_DIMENSION> tc;

  return 0;
}