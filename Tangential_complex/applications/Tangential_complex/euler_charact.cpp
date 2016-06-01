#include <CGAL/Tangential_complex/Simplicial_complex.h>

using namespace CGAL::Tangential_complex_;

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Usage: euler_charact <input.off>\n";
    return 0;
  }
  
  Simplicial_complex complex;
  complex.load_simplices_from_OFF(argv[1]);

  complex.euler_characteristic(true);

  return 0;
}
