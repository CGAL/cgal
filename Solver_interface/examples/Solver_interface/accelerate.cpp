#include <CGAL/Accelerate_sparse_matrix.h>

int main()
{
  {
    CGAL::Accelerate_sparse_matrix<double> sm(3);

    sm.set_coef(0,0,1);
    sm.set_coef(0,1,2);
    sm.set_coef(0,2,3);
    sm.set_coef(1,1,5);
    sm.set_coef(2,2,9);

    sm.assemble_matrix();

    assert(sm.get_coef(1,0) == 0);
  }
  std::cout << std::endl;
  {
    CGAL::Accelerate_sparse_matrix<double> sm(3, true);// symmetric only set lower left

    sm.set_coef(0,0,1);
    sm.set_coef(1,0,2);
    sm.set_coef(2,0,3);
    sm.set_coef(1,1,5);
    sm.set_coef(2,2,9);

    sm.assemble_matrix();

    assert(sm.get_coef(1,0) == sm.get_coef(0, 1));
  }

  return 0;
}
