#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
  MatrixXf m = MatrixXf::Random(3,3);
  m = (m + MatrixXf::Constant(3,3,1.2)) * 50;
  cout << "m =" << endl << m << endl;
  VectorXf v(3);
  v << 1, 2, 3;
  cout << "m * v =" << endl << m * v << endl;
}
