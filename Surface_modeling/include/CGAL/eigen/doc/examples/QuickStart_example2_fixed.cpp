#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
  Matrix3f m = Matrix3f::Random();
  m = (m + Matrix3f::Constant(1.2)) * 50;
  cout << "m =" << endl << m << endl;
  Vector3f v(1,2,3);
  
  cout << "m * v =" << endl << m * v << endl;
}
