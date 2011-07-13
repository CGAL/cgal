#ifndef _MSC_VER
  #warning deprecated
#endif
/*
Matrix3d m = Matrix3d::Zero();
m.part<Eigen::StrictlyUpperTriangular>().setOnes();
cout << "Here is the matrix m:" << endl << m << endl;
cout << "And let us now compute m*m.adjoint() in a very optimized way" << endl
     << "taking advantage of the symmetry." << endl;
Matrix3d n;
n.part<Eigen::SelfAdjoint>() = (m*m.adjoint()).lazy();
cout << "The result is:" << endl << n << endl;
*/