#ifndef _MSC_VER
  #warning deprecated
#endif
/*
Matrix3d m = Matrix3d::Zero();
m.part<Eigen::UpperTriangular>().setOnes();
cout << "Here is the matrix m:" << endl << m << endl;
Matrix3d n = Matrix3d::Ones();
n.part<Eigen::LowerTriangular>() *= 2;
cout << "Here is the matrix n:" << endl << n << endl;
cout << "And now here is m.inverse()*n, taking advantage of the fact that"
        " m is upper-triangular:" << endl
     << m.marked<Eigen::UpperTriangular>().solveTriangular(n);
*/