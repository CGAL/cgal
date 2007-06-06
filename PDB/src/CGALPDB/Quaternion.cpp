#include <CGAL/PDB/Quaternion.h>
#include <CGAL/PDB/Transform.h>
CGAL_PDB_BEGIN_NAMESPACE

/*Quaternion::Quaternion(const double m[3][3]){

  }*/

Quaternion::Quaternion(Transform t) {
  int            i;
  NT        qs2, qx2, qy2, qz2, tmp;
  NT        norm;
            
  for (i=0; i<4; ++i){
    q_[i] = 0;
  }
      
  //first compute squared magnitudes of q components - at least one
  //will be greater than 0 since q is unit magnitude
      
  qs2 = 0.25 * (t.m(0,0) + t.m(1,1) 
		+ t.m(2,2) + 1);
  qx2 = qs2 - 0.5 * (t.m(1,1) + t.m(2,2));
  qy2 = qs2 - 0.5 * (t.m(2,2) + t.m(0,0));
  qz2 = qs2 - 0.5 * (t.m(0,0) + t.m(1,1));
      
  //compute signed q components using numerically stable method
  if ((qs2 >= qx2) & (qs2 >= qy2) & (qs2 >= qz2)){
    q_[0] = sqrt(qs2);
    tmp = 0.25 / q_[0];
    q_[1] = (t.m(2,1) - t.m(1,2)) * tmp;
    q_[2] = (t.m(0,2) - t.m(2,0)) * tmp;
    q_[3] = (t.m(1,0) - t.m(0,1)) * tmp;
  } else if ((qx2 >= qs2) & (qx2 >= qy2) & (qx2 >= qz2)){
    q_[1] = sqrt(qx2);
    tmp = 0.25 / q_[1];
    q_[0] = (t.m(2,1) - t.m(1,2)) * tmp;
    q_[2] = (t.m(0,1) + t.m(1,0)) * tmp;
    q_[3] = (t.m(0,2) + t.m(2,0)) * tmp;
  } else if ((qy2 >= qs2) & (qy2 >= qx2) & (qy2 >= qz2)){
    q_[2] = sqrt(qy2);
    tmp = 0.25 / q_[2];
    q_[0] = (t.m(0,2) - t.m(2,0)) * tmp;
    q_[3] = (t.m(1,2) + t.m(2,1)) * tmp;
    q_[1] = (t.m(1,0) + t.m(0,1)) * tmp;
  } else if ((qz2 >= qs2) & (qz2 >= qx2) & (qz2 >= qy2)){
    q_[3] = sqrt(qz2);
    tmp = 0.25 / q_[3];
    q_[0] = (t.m(1,0) - t.m(0,1)) * tmp;
    q_[1] = (t.m(2,0) + t.m(0,2)) * tmp;
    q_[2] = (t.m(2,1) + t.m(1,2)) * tmp;
  }
      
  // for consistency, force positive scalar component [ (s; v) = (-s; -v) ]
  if (q_[0] < 0){
    for (i=0; i<4; ++i){
      q_[i] = -q_[i];
    }
  }
      
  norm = sqrt(q_[0]*q_[0] + q_[1]*q_[1]
	      + q_[2]*q_[2] + q_[3]*q_[3]);
  for (i=0; i<4; ++i){
    q_[i] /= norm;
  }
      
  //the first element we can uniquely derive from others
  // part_quarternion = Vector3D(quaternion[1], quaternion[2], quaternion[3]);
  //for (i=0; i<3; ++i){
  //part_quarternion[i] =
}

Transform Quaternion::transform() const {
  return transform(Vector(0,0,0));
}

Transform Quaternion::transform(Vector translation) const {
  NT q00 = operator[](0)*operator[](0);
  NT q11 = operator[](1)*operator[](1);
  NT q22 = operator[](2)*operator[](2);
  NT q33 = operator[](3)*operator[](3);
  NT q03 = operator[](0)*operator[](3);
  NT q13 = operator[](1)*operator[](3);
  NT q23 = operator[](2)*operator[](3);
  NT q02 = operator[](0)*operator[](2);
  NT q12 = operator[](1)*operator[](2);
  NT q01 = operator[](0)*operator[](1);

  //TNT::Matrix<NT> m(4,4);
  NT rot[3][3];
  rot[0][0] = q00 + q11 - q22 - q33;
  rot[1][1] = q00 - q11 + q22 - q33;
  rot[2][2] = q00 - q11 - q22 + q33;
  rot[0][1] = 2.0*(q12-q03);
  rot[1][0] = 2.0*(q12+q03);
  rot[0][2] = 2.0*(q13+q02);
  rot[2][0] = 2.0*(q13-q02);
  rot[1][2] = 2.0*(q23-q01);
  rot[2][1] = 2.0*(q23+q01);
      
  Transform ret(rot[0][0], rot[0][1], rot[0][2],translation[0],
		rot[1][0], rot[1][1], rot[1][2],translation[1],
		rot[2][0], rot[2][1], rot[2][2],translation[2]);
		
  return ret;
}

CGAL_PDB_END_NAMESPACE
