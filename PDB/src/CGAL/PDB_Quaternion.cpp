#include <CGAL/PDB/Quaternion.h>
CGAL_PDB_BEGIN_NAMESPACE

Quaternion::Quaternion(const double m[3][3]){
  int            i;
  NT        qs2, qx2, qy2, qz2, tmp;
  NT        norm;
            
  for (i=0; i<4; ++i){
    q_[i] = 0;
  }
      
  //first compute squared magnitudes of q components - at least one
  //will be greater than 0 since q is unit magnitude
      
  qs2 = 0.25 * (m[0][0] + m[1][1] 
		+ m[2][2] + 1);
  qx2 = qs2 - 0.5 * (m[1][1] + m[2][2]);
  qy2 = qs2 - 0.5 * (m[2][2] + m[0][0]);
  qz2 = qs2 - 0.5 * (m[0][0] + m[1][1]);
      
  //compute signed q components using numerically stable method
  if ((qs2 >= qx2) & (qs2 >= qy2) & (qs2 >= qz2)){
    q_[0] = sqrt(qs2);
    tmp = 0.25 / q_[0];
    q_[1] = (m[2][1] - m[1][2]) * tmp;
    q_[2] = (m[0][2] - m[2][0]) * tmp;
    q_[3] = (m[1][0] - m[0][1]) * tmp;
  } else if ((qx2 >= qs2) & (qx2 >= qy2) & (qx2 >= qz2)){
    q_[1] = sqrt(qx2);
    tmp = 0.25 / q_[1];
    q_[0] = (m[2][1] - m[1][2]) * tmp;
    q_[2] = (m[0][1] + m[1][0]) * tmp;
    q_[3] = (m[0][2] + m[2][0]) * tmp;
  } else if ((qy2 >= qs2) & (qy2 >= qx2) & (qy2 >= qz2)){
    q_[2] = sqrt(qy2);
    tmp = 0.25 / q_[2];
    q_[0] = (m[0][2] - m[2][0]) * tmp;
    q_[3] = (m[1][2] + m[2][1]) * tmp;
    q_[1] = (m[1][0] + m[0][1]) * tmp;
  } else if ((qz2 >= qs2) & (qz2 >= qx2) & (qz2 >= qy2)){
    q_[3] = sqrt(qz2);
    tmp = 0.25 / q_[3];
    q_[0] = (m[1][0] - m[0][1]) * tmp;
    q_[1] = (m[2][0] + m[0][2]) * tmp;
    q_[2] = (m[2][1] + m[1][2]) * tmp;
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
CGAL_PDB_END_NAMESPACE
