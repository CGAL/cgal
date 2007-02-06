#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/internal/tnt/tnt_cmat.h>

CGAL_PDB_BEGIN_NAMESPACE
Transform::Transform(const Vector &trans, const Quaternion &q) {
  NT q00 = q[0]*q[0];
  NT q11 = q[1]*q[1];
  NT q22 = q[2]*q[2];
  NT q33 = q[3]*q[3];
  NT q03 = q[0]*q[3];
  NT q13 = q[1]*q[3];
  NT q23 = q[2]*q[3];
  NT q02 = q[0]*q[2];
  NT q12 = q[1]*q[2];
  NT q01 = q[0]*q[1];

  //TNT::Matrix<NT> m(4,4);
  rot_[0][0] = q00 + q11 - q22 - q33;
  rot_[1][1] = q00 - q11 + q22 - q33;
  rot_[2][2] = q00 - q11 - q22 + q33;
  rot_[0][1] = 2.0*(q12-q03);
  rot_[1][0] = 2.0*(q12+q03);
  rot_[0][2] = 2.0*(q13+q02);
  rot_[2][0] = 2.0*(q13-q02);
  rot_[1][2] = 2.0*(q23-q01);
  rot_[2][1] = 2.0*(q23+q01);
      
  //TNT::Matrix<NT> mt=TNT::transpose(m);
  //TNT::Matrix<NT> mm= TNT::matmult(mt, m);
  //std::cout << mm << std::endl;
      
  //TNT::Matrix<NT> mat= q.to_matrix();
  /*for (unsigned int i=0; i< 3; ++i){
    for (unsigned int j=0; j<3; ++j){
    rot_[i][j]=m[i][j];
    }
    }*/
  trans_[0]=trans.x();
  trans_[1]=trans.y();
  trans_[2]=trans.z();
}

void Transform::write(std::ostream &out) const {
  for (unsigned int i=0; i< 3; ++i){
    out << rot_[i][0] << "\t" <<  rot_[i][1] << "\t" << rot_[i][2] << "\t" << trans_[i] << std::endl;
  }
  out << 0 << "\t" <<  0 << "\t" << 0 << "\t" << 1 << std::endl;

}
CGAL_PDB_END_NAMESPACE
