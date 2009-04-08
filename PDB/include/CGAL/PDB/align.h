#ifndef CGAL_DSRPDB_ALIGN_H
#define CGAL_DSRPDB_ALIGN_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d.h>
#include <CGAL/PDB/internal/tnt/jama_svd.h>
#include <CGAL/PDB/Matrix.h>
#include <CGAL/PDB/range.h>

namespace CGAL { namespace PDB {

//! \cond
inline double structal_score(const CGAL::PDB::Point& a, 
			     const CGAL::PDB::Point& b) {
  return 20.0/(1.0+squared_distance(a,b)*5.0);
}

inline double gap_score() {
  return -10.0;
}

struct DpP {
  double first;
  int second;
  DpP(double d, int s): first(d), second(s){}
  DpP(): first(-1), second(-3){}
};

inline std::ostream &operator<<(std::ostream &o, DpP p){
  return o << "(" << p.first << " " << p.second << ")";
}
//! \endcond


/*!  This computes the optimal rigid transform minimizing the least
  squares distance between the point sets defined by the two
  ranges. The two point sets must have equal size (since the points
  are taken to be in one-to-one correspondence.

  The value_types of RA and RB must be both
  be convertible to CGAL::PDB::Point. RA and RB must both be
  Boost.Range ranges.
*/
template <class RP, class RQ>
Transform transform_taking_first_to_second(RP rp, RQ rq){
  CGAL_precondition(CGAL::PDB::distance(rp)== CGAL::PDB::distance(rq));
  CGAL_precondition(CGAL::PDB::distance(rp) != 0);
  // compute the centroid of the points and transform
  // pointsets so that their centroids coinside
      
  typedef typename std::iterator_traits<typename RP::iterator>::value_type Point;
  typedef double RT;

  Vector center_p(0,0,0), center_q(0,0,0);
  int num_p = 0;
  int num_q = 0;
  CGAL_PDB_FOREACH (Point p, rp){
    //double x= p_it->x();
    center_p = center_p+(p- ORIGIN);//((*p_it)-CGAL::ORIGIN);
    num_p++;
  }
  center_p = center_p/num_p;
      
  CGAL_PDB_FOREACH(Point q, rq) {
    center_q = center_q + (q-ORIGIN);
    num_q++;
  }
  center_q = center_q/num_q;
  std::cout << center_p << " " << center_q << std::endl;
  CGAL_assertion(num_p == num_q);
      
  std::vector<Point> p_shifted, q_shifted;
  p_shifted.reserve(num_p);
  q_shifted.reserve(num_q);
  CGAL_PDB_FOREACH (Point p, rp){
    p_shifted.push_back(p - center_p);
  }
  CGAL_PDB_FOREACH(Point q, rq) {
    q_shifted.push_back(q - center_q);
  }
      
      
  // covariance matrix
  CGAL::PDB::TNT::Array2D<RT> H(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      H[i][j] = 0;
    }
  }
  for (int i = 0; i < num_p; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
	H[j][k] += p_shifted[i][j]*q_shifted[i][k];
      }
    }
  }

  CGAL::PDB::JAMA::SVD<RT> svd(H);
  CGAL::PDB::TNT::Array2D<RT> U(3, 3), V(3, 3);
  svd.getU(U);
  svd.getV(V);

  // the rotation matrix is R = VU^T
  CGAL::PDB::TNT::Array2D<RT> UT = transpose(U);
  CGAL::PDB::TNT::Array2D<RT> rot(3, 3);
  rot = matmult(V, UT);

  // check for reflection
  if (det(rot) < 0) {
    CGAL::PDB::TNT::Array2D<RT> VT = transpose(V);
    CGAL::PDB::TNT::Array2D<RT> UVT = matmult(U, VT);
    CGAL::PDB::TNT::Array2D<RT> S(3, 3);
    S[0][0] = S[1][1] = 1;
    S[2][2] = det(UVT);
    S[0][1] = S[0][2] = S[1][0] = S[1][2] = S[2][0] = S[2][1] = 0;
    rot = matmult(matmult(U, S), VT);
  }
  std::cout << rot << std::endl;
  Transform xf(rot[0][0], rot[0][1], rot[0][2],
	       rot[1][0], rot[1][1], rot[1][2],
	       rot[2][0], rot[2][1], rot[2][2]);
  Vector translation=center_q - xf(center_p);
  Transform trans(CGAL::Translation(), translation);
  return xf*trans;
}

class Transform_to_points{
  std::vector<Point> base_;
public:
  template <class It>
  Transform_to_points(It b, It e): base_(b,e) {
  }
  template <class It>
  Transform operator()(It b, It e) {
    return transform_taking_first_to_second(b,e, base_.begin(), base_.end());
  }
};



}}
#endif
