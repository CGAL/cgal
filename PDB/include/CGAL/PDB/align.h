#ifndef CGAL_DSRPDB_ALIGN_H
#define CGAL_DSRPDB_ALIGN_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/Transform.h>
#include <vector>

CGAL_PDB_BEGIN_NAMESPACE

//! \cond
inline double structal_score(const CGAL_PDB_NS::Point& a, 
			     const CGAL_PDB_NS::Point& b) {
  CGAL_PDB_NS::Squared_distance sd;
  return 20.0/(1.0+sd(a,b)*5.0);
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

  The value_types of ItA and ItB must be both
  be convertible to CGAL::PDB::Point.
*/
template <class ItA, class ItB>
Transform transform_taking_first_to_second(ItA pBegin,
					   ItA pEnd,
					   ItB qBegin, 
					   ItB qEnd) {
  CGAL_assertion(std::distance(pBegin, pEnd) == base_.size());
  // compute the centroid of the points and transform
  // pointsets so that their centroids coinside
      
  typedef typename std::iterator_traits<InputIteratorP>::value_type Point;
  typedef double RT;

  Vector center_p(0,0,0), center_q(0,0,0);
  int num_p = 0;
  int num_q = 0;
  for (InputIteratorP p_it = pBegin; p_it != pEnd; ++p_it) {
    //double x= p_it->x();
    center_p = center_p+Vector(p_it->x(), p_it->y(), p_it->z());//((*p_it)-CGAL::ORIGIN);
    num_p++;
  }
  center_p = center_p/num_p;
      
  for (InputIteratorQ q_it = qBegin; q_it != qEnd; ++q_it) {
    center_q = center_q + Vector(q_it->x(), q_it->y(), q_it->z()); //((*q_it)-CGAL::ORIGIN);
    num_q++;
  }
  center_q = center_q/num_q;
      
  CGAL_assertion(num_p == num_q);
      
  std::vector<Point> p_shifted, q_shifted;
  p_shifted.reserve(num_p);
  q_shifted.reserve(num_q);
  for (InputIteratorP p_it = pBegin; p_it != pEnd; ++p_it) {
    p_shifted.push_back((*p_it) - center_p);
  }
  for (InputIteratorQ q_it = qBegin; q_it != qEnd; ++q_it) {
    q_shifted.push_back((*q_it) - center_q);
  }
      
      
  // covariance matrix
  CGAL_TNT_NS::Array2D<RT> H(3, 3);
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

  CGAL_JAMA_NS::SVD<RT> svd(H);
  CGAL_TNT_NS::Array2D<RT> U(3, 3), V(3, 3);
  svd.getU(U);
  svd.getV(V);

  // the rotation matrix is R = VU^T
  CGAL_TNT_NS::Array2D<RT> UT = transpose(U);
  CGAL_TNT_NS::Array2D<RT> rot(3, 3);
  rot = matmult(V, UT);

  // check for reflection
  if (det(rot) < 0) {
    CGAL_TNT_NS::Array2D<RT> VT = transpose(V);
    CGAL_TNT_NS::Array2D<RT> UVT = matmult(U, VT);
    CGAL_TNT_NS::Array2D<RT> S(3, 3);
    S[0][0] = S[1][1] = 1;
    S[2][2] = det(UVT);
    S[0][1] = S[0][2] = S[1][0] = S[1][2] = S[2][0] = S[2][1] = 0;
    rot = matmult(matmult(U, S), VT);
  }
  
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



CGAL_PDB_END_NAMESPACE
#endif
