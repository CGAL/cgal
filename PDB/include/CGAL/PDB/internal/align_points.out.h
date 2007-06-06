#ifndef CGAL_DSR_ALIGN_POINTS_H_
#define CGAL_DSR_ALIGN_POINTS_H_
#include <CGAL/PDB/basic.h>
#include <tnt/tnt_cmat.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d_utils.h>
#include <CGAL/PDB/internal/tnt/jama_svd.h>
#include <cassert>
#include <vector>
#include <iterator>

CGAL_PDB_BEGIN_INTERNAL_NAMESPACE


  template <class T>
  T det(const CGAL_TNT_NS::Array2D<T>& m) {
    assert(m.dim1() == 3);
    assert(m.dim2() == 3);

    return (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
            m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
            m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]));
  }

  // Transpose a matrix
  template <class T>
  CGAL_TNT_NS::Array2D<T> transpose(const CGAL_TNT_NS::Array2D<T>& m) {
    CGAL_TNT_NS::Array2D<T> mt(m.dim2(), m.dim1());
    for (int i = 0; i < m.dim1(); i++) {
      for (int j = 0; j < m.dim2(); j++) {
        mt[j][i] = m[i][j];
      }
    }
    return mt;
  }

  template <class InputIterator>
  CGAL_PDB_NS::Transform transform_taking_first_to_second(InputIterator pBegin, InputIterator pEnd,
						     InputIterator qBegin, InputIterator qEnd) {
    
    // compute the centroid of the points and transform
    // pointsets so that their centroids coinside
      
    typedef typename InputIterator::value_type Point;
    typedef double RT;

    CGAL_PDB_NS::Vector center_p(0,0,0), center_q(0,0,0);
    int num_p = 0;
    int num_q = 0;
    for (InputIterator p_it = pBegin; p_it != pEnd; ++p_it) {
      //double x= p_it->x();
      center_p = center_p+CGAL_PDB_NS::Vector(p_it->x(), p_it->y(), p_it->z());//((*p_it)-CGAL::ORIGIN);
      num_p++;
    }
    center_p = center_p/num_p;
      
    for (InputIterator q_it = qBegin; q_it != qEnd; ++q_it) {
      center_q = center_q + CGAL_PDB_NS::Vector(q_it->x(), q_it->y(), q_it->z()); //((*q_it)-CGAL::ORIGIN);
      num_q++;
    }
    center_q = center_q/num_q;
      
    assert(num_p == num_q);
      
    std::vector<Point> p_shifted, q_shifted;
    p_shifted.reserve(num_p);
    q_shifted.reserve(num_q);
    for (InputIterator p_it = pBegin; p_it != pEnd; ++p_it) {
      p_shifted.push_back((*p_it) - center_p);
    }
    for (InputIterator q_it = qBegin; q_it != qEnd; ++q_it) {
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

    CGAL_PDB_NS::Transform xf(rot, Point(0,0,0));
   
    xf.set_translation(center_q - xf(center_p));
   
    return xf;
  }


  inline double structal_score(const dsrpdb::Point& a, const dsrpdb::Point& b) {
    dsrpdb::Squared_distance sd;
    return 20.0/(1.0+sd(a,b)/5.0);
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

  template <class ItP, class ItM>
  void greedy_matching(ItP b_0, ItP e_0,
		       ItP b_1, ItP e_1,
		       ItM b_m, ItM e_m) {
    const double gap_cost=10.0;
    unsigned int d_0= std::distance(b_0, e_0);
    unsigned int d_1= std::distance(b_1, e_1);
    //typedef std::pair<double, int> DpP;
    CGAL_TNT_NS::Array2D<DpP> dp_matrix(d_0+2, d_1+2);


    std::vector<int> best_i(d_1+2);
    int best_j(0);

    for (unsigned int i=0; i< d_0+1; ++i){
      dp_matrix[i][0]=DpP(-gap_cost*i,0);
      best_i[i]= 0;
    }
    best_i.back()= 0;

    for (unsigned int i=0; i< d_1+1; ++i){
      dp_matrix[0][i]=DpP(-gap_cost*i,0);
    }

    dp_matrix[0][d_1+1]=DpP(-gap_cost*d_1, 0);
    dp_matrix[d_0+1][0]=DpP(-gap_cost*d_0, 0);
    
    {
      ItP c_0=b_0;
      for (unsigned int i=1; i< d_0+2; ++i, ++c_0){
	ItP c_1= b_1;
	best_j= 0;
	for (unsigned int j=1; j< d_1+2; ++j, ++c_1){
	  double ss;
	  if (c_0== e_0 || c_1==e_1) {
	    ss= 0;
	  } else {
	    ss= structal_score(*c_0, *c_1);
	  }
	  
	  
	  double v_d=  dp_matrix[i-1][j-1].first;
	  double v_i= -gap_cost*best_i[j-1] + dp_matrix[i-best_i[j-1]-1][j-1].first;
	  double v_j= -gap_cost*best_j + dp_matrix[i-1][j-best_j-1].first;

	  if (v_d >= v_i && v_d >= v_j){
	    dp_matrix[i][j]= DpP(ss+v_d, 0);
	    best_i[j-1]=0;
	    best_j=0;
	  } else if (v_i >= v_j){
	    dp_matrix[i][j] = DpP(ss+v_i, best_i[j-1]);
	    assert(v_i > v_d);
	    if (v_d >= v_j) best_j=0;
	  } else {
	    dp_matrix[i][j] = DpP(ss+v_j, -best_j);
	    if (v_d >= v_i) best_i[j-1]=0;
	  }

	  const double eps=.00001;
	  
	  for (unsigned int k=0; k< i; ++k){
	    assert(dp_matrix[i-k-1][j-1].first -gap_cost*k +ss <= dp_matrix[i][j].first + eps);
	    if (dp_matrix[i][j].second == static_cast<int>(k)) {
	      assert(std::abs(dp_matrix[i-k-1][j-1].first - gap_cost*k + ss - dp_matrix[i][j].first) < eps);
	    }
	  }

	  for (unsigned int k=0; k< j; ++k){
	    assert(dp_matrix[i-1][j-k-1].first -gap_cost*k + ss <= dp_matrix[i][j].first + eps);
	    if (dp_matrix[i][j].second == -static_cast<int>(k)) {
	      assert( std::abs(ss - gap_cost*k + dp_matrix[i-1][j-k-1].first - dp_matrix[i][j].first) < eps);
	    }
	  }

	  ++best_i[j-1];
	  ++best_j;
	}
      }
    }

    /*std::cout.precision(2);
      std::cout << dp_matrix << std::endl;*/

    // backtrace
    {
      unsigned int i=d_0+1, j= d_1+1;
     
      ItM c_m= e_m; //--c_m;
      //bool set=false;
      while (i >0 && j > 0) {
	if (dp_matrix[i][j].second ==0){

	} else if (dp_matrix[i][j].second < 0) {
	  j-= std::abs(dp_matrix[i][j].second);
	} else {
	  for (int k =0; k< dp_matrix[i][j].second; ++k){
	    --c_m;
	    *c_m=-1;
	  }
	  i-= dp_matrix[i][j].second;
	}
	--i;
	--j;

	if (i==0 || j==0) break;
	assert(c_m != b_m);
	--c_m;
	assert(j < d_1+1);
	assert(i < d_0+1);
	*c_m=j-1;
	//std::cout << i << " " << j << std::endl;
      }
    }
  }



  template <class ItP, class OItM>
  CGAL_PDB_NS::Transform refine_alignment(ItP b_0, ItP e_0,
				     ItP b_1, ItP e_1,
				     double motion,
				     OItM m) {
    
    std::vector<CGAL_PDB_NS::Point> p_0(b_0, e_0), p_1(b_1, e_1);
    double dist= std::numeric_limits<double>::infinity();
    CGAL_PDB_NS::Squared_distance sd;

    typedef std::vector<int>::const_iterator ItM;
    std::vector<int> cur_matching(std::distance(b_0, e_0), -1);
    CGAL_PDB_NS::Transform tr;
    do {
      greedy_matching(p_0.begin(), p_0.end(),
		      p_1.begin(), p_1.end(),
		      cur_matching.begin(), cur_matching.end());

      /*std::copy(cur_matching.begin(), cur_matching.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;*/

      std::vector<CGAL_PDB_NS::Point> cp_0, cp_1;
      {
	ItP c_0= b_0;
	for (ItM c_m= cur_matching.begin(); c_m != cur_matching.end(); ++c_m, ++c_0 ){
	  if (*c_m != -1){
	    cp_0.push_back(*c_0);
	    cp_1.push_back(p_1[*c_m]);
	  }
	}
      }
      
      tr = transform_taking_first_to_second(cp_0.begin(),cp_0.end(),
					cp_1.begin(), cp_1.end());
      //std::cout << tr << std::endl;

      double sum_dist=0;
      unsigned int i=0;
      for (ItP c_0=b_0; c_0 != e_0; ++c_0, ++i){
	dsrpdb::Point pn= tr(*c_0);
	sum_dist+= sd(pn, p_0[i]);
	p_0[i]=pn;
      }

      dist=std::sqrt(sum_dist)/std::distance(b_0, e_0);

      std::cout << "Moved " << dist << std::endl;

    } while (dist > motion);
    std::copy(cur_matching.begin(), cur_matching.end(), m);
    return tr;
  }

CGAL_PDB_END_INTERNAL_NAMESPACE
#endif
