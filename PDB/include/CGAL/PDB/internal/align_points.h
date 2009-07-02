#ifndef CGAL_DSR_ALIGN_POINTS_H_
#define CGAL_DSR_ALIGN_POINTS_H_
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/internal/tnt/tnt_cmat.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/align.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d.h>
#include <CGAL/PDB/internal/tnt/tnt_array2d_utils.h>
#include <CGAL/PDB/internal/tnt/jama_svd.h>
#include <cassert>
#include <vector>
#include <iterator>
#include <boost/range.hpp>

namespace CGAL { namespace PDB { namespace internal {


 
 template <class Range0, class Range1, class RangeM>
  void greedy_matching(const Range0& r0,
                       const Range1& r1,
		       RangeM& rm) {
    const double gap_cost=10.0;
    unsigned int d_0= std::distance(r0.begin(), r0.end());
    unsigned int d_1= std::distance(r1.begin(), r1.end());
    CGAL_assertion(!r0.empty() && !r1.empty());
    //typedef std::pair<double, int> DpP;
    CGAL::PDB::TNT::Array2D<DpP> dp_matrix(d_0+2, d_1+2);


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
      typename Range0::const_iterator c_0=r0.begin();
      for (unsigned int i=1; i< d_0+2; ++i, ++c_0){
        typename Range1::const_iterator c_1=r1.begin();
	best_j= 0;
	for (unsigned int j=1; j< d_1+2; ++j, ++c_1){
	  double ss;
	  if (c_0== r0.end() || c_1==r1.end()) {
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
	    CGAL_assertion(v_i > v_d);
	    if (v_d >= v_j) best_j=0;
	  } else {
	    dp_matrix[i][j] = DpP(ss+v_j, -best_j);
	    if (v_d >= v_i) best_i[j-1]=0;
	  }

	  const double eps=.00001;
	  
	  for (unsigned int k=0; k< i; ++k){
	    CGAL_assertion(dp_matrix[i-k-1][j-1].first -gap_cost*k +ss <= dp_matrix[i][j].first + eps);
	    if (dp_matrix[i][j].second == static_cast<int>(k)) {
	      CGAL_assertion(std::abs(dp_matrix[i-k-1][j-1].first - gap_cost*k + ss - dp_matrix[i][j].first) < eps);
	    }
	  }

	  for (unsigned int k=0; k< j; ++k){
	    CGAL_assertion(dp_matrix[i-1][j-k-1].first -gap_cost*k + ss <= dp_matrix[i][j].first + eps);
	    if (dp_matrix[i][j].second == -static_cast<int>(k)) {
	      CGAL_assertion( std::abs(ss - gap_cost*k + dp_matrix[i-1][j-k-1].first - dp_matrix[i][j].first) < eps);
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
     
      typename RangeM::iterator c_m= rm.end(); //--c_m;
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
	CGAL_assertion(c_m != rm.begin());
	--c_m;
	CGAL_assertion(j < d_1+1);
	CGAL_assertion(i < d_0+1);
	*c_m=j-1;
      }
    }
  }



  template <class Range, class OItM>
  CGAL::PDB::Transform refine_alignment(Range r0,
                                        Range r1,
                                        double motion,
                                        OItM m) {
    CGAL_precondition(!r0.empty());
    CGAL_precondition(!r1.empty());
    std::vector<CGAL::PDB::Point> p_0(r0.begin(), r0.end()),
      p_1(r1.begin(), r1.end());
    double dist= std::numeric_limits<double>::infinity();

    typedef std::vector<int>::const_iterator ItM;
    std::vector<int> cur_matching(std::distance(r0.begin(), r0.end()), -1);
    CGAL::PDB::Transform tr;
    do {
      greedy_matching(p_0, p_1,
		      cur_matching);

      /*std::copy(cur_matching.begin(), cur_matching.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;*/

      std::vector<CGAL::PDB::Point> cp_0, cp_1;
      {
	typename Range::iterator c_0= r0.begin();
	for (ItM c_m= cur_matching.begin(); c_m != cur_matching.end(); ++c_m, ++c_0 ){
	  if (*c_m != -1){
	    cp_0.push_back(*c_0);
	    cp_1.push_back(p_1[*c_m]);
	  }
	}
      }
      
      tr = transform_taking_first_to_second(cp_0,
                                            cp_1);
      //std::cout << tr << std::endl;

      double sum_dist=0;
      unsigned int i=0;
      for (typename Range::iterator c_0=r0.begin(); c_0 != r0.end(); ++c_0, ++i){
	Point pn= tr(*c_0);
	double csd= squared_distance(pn, p_0[i]);
        std::cout << csd << " " << pn << " " << p_0[i] 
                  << " " << *c_0 << std::endl;
        CGAL_assertion(csd >=0);
        sum_dist+=csd;
	p_0[i]=pn;
      }

      dist=std::sqrt(sum_dist)/std::distance(r0.begin(), r0.end());

      std::cout << "Moved " << dist << " (" << sum_dist << ")" << std::endl;

    } while (dist > motion);
    std::copy(cur_matching.begin(), cur_matching.end(), m);
    return tr;
  }

}}}
#endif
