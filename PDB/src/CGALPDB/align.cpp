#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/iterator.h>
#include <CGAL/PDB/align.h>
#include <cassert>


CGAL_PDB_BEGIN_NAMESPACE

#if 0

  
/*!  This function computes the optimal fixed matching between the
  two point sets defined by the two ranges. The matching is put into
  the range (b_m...e_m) which must have the same size as
  (b_0...e_0). The matching is represented as an iterator in the
  range (b_1...e_1) for each element in (b_0...e_0) with e_1 used to
  designate that there is no matching (i.e. that there is a gap).

  The value_types of It0 and It1 must both be convertible to
  CGAL::PDB::Point.

  The run time is the product of distance(b_0,e_0)*distance(b_1,e_1).

  Return the structal score.
*/
template <class It0, class It1, class ItM>
double optimal_fixed_matching(It0 b_0, It0 e_0,
			      It1 b_1, It1 e_1,
			      ItM b_m, ItM e_m) {
  using namespace CGAL_PDB_INTERNAL_NS::Align;
  assert(std::distance(b_0, e_0) == std::distance(b_m, e_m));
  const double gap_cost=10.0;
  unsigned int d_0= std::distance(b_0, e_0);
  unsigned int d_1= std::distance(b_1, e_1);
  //typedef std::pair<double, int> DpP;
  CGAL_TNT_NS::Array2D<DpP> dp_matrix(d_0+2, d_1+2);


  std::vector<int> best_i(d_1+2);
  int best_j(0);

  for (unsigned int i=0; i< d_0+1; ++i){
    dp_matrix[i][0]=DpP(-gap_cost*i,0);
  }
  for (unsigned int j=0; j< d_1+2; ++j){
    best_i[j]= 0;
  }
  //best_i.back()= 0;

  for (unsigned int j=0; j< d_1+1; ++j){
    dp_matrix[0][j]=DpP(-gap_cost*j,0);
  }

  dp_matrix[0][d_1+1]=DpP(-gap_cost*d_1, 0);
  dp_matrix[d_0+1][0]=DpP(-gap_cost*d_0, 0);
    
  {
    It0 c_0=b_0;
    for (unsigned int i=1; i< d_0+2; ++i, ++c_0){
      It1 c_1= b_1;
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

	/*const double eps=.00001;
	  
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
	  }*/

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
    It1 c_1= e_1;
    ItM c_m= e_m; //--c_m;
    //int d= c_1-b_1;
    assert(c_1-b_1 == static_cast<int>(j)-1);

    //bool set=false;
    while (i >0 && j > 0) {
      if (dp_matrix[i][j].second ==0){

      } else if (dp_matrix[i][j].second < 0) {
	for (int k=0; k< std::abs(dp_matrix[i][j].second); ++k) --c_1;
	j-= std::abs(dp_matrix[i][j].second);
      } else {
	for (int k =0; k< dp_matrix[i][j].second; ++k){
	  --c_m;
	  *c_m=e_1;
	}
	i-= dp_matrix[i][j].second;
      }
      --i;
      --j;
      --c_1;

      if (i==0 || j==0) break;
      assert(c_m != b_m);
      --c_m;
      assert(j < d_1+1);
      assert(i < d_0+1);
      *c_m=c_1; //j-1;

      //int d= c_1-b_1;
      assert(c_1-b_1 == static_cast<int>(j)-1);
      //std::cout << i << " " << j << std::endl;
    }
  }
  return dp_matrix[d_0+1][d_1+1].first + gap_cost*(std::max  BOOST_PREVENT_MACRO_SUBSTITUTION(d_0, d_1)- std::min  BOOST_PREVENT_MACRO_SUBSTITUTION (d_0, d_1));
}



/*!  This function refines the alignment between the two point sets
  represented by the ranges by alternating dynamic programming based
  sequence alignement and then rigid transformation of the point
  set. The final correspondence is put into m as a sequence of
  iterators in the range (b_1...e_1). motion is the minimal amount
  that the points must move in an iteration for the process to
  continue (if they move less than that amount the function
  returns).

  The run time is the
  number_of_iterators*distance(b_0,e_0)*distance(b_1,e_1).
    
*/
template <class It0, class It1, class OItM>
Transform refine_alignment(It0 b_0, It0 e_0,
			   It1 b_1, It1 e_1,
			   double motion,
			   OItM m) {
    
  std::vector<Point> p_0(b_0, e_0);
  double dist= std::numeric_limits<double>::infinity();
  Squared_distance sd;

  typedef typename std::vector<It1>::const_iterator ItM;
  std::vector<It1> cur_matching(std::distance(b_0, e_0), e_1);
  Transform tr;

  double last_score=-std::numeric_limits<double>::infinity();
  do {
    double score = optimal_fixed_matching(p_0.begin(), p_0.end(),
					  b_1, e_1,
					  cur_matching.begin(), cur_matching.end());
    if (score < last_score) break;
    last_score=score;
    std::cout << score << std::endl;
    /*for (unsigned int i=0; i< cur_matching.size(); ++i){
      if (cur_matching[i] != e_1) {
      std::cout << cur_matching[i]- b_1 << " ";
      } else {
      std::cout << -1 << " ";
      }
      }
      std::cout << std::endl;*/

    std::vector<Point> cp_0, cp_1;
    {
      It0 c_0= b_0;
      for (ItM c_m= cur_matching.begin(); c_m != cur_matching.end(); ++c_m, ++c_0 ){
	if (*c_m != e_1){
	  cp_0.push_back(*c_0);
	  cp_1.push_back(**c_m);
	}
      }
    }
      
    tr = transform_taking_first_to_second(cp_0.begin(),cp_0.end(),
					  cp_1.begin(), cp_1.end());
    //std::cout << tr << std::endl;

    double sum_dist=0;
    unsigned int i=0;
    for (It0 c_0=b_0; c_0 != e_0; ++c_0, ++i){
      Point pn= tr(*c_0);
      sum_dist+= std::sqrt(sd(pn, p_0[i]));
      p_0[i]=pn;
    }

    dist=sum_dist/(std::distance(b_0, e_0));

    std::cout << "Moved " << dist << std::endl;

  } while (dist > motion);
    
  It1 c_1= b_1;
  int o_1=0;
  for (typename std::vector<It1>::const_iterator it= cur_matching.begin();
       it != cur_matching.end(); ++it){
    if (*it == e_1){
      *m= e_1;
      ++m;
    } else {
      while (o_1 < std::distance(b_1, *it)){
	++o_1;
	++c_1;
      }
      *m= c_1;
      ++m;
    }
  }

  return tr;
}


#endif

#if 0
void align_second_to_first(const Chain &base, 
				 Chain &o) {
  assert(base.number_of_monomers()==o.number_of_monomers());
   
  //bp= base.backbone();
  //op= o.backbone();
  //get_backbone_coordinates(base, std::back_inserter(bp));
  //get_backbone_coordinates(o, std::back_inserter(op));
  Transform tr
    = transform_taking_first_to_second(backbone_coordinates_begin(o),
				       backbone_coordinates_end(o),
				       backbone_coordinates_begin(base),
				       backbone_coordinates_end(base));
    
  if (true) {
    std::cout << tr << std::endl;
  }
     

  for (Chain::Monomer_iterator it= o.monomers_begin(); 
       it != o.monomers_end(); ++it){
    Monomer &aa= it->monomer();
    //dsr::vector<Residue::Atom_label> als= aa->atoms();
       
    for (Monomer::Atom_iterator it= aa.atoms_begin(); 
	 it != aa.atoms_end(); ++it){
      Atom a= it->atom();
      Point pt= Point(a.cartesian_coords().x(), a.cartesian_coords().y(),
		      a.cartesian_coords().z());
      Point tpt(tr(pt));
      a.set_cartesian_coords(Point(tpt.x(), tpt.y(), tpt.z()));

      it->atom()=a;

    }
    
  }
}
#endif

#if 0
Transform refine_alignment_of_second_to_first(const Chain &base,
						      Chain &o) {
  std::vector<Point> sheet_points, model_points;
  ca_coordinates(base.atoms_begin(), base.atoms_end(),
			 std::back_inserter(sheet_points));
  ca_coordinates(o.atoms_begin(), o.atoms_end(),
			 std::back_inserter(model_points));
  std::vector<std::vector<Point>::const_iterator> match;
  Transform t= refine_alignment(model_points.begin(), model_points.end(),
					sheet_points.begin(), sheet_points.end(),
					.01,
					std::back_inserter(match));
  for (Chain::Monomer_iterator rit= o.monomers_begin(); rit != o.monomers_end(); ++rit){
    for (Monomer::Atom_iterator ait= rit->monomer().atoms_begin(); 
	 ait != rit->monomer().atoms_end(); ++ait){
      ait->atom().set_cartesian_coords(t(ait->atom().cartesian_coords()));
    }
  }
  return t;
}
#endif

CGAL_PDB_END_NAMESPACE
