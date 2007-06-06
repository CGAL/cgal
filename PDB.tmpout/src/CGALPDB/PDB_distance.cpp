#include <CGAL/PDB/distance.h>
#include <iterator>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/align_generic.h>
#include <CGAL/PDB/iterator.h>
#include <CGAL/PDB/Model.h>
#include <CGAL/PDB/internal/Error_logger.h>

#ifdef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
// need std::copy
#include <algorithm>
#endif

CGAL_PDB_BEGIN_NAMESPACE

  struct Identity {
    const Point &operator()(const Point &a) const {
      return a;
    }
  };

  template <class It, class F>
  double no_align_cRMS(It ab, It ae, It bb, It be, F f){
    Squared_distance sd;
    if (std::distance(ab, ae) != std::distance(bb, be)){
      CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error("Protein chains used for computing cRMS must have equal lengths.\n");
    }
    double ret=0;
    int num=0;
    for (It bc= bb, ac= ab; bc != be; ++bc, ++ac){

      Point pt= Point(*bc);
      Point tpt= f(pt);
      ret += sd(*ac, tpt);
      ++num;
    }
    return std::sqrt(ret)/ num;
  }

  // base, o
  template <class It>
  double cRMS(It ab, It ae, It bb, It be){
    assert(std::distance(ab, ae) == std::distance(bb, be));

    CGAL::PDB::Transform tr
      = CGAL::PDB::transform_taking_first_to_second( bb, be, ab, ae);
    /*{
      std::cout << "First write\n";
      std::ofstream outtest("/tmp/foo.pdb");
      o.write_pdb(outtest);
      }*/

    if (false) {
      std::cout << tr << std::endl;
    }
     
    return no_align_cRMS(bb, be, ab, ae, tr);
   
  }

  double no_align_cRMS(const Protein &a, const Protein &b) {
    return no_align_cRMS(all_coordinates_begin(a), all_coordinates_end(a),
			 all_coordinates_begin(b), all_coordinates_end(b), Identity());
      /*double err=0;
    int n=0;
    assert(a.number_of_atoms() == b.number_of_atoms());
    for (Protein::Const_atoms_iterator ait= a.atoms_begin(), bit = b.atoms_begin(); 
	 ait != a.atoms_end(); ++ait, ++bit){
      ++n;
      Vector v= ait->second.cartesian_coords() - bit->second.cartesian_coords();
      err += v*v;
    }
    
    return std::sqrt(err)/n;*/
  }

  double no_align_ca_cRMS(const Protein &a, const Protein &b) {
    return no_align_cRMS(ca_coordinates_begin(a), ca_coordinates_end(a), 
			 ca_coordinates_begin(b), ca_coordinates_end(b), Identity());
  }

  double cRMS(const Protein &a, const Protein &b) {
    return cRMS(all_coordinates_begin(a), all_coordinates_end(a),
		all_coordinates_begin(b), all_coordinates_end(b));
  }

  double ca_cRMS(const Protein &a, const Protein &b) {
    return cRMS(ca_coordinates_begin(a), ca_coordinates_end(a), 
		ca_coordinates_begin(b), ca_coordinates_end(b));
  }

  template <class It>
  double dRMS(It ab, It ae, It bb, It be) {
    assert(std::distance(ab, ae) == std::distance(bb, be));
    double ret=0;
    int count=0;
    for (It ac= ab, bc= bb; ac != ae; ++ac, ++bc) {
      for (It ac2= ab, bc2= bb; ac2 != ac; ++ac2, ++bc2) {
	Vector va= *ac- *ac2;
	Vector vb= *bc- *bc2;
	double da= std::sqrt(va*va);
	double db= std::sqrt(vb*vb);
	ret+= (da-db)*(da-db);
	++count;
      }
    }
    return ret/count;
  }

  double dRMS(const Protein &a, const Protein &b) {
    std::vector<Point> aa;
    std::vector<Point> ba;
    
    coordinates(a.atoms_begin(), a.atoms_end(), std::back_inserter(aa));
    coordinates(b.atoms_begin(), b.atoms_end(), std::back_inserter(ba));
    return dRMS(aa.begin(), aa.end(), ba.begin(), ba.end());
  }

  double ca_dRMS(const Protein &a, const Protein &b){
    std::vector<Point> aa;
    std::vector<Point> ba;
    
    ca_coordinates(a.atoms_begin(), a.atoms_end(), std::back_inserter(aa));
    ca_coordinates(b.atoms_begin(), b.atoms_end(), std::back_inserter(ba));
    return dRMS(aa.begin(), aa.end(), ba.begin(), ba.end());
  }
  template <class It>
  Matrix distance_matrix(It ab, It ae){
    int dist= std::distance(ab, ae);
    Matrix ret(dist,dist);
    
    int indi=0;
    for (It ac= ab; ac != ae; ++ac, ++indi){
      int indj=0;
      for (It ac2= ab; ac2 != ac; ++ac2, ++indj){
	Vector dif= *ac- *ac2;
	//std::cout << dif << std::endl;
	double d= std::sqrt(dif*dif);
	//std::cout << d << std::endl;
	ret[indi][indj]=d;
	ret[indj][indi]=d;
      }
      ret[indi][indi]=0;
    }
    return ret;
  }

  Matrix distance_matrix(const Protein &a){
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
    std::vector<Point> aa(all_coordinates_begin(a), all_coordinates_end(a));
#else
    std::vector<Point> aa;
    std::copy(all_coordinates_begin(a), all_coordinates_end(a), std::back_inserter(aa));
#endif
    return distance_matrix(aa.begin(), aa.end());
  }

  Matrix ca_distance_matrix(const Protein &a) {
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
    std::vector<Point> aa(ca_coordinates_begin(a), ca_coordinates_end(a));
#else
    std::vector<Point> aa;
    std::copy(ca_coordinates_begin(a), ca_coordinates_end(a), std::back_inserter(aa));
#endif
    return distance_matrix(aa.begin(), aa.end());
  }
  Matrix backbone_distance_matrix(const Protein &a) {
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
    std::vector<Point> aa(backbone_coordinates_begin(a), backbone_coordinates_end(a));
#else
    std::vector<Point> aa;
    std::copy(backbone_coordinates_begin(a), backbone_coordinates_end(a), std::back_inserter(aa));
#endif
    return distance_matrix(aa.begin(), aa.end());
  }

  Matrix distance_matrix(const Model &a){
    std::vector<Point> aa;
    for (unsigned int j=0; j< a.number_of_chains(); ++j){
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
      aa.insert(aa.end(), all_coordinates_begin(a.chain(j)), all_coordinates_end(a.chain(j)));
#else
      std::copy(all_coordinates_begin(a.chain(j)), all_coordinates_end(a.chain(j)), std::back_inserter(aa));
#endif
    }
    return distance_matrix(aa.begin(), aa.end());
  }

  Matrix ca_distance_matrix(const Model &a) {
    std::vector<Point> aa;
    for (unsigned int j=0; j< a.number_of_chains(); ++j){
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
      aa.insert(aa.end(), ca_coordinates_begin(a.chain(j)), ca_coordinates_end(a.chain(j)));
#else
      std::copy(ca_coordinates_begin(a.chain(j)), ca_coordinates_end(a.chain(j)), std::back_inserter(aa));
#endif
    }
    return distance_matrix(aa.begin(), aa.end());
  }

  Matrix backbone_distance_matrix(const Model &a) {
    std::vector<Point> aa;
    for (unsigned int j=0; j< a.number_of_chains(); ++j){
#ifndef CGAL_CFG_MISSING_TEMPLATE_VECTOR_CONSTRUCTORS_BUG
      aa.insert(aa.end(), backbone_coordinates_begin(a.chain(j)), backbone_coordinates_end(a.chain(j)));
#else
      std::copy(backbone_coordinates_begin(a.chain(j)), backbone_coordinates_end(a.chain(j)), std::back_inserter(aa));
#endif
    }
    return distance_matrix(aa.begin(), aa.end());
  }


CGAL_PDB_END_NAMESPACE
