#include <dsrpdb/config.h>
#include <dsrpdb/align.h>
#include <dsrpdb/Point.h>
#include <dsrpdb/Transform.h>
#include <dsrpdb_internal/align_points.h>
#include <cassert>
#include <dsrpdb/iterator.h>

namespace dsrpdb {

 

  Transform compute_transform_taking_first_to_second(const std::vector<Point> &a,
						     const std::vector<Point> &b) {
    Transform tr
      = dsrpdb_internal::transform_taking_first_to_second(a.begin(),
						   a.end(),
						   b.begin(),
						   b.end());
    return tr;
  }


  void align_second_protein_to_first(const dsrpdb::Protein &base, 
				     dsrpdb::Protein &o) {
    assert(base.number_of_residues()==o.number_of_residues());
   
    //bp= base.backbone();
    //op= o.backbone();
    //get_backbone_coordinates(base, std::back_inserter(bp));
    //get_backbone_coordinates(o, std::back_inserter(op));
    Transform tr
      = dsrpdb_internal::transform_taking_first_to_second(backbone_coordinates_begin(o),
							  backbone_coordinates_end(o),
							  backbone_coordinates_begin(base),
							  backbone_coordinates_end(base));
    
    if (false) {
      std::cout << tr << std::endl;
    }
     

    for (dsrpdb::Protein::Residues_iterator it= o.residues_begin(); it != o.residues_end(); ++it){
      dsrpdb::Residue &aa= *it;
      //dsr::vector<Residue::Atom_label> als= aa->atoms();
       
      for (dsrpdb::Residue::Atoms_iterator it= aa.atoms_begin(); it != aa.atoms_end(); ++it){
	dsrpdb::Atom a= it->second;
	Point pt= Point(a.cartesian_coords().x(), a.cartesian_coords().y(), a.cartesian_coords().z());
	Point tpt(tr(pt));
	a.set_cartesian_coords(dsrpdb::Point(tpt.x(), tpt.y(), tpt.z()));

	it->second=a;

      }
    
    }
  }
}
