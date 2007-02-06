#include <CGAL/PDB/Residue.h>
#include <CGAL/PDB/iterator.h>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/align.h>
#include <CGAL/PDB/align_generic.h>
#include <cassert>


CGAL_PDB_BEGIN_NAMESPACE

 

Transform transform_taking_first_to_second(const std::vector<Point> &a,
					   const std::vector<Point> &b) {
  Transform tr
    = transform_taking_first_to_second(a.begin(),
					       a.end(),
					       b.begin(),
					       b.end());
  return tr;
}


void align_second_protein_to_first(const Protein &base, 
				   Protein &o) {
  assert(base.number_of_residues()==o.number_of_residues());
   
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
     

  for (Protein::Residues_iterator it= o.residues_begin(); it != o.residues_end(); ++it){
    Residue &aa= *it;
    //dsr::vector<Residue::Atom_label> als= aa->atoms();
       
    for (Residue::Atoms_iterator it= aa.atoms_begin(); it != aa.atoms_end(); ++it){
      Atom a= it->second;
      Point pt= Point(a.cartesian_coords().x(), a.cartesian_coords().y(), a.cartesian_coords().z());
      Point tpt(tr(pt));
      a.set_cartesian_coords(Point(tpt.x(), tpt.y(), tpt.z()));

      it->second=a;

    }
    
  }
}

Transform refine_alignment_of_second_protein_to_first(const Protein &base,
						      Protein &o) {
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
  for (Protein::Residues_iterator rit= o.residues_begin(); rit != o.residues_end(); ++rit){
    for (Residue::Atoms_iterator ait= rit->atoms_begin(); ait != rit->atoms_end(); ++ait){
      ait->second.set_cartesian_coords(t(ait->second.cartesian_coords()));
    }
  }
  return t;
}


CGAL_PDB_END_NAMESPACE
