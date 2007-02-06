#include <CGAL/PDB/align_generic.h>
#include <CGAL/PDB/align.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/geometry.h>
#include <CGAL/PDB/transforms.h>

int main(int, char *[]){
  typedef CGAL_PDB_NS::Point P;
  std::vector<P> hairpin;
  typedef std::vector<P>::const_iterator H;
  hairpin.push_back(P(0, 0,0));
  hairpin.push_back(P(1, 0,0));
  hairpin.push_back(P(2, 0,0));
  hairpin.push_back(P(3, 0,0));
  hairpin.push_back(P(4, 0,0));
  hairpin.push_back(P(4,.5,0));
  hairpin.push_back(P(4, 1,0));
  hairpin.push_back(P(3, 1,0));
  hairpin.push_back(P(2, 1,0));
  hairpin.push_back(P(1, 1,0));
  hairpin.push_back(P(0, 1,0));
  
 
  {
    CGAL_PDB_NS::Transform tr;
    tr.set_translation(P(0,0,1));

    std::vector<P> trp;
    for (unsigned int i=0; i< hairpin.size(); ++i){
      trp.push_back(tr(hairpin[i]));
    }

    std::vector<H> match;
    CGAL_PDB_NS::Transform tro= CGAL_PDB_NS::refine_alignment(hairpin.begin(), hairpin.end(),
						    trp.begin(), trp.end(),
						    .01,
						    std::back_inserter(match));
    /*std::cout << tr << std::endl;
    std::cout << tro << std::endl;
    std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;*/
    assert(match.size()==hairpin.size());
    for (unsigned int i=0; i< match.size(); ++i){
      //int d= match[i]-trp.begin();
      assert(match[i]- trp.begin() == static_cast<int>(i));
    }
    std::cout << tr << std::endl << tro << std::endl;
  }

  {
    CGAL_PDB_NS::Transform tr;
    tr.set_translation(P(0,.5,1));

    std::vector<P> trp;
    for (unsigned int i=0; i< hairpin.size(); ++i){
      trp.push_back(tr(hairpin[i]));
    }

    std::vector<H> match;
    CGAL_PDB_NS::Transform tro= CGAL_PDB_NS::refine_alignment(hairpin.begin(), hairpin.end(),
						    trp.begin(), trp.end(),
						    .01,
						    std::back_inserter(match));
    /*std::cout << tr << std::endl;
    std::cout << tro << std::endl;
    std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;*/
    assert(match.size()==hairpin.size());
    for (unsigned int i=0; i< match.size(); ++i){
      assert(match[i] - trp.begin() == static_cast<int>(i));
    }
    std::cout << tr << std::endl << tro << std::endl;
  }


  {
    std::vector<P> hairpin_2;
    hairpin_2.push_back(P(-1,0,0));
    hairpin_2.push_back(P(0,0,0));
    for (unsigned int i=0; i< hairpin.size(); ++i){
      hairpin_2.push_back(P(hairpin[i].x()+1, hairpin[i].y(), hairpin[i].z()));
    }
    hairpin_2.push_back(P(0,1,0));
    hairpin_2.push_back(P(-1,1,0));

    std::vector<H> match;
    CGAL_PDB_NS::refine_alignment(hairpin.begin(), hairpin.end(),
			     hairpin_2.begin(), hairpin_2.end(),
			     .01,
			     std::back_inserter(match));
    //std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
    for (unsigned int i=0; i< match.size(); ++i){
      std::cout << match[i]- hairpin_2.begin() << " ";
    }
    std::cout << std::endl;
    assert(match.size() == hairpin.size());
    for (unsigned int i=0; i< 5; ++i){
      //assert(match[i]-hairpin_2.begin()==static_cast<int>(i));
    }
    for (unsigned int i=0; i< 6; ++i){
      //assert(match[i+5]-hairpin_2.begin()== static_cast<int>(i+9));
    }
  }

  {
    std::ifstream sheetstream("check_bonds.pdb");
    CGAL_PDB_NS::Protein sheet(sheetstream);
    std::ifstream modelstream("check_refine_alignment_0.pdb");
    CGAL_PDB_NS::Protein model(modelstream);

    std::vector<CGAL_PDB_NS::Point> sheet_points, model_points;
    CGAL_PDB_NS::backbone_coordinates(sheet.atoms_begin(), sheet.atoms_end(),
				 std::back_inserter(sheet_points));
    CGAL_PDB_NS::backbone_coordinates(model.atoms_begin(), model.atoms_end(),
				 std::back_inserter(model_points));

    std::vector<H> match;
    CGAL_PDB_NS::refine_alignment(sheet_points.begin(), sheet_points.end(),
			     model_points.begin(), model_points.end(),
			     .01,
			     std::back_inserter(match));
    /*std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
      std::cout << std::endl;*/
    // 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 
    for (unsigned int i=0; i < 15; ++i){
      assert(match[i]-model_points.begin() == 84+static_cast<int>(i));
    }
    for (unsigned int i=0; i< 14; ++i){
      //int d= match[15+i]-model_points.begin();
      assert(match[15+i]-model_points.begin() == 138+ static_cast<int>(i));
    }

    {
      CGAL_PDB_NS::Transform tr;
      tr.set_translation(P(1,1,1));
      for (unsigned int i=0; i< sheet_points.size(); ++i){
	sheet_points[i] = tr(sheet_points[i]);
      }
      std::vector<H> match2;
      CGAL_PDB_NS::refine_alignment(sheet_points.begin(), sheet_points.end(),
			       model_points.begin(), model_points.end(),
			       .01,
			       std::back_inserter(match2));
      for (unsigned int i=0; i < match.size(); ++i){
	assert(match2[i] == match[i]);
      }
    }
    
  }

  {

    std::ifstream t1("check_refine_alignment_1.pdb");
    CGAL_PDB_NS::Protein pt1(t1);
    std::ifstream t2("check_refine_alignment_2.pdb");
    CGAL_PDB_NS::Protein pt2(t2);
    

    CGAL_PDB_NS::Transform tr;
    tr.set_translation(P(1,1,1));
    CGAL_PDB_NS::Protein cpt1(pt1);
    CGAL_PDB_NS::transform_protein(tr, cpt1);
    {
      CGAL_PDB_NS::Protein lpt2(pt2);
      CGAL_PDB_NS::Transform otr0= CGAL_PDB_NS::refine_alignment_of_second_protein_to_first(pt1, lpt2);
      std::cout << otr0 << std::endl;
    }
    {
      CGAL_PDB_NS::Protein lpt1(pt1);
      CGAL_PDB_NS::Transform otr0= CGAL_PDB_NS::refine_alignment_of_second_protein_to_first(pt2, lpt1);
      std::cout << otr0 << std::endl;
    }
    {
      CGAL_PDB_NS::Protein lcpt1(cpt1);
      CGAL_PDB_NS::Transform otr1 = CGAL_PDB_NS::refine_alignment_of_second_protein_to_first(pt1, lcpt1);
      std::cout << otr1 << std::endl;
    }
    {
      CGAL_PDB_NS::Protein lcpt1(cpt1);
      CGAL_PDB_NS::Transform otr2 = CGAL_PDB_NS::refine_alignment_of_second_protein_to_first(pt2, cpt1);
      std::cout << otr2 << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
