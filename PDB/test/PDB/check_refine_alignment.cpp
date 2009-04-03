#include <CGAL/PDB/internal/align_points.h>
#include <vector>
#include <iostream>
#include <iterator>

int main(int, char *[]){
  typedef CGAL::PDB::Point P;
  std::vector<P> hairpin;
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
    CGAL::PDB::Transform tr(1,0,0,0,
                            0,1,0,0,
                            0,0,1,1);

    std::vector<P> trp;
    for (unsigned int i=0; i< hairpin.size(); ++i){
      trp.push_back(tr(hairpin[i]));
    }

    std::vector<int> match;
    CGAL::PDB::Transform tro
      = CGAL::PDB::internal::refine_alignment(hairpin,
                                              trp,
                                              .25,
                                              std::back_inserter(match));
    //std::cout << tr << std::endl;
    //std::cout << tro << std::endl;
    std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
  }


  {
    std::vector<P> hairpin_2(hairpin);
    hairpin_2.insert(hairpin_2.begin(), P(-1,0,0));
    hairpin_2.push_back(P(-1,1,0));
    std::cout << hairpin_2[0] << std::endl;

    std::vector<int> match(hairpin.size(), -1);
    CGAL::PDB::internal::greedy_matching(hairpin,
                                         hairpin_2,
                                         match);
    std::copy(match.begin(), match.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
