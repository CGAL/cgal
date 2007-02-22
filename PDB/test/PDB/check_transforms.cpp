#include <CGAL/PDB/Transform.h>
#include <CGAL/PDB/Quaternion.h>
#include <iostream>

int main(int, char *[]){
  CGAL::PDB::Quaternion q(CGAL::PDB::Vector(0,1,1), .3);
  //std::cout << q << std::endl;
  CGAL::PDB::Transform t(CGAL::PDB::Vector(0,0,1), q);
  //std::cout << t << std::endl;
  CGAL:PDB::Quaternion qt=t.quaternion();
  //std::cout << qt << std::endl;

  double diff=0;
  for (unsigned int i=0; i< 4; ++i){
    diff += (q[i]-qt[i])*(q[i]-qt[i]);
  }
  CGAL_assertion(diff < .1);

  return EXIT_SUCCESS;
}
