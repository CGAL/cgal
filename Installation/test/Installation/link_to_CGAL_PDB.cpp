// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/PDB/PDB.h>

int main()
{
  std::ifstream in("");  
  
  volatile CGAL::PDB::PDB pdb(in,false);
  
  return 0;
}
