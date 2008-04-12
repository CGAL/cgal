#include <CGAL/PDB/small_map.h>
#include <CGAL/Tools/Label.h>
#include <CGAL/Random.h>
#include <CGAL/PDB/Monomer.h>
#include <CGAL/Gmpq.h>
#include <set>
#include <cassert>

typedef CGAL::Gmpq Data;
typedef CGAL::Label<int> Key;
CGAL_SMALL_MAP_VALUE_TYPE(SM_vt, Key, Data, data);
typedef CGAL::PDB::small_map<SM_vt> SM;


int main(int, char *[]) {
  //typedef CGAL::Label<int> Data;
  SM sm;
  CGAL::Random rand;
  for (unsigned int i=0; i< 100; ++i) {
    unsigned int cur= rand.get_int(0, 500);
    if (sm.find(Key(cur)) == sm.end()) {
      sm.insert(SM::value_type(Key(cur), Data(cur,1)));
    }
  }
  for (unsigned int i=0; i< 100; ++i) {
    unsigned int cur= rand.get_int(0, 500);
    sm[Key(cur)]= Data(static_cast<int>(cur));
  }
  {
    std::set<unsigned int> vals;
    for (SM::const_iterator it= sm.begin(); it != sm.end(); ++it) {
      assert(Data(it->key().index(), 1) == it->data());
      assert(vals.find(it->key().index()) == vals.end());
      vals.insert(it->key().index());
    }
    std::cout << "size is " << sm.size() << std::endl;
  }

  for (unsigned int i=0; i< 100; ++i) {
    int cur= 700-i;
    sm.insert(SM::value_type(Key(cur), Data(cur)));
  }
 {
    std::set<unsigned int> vals;
    for (SM::const_iterator it= sm.begin(); it != sm.end(); ++it) {
      assert(Data(it->key().index(),1) == it->data());
      assert(vals.find(it->key().index()) == vals.end());
      vals.insert(it->key().index());
    }
    std::cout << "size is " << sm.size() << std::endl;
  }
 
 {
   typedef CGAL::PDB::small_map<CGAL::PDB::small_map_value_type<int, CGAL::PDB::Monomer> > MM;
   MM mm;
   CGAL::PDB::Monomer m(CGAL::PDB::Monomer::MET);
   mm.insert(MM::value_type(0, m));
   std::cout << "Outside size is " << sizeof(CGAL::PDB::Monomer) << std::endl;
   CGAL::PDB::Monomer mt= mm[0];
   CGAL_postcondition(mm[0].type() == CGAL::PDB::Monomer::MET);
 }
  return EXIT_SUCCESS;
}
