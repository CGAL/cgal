#include <list>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/test_macros.h>

using namespace std;
typedef list<int>::iterator Iterator;

struct Integer_hash_function
{
  int operator()(int i) const { return i; }
};

int main() {
    CGAL_TEST_START;
    list<int> L;
    L.push_back(1);
    L.push_back(2);

    Iterator it1 = L.begin();
    CGAL::Unique_hash_map<Iterator,int> H1;
    CGAL::Unique_hash_map<Iterator,int> H2(-1);
    CGAL_TEST( H1.default_value() ==  0);
    CGAL_TEST( H2.default_value() == -1);
    H1[it1] = 2;
    CGAL_TEST(H1[it1]==2);
    CGAL_TEST(H2[it1]==-1);
    H1.clear();
    H2.clear(-2);
    H2[it1] = 2;
    CGAL_TEST(H1[it1]==0);
    CGAL_TEST(H2[it1]==2);
    Iterator it2 = L.end();
    --it2;
    const CGAL::Unique_hash_map<Iterator,int>* pH = &H2;
    CGAL_TEST((*pH)[it2]==-2);

    H1.clear(-1);
    L.push_back(3);
    L.push_back(4);
    L.push_back(5);
    CGAL_TEST( L.size() == 5);
    H1.insert( L.begin(), L.end(), 1);
    for ( Iterator i = L.begin(); i != L.end(); ++i) {
        CGAL_TEST( H1[i] == *i);
    }
    CGAL::Unique_hash_map<Iterator,int> H3( L.begin(), L.end(), 1);
    for ( Iterator j = L.begin(); j != L.end(); ++j) {
        CGAL_TEST( H3[j] == *j);
    }
    CGAL::Handle_hash_function hash;
    CGAL::Unique_hash_map<Iterator,int,CGAL::Handle_hash_function>
        H4( L.begin(), L.end(), 1, -1, 512, hash);
    for ( Iterator k = L.begin(); k != L.end(); ++k) {
        CGAL_TEST( H4[k] == *k);
    }
    hash = H4.hash_function();

    // test the overload of boost::associative_property_map for CGAL::Unique_hash_map
    typedef CGAL::Unique_hash_map<Iterator, int> Iterator_hmap;
    typedef boost::associative_property_map<Iterator_hmap> Iterator_pmap;
    Iterator_pmap H4_pmap = boost::make_assoc_property_map(H4);
    for(Iterator k=L.begin(); k!=L.end(); ++k){
      CGAL_TEST(get(H4_pmap, k) == *k);
    }
    L.push_front(0);
    put(H4_pmap, L.begin(), 0);
    CGAL_TEST(get(H4_pmap, L.begin()) == 0);

    typedef CGAL::Unique_hash_map<int, int, Integer_hash_function>  Int_hmap;
    typedef boost::associative_property_map<Int_hmap>               Int_pmap;
    Int_hmap H5(-1);
    Int_pmap H5_pmap(H5);
    put(H5_pmap, -1, 1);
    CGAL_TEST(get(H5_pmap, -1) == 1);
    CGAL_TEST(H5_pmap[0] == -1);

    std::cerr << "done" << std::endl;
    CGAL_TEST_END;
}


