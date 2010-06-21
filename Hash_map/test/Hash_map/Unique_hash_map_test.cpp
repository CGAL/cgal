#include <CGAL/Cartesian.h>
#include <list>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/test_macros.h>

using namespace std;
typedef list<int>::iterator Iterator;

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
    std::cerr << "done" << std::endl;
    CGAL_TEST_END;
}


