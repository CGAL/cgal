#include <CGAL/Kinetic/Active_objects_vector.h>
//#include <CGAL/Kinetic/Active_objects_set.h>

template <class AOT>
void test_erase(AOT &aot, CGAL::Tag_true) {
  std::set<typename AOT::Key> erased;
  std::set<typename AOT::Key> remaining;
  int i=0;
  aot.set_is_editing(true);
  for (typename AOT::Key_iterator it= aot.keys_begin(); it != aot.keys_end(); ++it){
    if (i%2 ==0) {
      erased.insert(*it);
      aot.erase(*it);
    } else {
      remaining.insert(*it);
    }
    ++i;
  }
  aot.set_is_editing(false);
  for (typename AOT::Key_iterator it= aot.keys_begin(); it != aot.keys_end(); ++it){
    assert(erased.find(*it) == erased.end());
    assert(remaining.find(*it) != remaining.end());
    ++i;
  }
  assert(aot.size() == remaining.size());
}

template <class AOT>
void test_erase(AOT &, CGAL::Tag_false) {
  
}

template <class AOT, class Erase_tag> 
void test(AOT &aot, Erase_tag et) {
  std::vector<typename AOT::Key> keys;
  for (unsigned int i=0; i< 100; ++i){
    keys.push_back(aot.insert(i));
  }
  assert(static_cast<unsigned int> (std::distance(aot.keys_begin(), 
							  aot.keys_end()))
		 == keys.size());

  for (int i=0; i< 100; ++i){
    assert(aot[keys[i]]==i);
  }
  assert(aot.size()==100);
  
  test_erase(aot, et);
}

int main(int, char *[]) {
  typedef CGAL::Kinetic::Active_objects_vector<int> AOV;
  // typedef CGAL::Kinetic::Active_objects_set<int> AOS;
  
  AOV aov;
  test(aov, CGAL::Tag_true());
  
  /*AOS aos;
    test(aos, CGAL::Tag_true());*/

  return EXIT_SUCCESS;
}
