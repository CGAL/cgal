#include <iostream>
#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Exponent_vector.h>
#include <boost/concept_check.hpp>

int main() {   
  typedef CGAL::Exponent_vector T; 

  boost::function_requires< boost::DefaultConstructibleConcept<T> >();
  boost::function_requires< boost::AssignableConcept<T> >();
  boost::function_requires< boost::EqualityComparableConcept<T> >();
  boost::function_requires< boost::ComparableConcept<T> >();
  boost::function_requires< boost::RandomAccessContainerConcept<T> >();

  CGAL::Exponent_vector ev1,ev2;
  assert(ev1.size()==0);
  assert(ev2.size()==0);
  assert( (ev1==ev2));
  assert(!(ev1!=ev2));
  assert(!(ev1<ev2));
  assert(!(ev1>ev2));
    
  ev1.push_back(1); 
  ev1.push_back(0); 
  ev1.push_back(0);
 
  ev2.push_back(0); 
  ev2.push_back(1); 
  ev2.push_back(0); 
 
  assert(!(ev2 == ev1)) ; 
  assert( (ev2 != ev1)) ; 
  assert(!(ev2 < ev1)) ; 
  assert( (ev2 > ev1)) ;  
    
  assert(CGAL::is_valid(ev1));

  std::vector<int> vec;
  vec.push_back(0); vec.push_back(1); vec.push_back(5); 
    
  ev1 = CGAL::Exponent_vector(vec.begin(),vec.end());
  assert(ev1[0] == 0);  
  assert(ev1[1] == 1);  
  assert(ev1[2] == 5);
   
  // new constructors
  // univariate
  {
    CGAL::Exponent_vector ev(-1);
    assert(!CGAL::is_valid(ev));
  }
  {
    CGAL::Exponent_vector ev(0);
    assert(CGAL::is_valid(ev));
    assert(ev[0] == 0);
  }
  {
    CGAL::Exponent_vector ev(1);
    assert(CGAL::is_valid(ev));
    assert(ev[0] == 1);
  }

  // bivariate
  {
    CGAL::Exponent_vector ev(0,-1);
    assert(!CGAL::is_valid(ev));
  }
  {
    CGAL::Exponent_vector ev(2,1);
    assert(CGAL::is_valid(ev));
    assert(ev[0] == 2);
    assert(ev[1] == 1);
        
  }
    

  // trivariate
  {
    CGAL::Exponent_vector ev(0,0,-1);
    assert(!CGAL::is_valid(ev));
  }
  {
    CGAL::Exponent_vector ev(3,2,1);
    assert(CGAL::is_valid(ev));
    assert(ev[0] == 3);
    assert(ev[1] == 2);
    assert(ev[2] == 1);
  }

  // four-variate
  {
    CGAL::Exponent_vector ev(0,0,0,-1);
    assert(!CGAL::is_valid(ev));
  }
  {
    CGAL::Exponent_vector ev(4,3,2,1);
    assert(CGAL::is_valid(ev));
    assert(ev[0] == 4);
    assert(ev[1] == 3);
    assert(ev[2] == 2);
    assert(ev[3] == 1);
  }
  
  using std::swap; swap(ev1,ev2);

  return 0;
}
