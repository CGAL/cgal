#include <iostream>
#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Exponent_vector.h>

int main() {   
  
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
   
    
    return 0;
}
