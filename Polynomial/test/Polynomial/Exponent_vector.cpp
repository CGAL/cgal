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
    
    ev1 = CGAL::Exponent_vector(3);
    ev2 = CGAL::Exponent_vector(3);
    
    ev1[0]=1;
    ev2[1]=1;
 
    assert(!(ev2 == ev1)) ; 
    assert( (ev2 != ev1)) ; 
    assert(!(ev2 < ev1)) ; 
    assert( (ev2 > ev1)) ;  
    
    
    assert(CGAL::is_valid(ev1));
    
    return 0;
}
