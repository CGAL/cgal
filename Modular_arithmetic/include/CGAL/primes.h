
#ifndef CGAL_PRIMES_H
#define CGAL_PRIMES_H

#include <CGAL/basic.h>

namespace CGAL { 
namespace CGALi {
extern int primes[2000];

static inline 
int get_next_lower_prime(int current_prime){
    bool is_prime = false;
    
    int i;
    CGAL_precondition_msg(current_prime != 2 ," primes definitely exhausted ");

    if((current_prime <= 7) && (current_prime > 2)){
        if(current_prime <= 5){
            if(current_prime == 3)
                return 2;
            return 3;
        }
        return 5;
    }                
    for(i=current_prime-2;(i>1 && !is_prime);i=i-2){
        int r = 1;
        for(int j=3; (j <= i/2 && (r != 0)); j++){
            r = i % j;
//                std::cout<<"i " <<i<<std::endl;
//                std::cout<<"j " <<j<<std::endl;
//                std::cout<<"i%j " <<i%j<<std::endl;
            if(j==i/2 && r != 0)
                is_prime = true;
        }
    }
//    CGAL_precondition_msg(is_prime," primes definitely exhausted ");
    return i+2;
}


}
}

#endif // CGAL_PRIMES_H
