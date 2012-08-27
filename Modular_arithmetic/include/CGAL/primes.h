// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_PRIMES_H
#define CGAL_PRIMES_H

#include <CGAL/basic.h>

namespace CGAL { 
namespace internal {
CGAL_EXPORT extern int primes[2000];

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
