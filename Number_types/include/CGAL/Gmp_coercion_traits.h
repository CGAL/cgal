// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//
// ============================================================================

/*! \file NiX/Gmp/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the Gmp types.
 */

#ifndef CGAL_GMP_COERCION_TRAITS_H
#define CGAL_GMP_COERCION_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/GMP/Gmpz_type.h>
#include <CGAL/GMP/Gmpq_type.h>
#include <CGAL/Coercion_traits.h>

#ifdef CGAL_USE_GMP

CGAL_BEGIN_NAMESPACE

//Gmp internal coercions:    
   
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpz);  
CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Gmpq);  
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(Gmpz,Gmpq);  

// The following definitions reflect the interaction of the Gmp number types
// with the built in types, 
// Gmp BigInt:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,Gmpz);        
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,Gmpz);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,Gmpz);
   
// Gmp BigRat:    
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,Gmpq);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,Gmpq);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,Gmpq);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,Gmpq);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,Gmpq);   

CGAL_END_NAMESPACE


#endif // CGAL_USE_GMP
#endif //CGAL_GMP_COERCION_TRAITS_H 1
//EOF
