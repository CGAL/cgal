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

/*! \file NiX/CORE/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the CORE types.
 */

#ifndef CGAL_CORE_COERCION_TRAITS_H
#define CGAL_CORE_COERCION_TRAITS_H 1

#include <CGAL/number_type_basic.h>

#ifdef CGAL_USE_CORE
#include <CORE/BigInt.h>
#include <CORE/BigRat.h>
#include <CORE/BigFloat.h>
#include <CORE/Expr.h>

//#include <NiX/CORE/Coercion_traits.h>

CGAL_BEGIN_NAMESPACE

//CORE internal coercions:    
   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigInt,::CORE::BigFloat);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigInt,::CORE::BigRat);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigInt,::CORE::Expr);

    // CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigFloat,::CORE::BigRat);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigFloat,::CORE::Expr);

    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigRat,::CORE::Expr);

// The following definitions reflect the interaction of the CORE number types
// with the built in types, 
// CORE BigInt:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,::CORE::BigInt);        
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,::CORE::BigInt);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,::CORE::BigInt);
   

// CORE BigRat:    
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::CORE::BigRat);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::CORE::BigRat);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::CORE::BigRat);
    
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::CORE::BigRat);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::CORE::BigRat);
   

// CORE BigFloat:      :
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::CORE::BigFloat);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::CORE::BigFloat);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::CORE::BigFloat); 
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::CORE::BigFloat);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::CORE::BigFloat);

// CORE Expr:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::CORE::Expr);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::CORE::Expr);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::CORE::Expr);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::CORE::Expr);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::CORE::Expr);

// not provieded by CORE 
// Note that this is not symmetric to LEDA
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigInt);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigRat);  
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::BigRat);  
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigFloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::BigFloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::Expr); 
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::Expr);

CGAL_END_NAMESPACE

#endif // CGAL_USE_CORE
#endif //CGAL_CORE_COERCION_TRAITS_H 1
//EOF
