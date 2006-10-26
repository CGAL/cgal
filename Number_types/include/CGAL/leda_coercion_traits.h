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

/*! \file NiX/LEDA/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the LEDA number types.
 */

#ifndef CGAL_LEDA_COERCION_TRAITS_H
#define CGAL_LEDA_COERCION_TRAITS_H

#include <CGAL/Coercion_traits.h>
#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/integer.h>
#include <LEDA/bigfloat.h>
#include <LEDA/rational.h>
#include <LEDA/real.h>
#else
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/bigfloat.h>
#include <LEDA/numbers/rational.h>
#include <LEDA/numbers/real.h>
#endif


CGAL_BEGIN_NAMESPACE


//LEDA internal coercions:
   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::bigfloat); 
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::rational);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::real);

    //NiX_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat,::leda::rational);
//    NiX_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat,::leda::real);

    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::rational,::leda::real);

// The following definitions reflect the interaction of the LEDA number types
// with the built in types, 
// leda integer:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,::leda::integer);        
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,::leda::integer);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,::leda::integer);

// leda rational:    
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::rational);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::rational);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::rational);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::leda::rational);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::rational);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::rational);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::rational);

// leda bigfloat:      :
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::bigfloat);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::bigfloat);
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::bigfloat);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::bigfloat);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::bigfloat);

// leda real:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::real);   
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::real);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::real);  
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::real);

//not provided by LEDA
//Note that this is not symmetric to CORE 
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,::leda::integer);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::leda::bigfloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::bigfloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::real);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::leda::real); 
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::real);

CGAL_END_NAMESPACE
#endif // CGAL_USE_LEDA
#endif //CGAL_LEDA_COERCION_TRAITS_H
//EOF
