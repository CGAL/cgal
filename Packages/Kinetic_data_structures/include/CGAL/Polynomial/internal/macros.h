#ifndef CGAL_POLYNOMIAL_INTERNAL_MACROS_H
#define CGAL_POLYNOMIAL_INTERNAL_MACROS_H

#include <CGAL/Polynomial/internal/config.h>

#ifdef POLYNOMIAL_USE_CGAL
/*
  When CGAL is present
*/
#include <CGAL/basic.h>


#define CGAL_POLYNOMIAL_BEGIN_NAMESPACE CGAL_BEGIN_NAMESPACE \
namespace POLYNOMIAL {

#define CGAL_POLYNOMIAL_END_NAMESPACE } \
CGAL_END_NAMESPACE


#define POLYNOMIAL_NS CGAL::POLYNOMIAL
#define Polynomial_assertion(x) CGAL_assertion(x)
#define Polynomial_assertion_code(x) CGAL_assertion_code(x)
#define Polynomial_precondition(x) CGAL_precondition(x)
#define Polynomial_precondition_code(x) CGAL_precondition_code(x)
#define Polynomial_postcondition(x) CGAL_postcondition(x)
#define Polynomial_expensive_precondition(x) CGAL_expensive_precondition(x)
#define Polynomial_expensive_assertion(x) CGAL_expensive_assertion(x)
#define Polynomial_expensive_postcondition(x) CGAL_expensive_postcondition(x)
#define Polynomial_exactness_assertion(x) CGAL_exactness_assertion(x)
#define Polynomial_exactness_postcondition(x) CGAL_exactness_postcondition(x)
#define Polynomial_exactness_precondition(x) CGAL_exactness_precondition(x)

#else
/*
  When no CGAL is present
*/


#define CGAL_POLYNOMIAL_BEGIN_NAMESPACE \
namespace Polynomial {

#define CGAL_POLYNOMIAL_END_NAMESPACE }

#define POLYNOMIAL_NS Polynomial


#include <assert.h>

#define Polynomial_assertion(x) assert(x)
// This does not work
#define Polynomial_assertion_code(x) x 
#define Polynomial_precondition(x) assert(x)
#define Polynomial_postcondition(x) assert(x)
#define Polynomial_expensive_precondition(x) 
#define Polynomial_expensive_assertion(x)
#define Polynomial_expensive_postcondition(x)
#define Polynomial_exactness_postcondition(x) 
#define Polynomial_exactness_precondition(x) 





#endif


#define CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE CGAL_POLYNOMIAL_BEGIN_NAMESPACE \
namespace internal {

#define CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE } \
CGAL_POLYNOMIAL_END_NAMESPACE


#endif

