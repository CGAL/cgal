// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//


#ifndef CGAL_TEST_BASIC_H
#define CGAL_TEST_BASIC_H

#include <cassert>
#include <CGAL/functional_base.h>

CGAL_BEGIN_NAMESPACE

template< class Functor >
class Test_functor_arity {
  public:
    void operator()( int arity ) {
      assert( CGAL::Arity_traits< Functor >::Arity::arity == arity );      
    }
};

template<>
class Test_functor_arity< CGAL::Null_functor > {
  public:
    // Test nothing for Null_functors
    void operator()( int ) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_TEST_BASIC
