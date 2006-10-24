// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : Support
// File          : test/Cache.C
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de>
//
// ============================================================================


#include <CGAL/basic.h>
#include <CGAL/Testsuite/assert.h>
#include <CGAL/Cache.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/function_objects.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

struct Int_rep {
    int val;
    Int_rep( int i = 0) : val(i) {}
    Int_rep( int i, int j) : val(i+j) {}
    Int_rep( int i, int j, int k) : val(i+j+k) {}
};

template < class Unify>
struct Int_t : public CGAL::Handle_with_policy< Int_rep, Unify > {
    typedef CGAL::Handle_with_policy< Int_rep, Unify > Base;
    Int_t( int i = 0) : Base( i) {}
    Int_t( int i, int j) : Base( i, j) {}     // test template constructors
    Int_t( int i, int j, int k) : Base( Base::USE_WITH_INITIALIZE_WITH) {
        // test initialize_with
        this->initialize_with( i, j + k);
    }
    int  value() const { return this->ptr()->val; }
    void set_value( int i) {
        this->copy_on_write();
        this->ptr()->val = i;
    }
    bool operator==( const Int_t<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            unify(i);
        return equal;
    }
};

void test_typedefs(){
    typedef CGAL::Cache<int,double> Cache;
    BOOST_STATIC_ASSERT(( ::boost::is_same< Cache::Input, int >::value ));
    BOOST_STATIC_ASSERT(( ::boost::is_same< Cache::Output,double>::value ));
    typedef CGAL::Creator_1<int,double> Creator_double; 
    BOOST_STATIC_ASSERT(( ::boost::is_same<Cache::Creator,Creator_double>::value ));
    typedef CGAL::Creator_1<int,int> Creator_int; 
    BOOST_STATIC_ASSERT(( ::boost::is_same<Cache::Canonicalizer,Creator_int>::value ));
    BOOST_STATIC_ASSERT(( ::boost::is_same<Cache::Compare,std::less<int> >::value ));
    BOOST_STATIC_ASSERT(( ::boost::is_same<Cache::Self,CGAL::Cache<int,double> >::value ));
}
int main(){
    {
        test_typedefs(); 
        {
            typedef CGAL::Cache<int,double> Cache;
            double d;
            Cache cache;
            CGAL_test_assert(cache.is_empty());
            CGAL_test_assert(cache.size()==0);
            d=cache(3);
            CGAL_test_assert(d==double(3));
            CGAL_test_assert(cache.size()==1);
            d=cache(4);
            CGAL_test_assert(d==double(4));
            CGAL_test_assert(cache.size()==2);
            d=cache(3);
            CGAL_test_assert(d==double(3));
            CGAL_test_assert(cache.size()==2);
            d=cache(2);
            CGAL_test_assert(d==double(2));
            CGAL_test_assert(cache.size()==3);

            typedef Cache::Iterator Iterator;
            typedef Cache::Const_iterator Const_iterator;
            typedef Cache::Reverse_iterator Reverse_iterator;
            typedef Cache::Const_reverse_iterator Const_reverse_iterator;
            typedef Cache::Size_type Size_type;
        
            Iterator it;
            d=0;
            for(it=cache.begin();it!=cache.end();it++){
                CGAL_test_assert(d<(*it).second);
                d=(*it).second;
            }
            cache.clear();
            CGAL_test_assert(cache.size()==0);
        }
        {
            typedef Int_t< CGAL::Handle_policy_no_union > Int;
            typedef CGAL::Cache<int,Int> Cache;
            Int hi,hi2;
            Cache cache;
            CGAL_test_assert(cache.is_empty());
            CGAL_test_assert(cache.size()==0);
            hi=cache(3);
            CGAL_test_assert(hi==Int(3));
            CGAL_test_assert(cache.size()==1);
            hi2=cache(4);
            CGAL_test_assert(hi2==Int(4));
            CGAL_test_assert(cache.size()==2);
            hi2=cache(3);
            CGAL_test_assert(hi.id()==hi2.id());       
        }
    }

return EXIT_SUCCESS;
}
