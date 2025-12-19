// ============================================================================
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// $URL$
// $Id$
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de>
//
// ============================================================================


#include <cassert>
#include <type_traits>
#include <CGAL/Cache.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/function_objects.h>
#include <CGAL/use.h>
#include <CGAL/assertions.h>

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

    // This is needed to prevent VC7.1 and VC8 to call
    // the explicit templated constructor in Base instead of its copy-ctor.
    Int_t( Int_t const& rhs ) : Base( static_cast<Base const&>(rhs) ) {}
    Int_t& operator=(Int_t const&)=default;
    int  value() const { return this->ptr()->val; }
    void set_value( int i) {
        this->copy_on_write();
        this->ptr()->val = i;
    }
    bool operator==( const Int_t<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            Base::unify(i);
        return equal;
    }
};

void test_typedefs(){
    typedef CGAL::Cache<int,double> Cache;
    CGAL_USE_TYPE(Cache);
    static_assert(::std::is_same< Cache::Input, int >::value );
    static_assert(::std::is_same< Cache::Output,double>::value );
    typedef CGAL::Creator_1<int,double> Creator_double;
    CGAL_USE_TYPE(Creator_double);
    static_assert(::std::is_same<Cache::Creator,Creator_double>::value );
    typedef CGAL::Creator_1<int,int> Creator_int;
    CGAL_USE_TYPE(Creator_int);
    static_assert(::std::is_same<Cache::Canonicalizer,Creator_int>::value );
    static_assert(::std::is_same<Cache::Compare,std::less<int> >::value );
    static_assert(::std::is_same<Cache::Self,CGAL::Cache<int,double> >::value );
}
int main(){
    {
        test_typedefs();
        {
            typedef CGAL::Cache<int,double> Cache;
            double d;
            Cache cache;
            assert(cache.is_empty());
            assert(cache.size()==0);
            d=cache(3);
            assert(d==double(3));
            assert(cache.size()==1);
            d=cache(4);
            assert(d==double(4));
            assert(cache.size()==2);
            d=cache(3);
            assert(d==double(3));
            assert(cache.size()==2);
            d=cache(2);
            assert(d==double(2));
            assert(cache.size()==3);

            typedef Cache::Iterator Iterator;
            CGAL_USE_TYPE(Cache::Const_iterator);
            CGAL_USE_TYPE(Cache::Reverse_iterator);
            CGAL_USE_TYPE(Cache::Const_reverse_iterator);
            CGAL_USE_TYPE(Cache::Size_type);

            Iterator it;
            d=0;
            for(it=cache.begin();it!=cache.end();it++){
                assert(d<(*it).second);
                d=(*it).second;
            }
            cache.clear();
            assert(cache.size()==0);
        }
        {
            typedef Int_t< CGAL::Handle_policy_no_union > Int;
            typedef CGAL::Cache<int,Int> Cache;
            Int hi,hi2;
            Cache cache;
            assert(cache.is_empty());
            assert(cache.size()==0);
            hi=cache(3);
            assert(hi==Int(3));
            assert(cache.size()==1);
            hi2=cache(4);
            assert(hi2==Int(4));
            assert(cache.size()==2);
            hi2=cache(3);
            assert(hi.id()==hi2.id());
        }
    }

return EXIT_SUCCESS;
}
