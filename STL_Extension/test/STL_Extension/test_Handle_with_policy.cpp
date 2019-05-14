// ============================================================================
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// $URL$
// $Id$
// Author(s)     : Michael Seel <seel@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Lutz Kettner <kettner@mpi-inf.mpg.de>
//
// ============================================================================

#define CGAL_HANDLE_WITH_POLICY_INTERNAL_TEST

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <cassert>

#include <CGAL/Handle_with_policy.h>
#include <cstdlib>

struct Int_rep {
    int val;
    Int_rep( int i = 0) : val(i) {}
    Int_rep( int i, int j) : val(i+j) {}
    Int_rep( int i, int j, int k) : val(i+j+k) {}
};

template < class Unify>
struct Int_t : public ::CGAL::Handle_with_policy< Int_rep, Unify > {
    typedef ::CGAL::Handle_with_policy< Int_rep, Unify > Base;
    Int_t( int i = 0) : Base( i) {}
    Int_t( int i, int j) : Base( i, j) {}     // test template constructors
    Int_t( int i, int j, int k) : Base( Base::USE_WITH_INITIALIZE_WITH) {
        // test initialize_with
        this->initialize_with( i, j + k);
    }

    // This is needed to prevent VC7.1 and VC8 to call
    // the explicit templated constructor in Base instead of its copy-ctor.
    Int_t( Int_t const& rhs ) : Base( static_cast<Base const&>(rhs) ) {}

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

void test_handle() {
    { // test template constructor and initialize_with
        typedef Int_t< ::CGAL::Handle_policy_no_union> Int;
        Int i(5,6);
        assert( i == Int(11));
        Int j(5,6,7);
        assert( j == Int(18));
    }
    {
        typedef Int_t< ::CGAL::Handle_policy_in_place> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
    }
    {
        typedef Int_t< ::CGAL::Handle_policy_no_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( (std::ptrdiff_t)(ID_Number(i)) == i.id() );
    }
    {
        typedef Int_t< ::CGAL::Handle_policy_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
    }
    {
        typedef Int_t< ::CGAL::Handle_policy_union_and_reset> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
    }
    {
        typedef Int_t< ::CGAL::Handle_policy_union> Int;
        Int i(5);
        Int j(5);
        Int k(5);
        Int l(5);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // pump up the union_size counter for j, k and l
        Int j1(j);
        assert( j1.test_identical_ptr( j));
        Int k1(k);
        Int k2(k);
        Int k3(k);
        Int l1(l);
        Int l2(l);
        Int l3(l);
        Int l4(l);
        Int l5(l);
        Int l6(l);
        Int l7(l);
        Int l8(l);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 2);
        assert( j.union_size() == 3);
        assert( k.union_size() == 5);
        assert( l.union_size() == 10);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link i to j
        assert( i == j);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 4);
        assert( j.union_size() == 4);
        assert( k.union_size() == 5);
        assert( l.union_size() == 10);
        assert( i.test_identical_ptr( j));
        assert( i.test_identical_ptr( j1));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link j to k
        assert( j == k);
        assert( i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 9);
        assert( k.union_size() == 9);
        assert( l.union_size() == 10);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link k to l
        assert( k == l);
        assert( i.is_forwarding());
        assert( j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 8);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
        // find j, links it to k and l
        assert( j.value() == 5);
        assert( i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 19);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
        // find i, links it to j, k and l
        assert( i.value() == 5);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 19);
        assert( j.union_size() == 19);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( i.test_identical_ptr( j));
        assert( i.test_identical_ptr( k));
        assert( i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
    }
}

// fully generic example to show, how allocator and policy can be transported
// from the templated handle to the templated rep. For hard-coded choices
// one could derive directly from one of the following base classes:
//  -- struct Int_vrep : public ::CGAL::Reference_counted_hierarchy_with_union<>
//  -- struct Int_vrep : public ::CGAL::Reference_counted_hierarchy<>

template <class Policy, class Alloc>
struct Int_vrep : public Policy::template Hierarchy_base< Alloc>::Type {
    int val;
    virtual CGAL::Reference_counted_hierarchy<Alloc>* clone() {
        return new Int_vrep( *this);
    }
    virtual int  get_val() const { return val; }
    virtual void set_val( int i) { val = i; }
    Int_vrep( int i = 0) : val(i) {}
};

template <class Policy, class Alloc>
struct Int_vrep2 : public Int_vrep<Policy,Alloc> {
    int val2;
    virtual ::CGAL::Reference_counted_hierarchy<Alloc>* clone() {
        return new Int_vrep2( *this);
    }
    virtual int get_val() const { return this->val + val2; }
    virtual void set_val( int i) { this->val = i - val2; }
    Int_vrep2( int i, int j) : Int_vrep<Policy,Alloc>(i), val2(j) {}
};

template <class Policy, class Alloc>
struct Int_vrep3 : public Int_vrep2<Policy,Alloc> {
    int val3;
    virtual ::CGAL::Reference_counted_hierarchy<Alloc>* clone() {
        return new Int_vrep3( *this);
    }
    virtual int get_val() const { return this->val + this->val2 + val3; }
    virtual void set_val( int i) { this->val = i - this->val2 - val3; }
    Int_vrep3( int i, int j, int k) : Int_vrep2<Policy,Alloc>(i,j), val3(k) {}
};

template < class Unify, class Alloc = CGAL_ALLOCATOR(char) >
struct Int_vt : public ::CGAL::Handle_with_policy< Int_vrep<Unify,Alloc>, Unify > {
    typedef ::CGAL::Handle_with_policy< Int_vrep<Unify,Alloc>, Unify > Base;
    Int_vt( int i = 0) : Base( new Int_vrep<Unify,Alloc>(i)) {}
    Int_vt( int i, int j) : Base( new Int_vrep2<Unify,Alloc>(i,j)) {}
    Int_vt( int i, int j, int k) : Base( new Int_vrep3<Unify,Alloc>(i,j,k)) {}

    // This is needed to prevent VC7.1 and VC8 to call
    // the explicit templated constructor in Base instead of its copy-ctor.
    Int_vt( Int_vt const& rhs ) : Base( static_cast<Base const&>(rhs) ) {}

    int  value() const { return this->ptr()->get_val(); }
    void set_value( int i) {
        this->copy_on_write();
        this->ptr()->set_val(i);
    }
    bool operator==( const Int_vt<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            Base::unify(i);
        return equal;
    }
};



void test_handle_with_class_hierarchy() {
    { // test template constructor and initialize_with
        typedef Int_vt< ::CGAL::Handle_policy_no_union> Int;
        Int i(5,6);
        assert( i == Int(11));
        Int j(5,6,7);
        assert( j == Int(18));
    }
    {
        //typedef Int_vt< ::CGAL::Handle_policy_in_place> Int; // That's supposed to fail
    }
    {
        typedef Int_vt< ::CGAL::Handle_policy_no_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( (std::ptrdiff_t)(ID_Number(i)) == i.id() );
    }
    {
        typedef Int_vt< ::CGAL::Handle_policy_union> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
    }
    {
        typedef Int_vt< ::CGAL::Handle_policy_union_and_reset> Int;
        Int i(5);
        Int j(5);
        Int k(6);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
        assert( i == j);
        assert( ! (i == k));
        assert( i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( k));
    }
    {
        typedef Int_vt< ::CGAL::Handle_policy_union> Int;
        Int i(5);
        Int j(5);
        Int k(5);
        Int l(5);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // pump up the union_size counter for j, k and l
        Int j1(j);
        assert( j1.test_identical_ptr( j));
        Int k1(k);
        Int k2(k);
        Int k3(k);
        Int l1(l);
        Int l2(l);
        Int l3(l);
        Int l4(l);
        Int l5(l);
        Int l6(l);
        Int l7(l);
        Int l8(l);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 2);
        assert( j.union_size() == 3);
        assert( k.union_size() == 5);
        assert( l.union_size() == 10);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link i to j
        assert( i == j);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 4);
        assert( j.union_size() == 4);
        assert( k.union_size() == 5);
        assert( l.union_size() == 10);
        assert( i.test_identical_ptr( j));
        assert( i.test_identical_ptr( j1));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link j to k
        assert( j == k);
        assert( i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 9);
        assert( k.union_size() == 9);
        assert( l.union_size() == 10);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( ! k.test_identical_ptr( l));
        // link k to l
        assert( k == l);
        assert( i.is_forwarding());
        assert( j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 8);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( ! j.test_identical_ptr( k));
        assert( ! j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
        // find j, links it to k and l
        assert( j.value() == 5);
        assert( i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 3);
        assert( j.union_size() == 19);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( ! i.test_identical_ptr( j));
        assert( ! i.test_identical_ptr( k));
        assert( ! i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
        // find i, links it to j, k and l
        assert( i.value() == 5);
        assert( ! i.is_forwarding());
        assert( ! j.is_forwarding());
        assert( ! k.is_forwarding());
        assert( ! l.is_forwarding());
        assert( i.union_size() == 19);
        assert( j.union_size() == 19);
        assert( k.union_size() == 19);
        assert( l.union_size() == 19);
        assert( i.test_identical_ptr( j));
        assert( i.test_identical_ptr( k));
        assert( i.test_identical_ptr( l));
        assert( j.test_identical_ptr( k));
        assert( j.test_identical_ptr( l));
        assert( k.test_identical_ptr( l));
    }
}

int main() {
    test_handle();
    test_handle_with_class_hierarchy();
    return EXIT_SUCCESS;
}
