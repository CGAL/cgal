// Copyright (c) 2001-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-inf.mpg.de>
//                 Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Lutz Kettner <kettner@mpi-inf.mpg.de>

#ifndef CGAL_HANDLE_WITH_POLICY_H
#define CGAL_HANDLE_WITH_POLICY_H

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <CGAL/type_traits.h>

#include <CGAL/assertions.h>
#include <CGAL/use.h>

#include <boost/mpl/if.hpp>

#include <cstddef>

#ifdef CGAL_USE_LEDA
#  if CGAL_LEDA_VERSION < 500
#    include <LEDA/memory.h>
#  else
#    include <LEDA/system/memory.h>
#  endif
#endif



namespace CGAL {

/*! \brief <tt>\#include <CGAL/Handle_with_policy.h></tt> for handles with policy
    parameter for reference counting and union-find strategy. Uses 
    \c LEDA_MEMORY if available. 

    There are two fundamentally different usages of this base class:

        - with a single representation class. In this case the handle
          manages allocation and deallocation and the type \c T can
          be an arbitrary type---the handle adds the necessary reference
          counter internally.

        - with a hierarchy of representation classes. Type \c T will be 
          the common base class of this hierarchy and it has to be derived
          itself from a specific base class, which can be accessed directly
          or generically from the policy class. The allocator in the 
          handle will not be used in this scenario, since the handle class
          does not allocate any representations. Instead, the handle class
          derived from this handle base class is allocating the different
          representations with the \c new operator. In this case, 
          the allocator in the base class of \c T is used.

    We give an example for each usage. See also the documentation
    of \c Handle_with_policy.

\b Example

We use a single representation class to store an integer. The second
constructor makes use of one of the forwarding template constructors
that simply forward their parameter list to the representation
constructors. They exist for up to ten parameters.  The third
constructor illustrates how the \c USE_WITH_INITIALIZE_WITH can be
used. It is useful if extensive computations are necessary before the
representation can be created.

\code
struct Int_rep {
    int val;
    Int_rep( int i = 0) : val(i) {}
    Int_rep( int i, int j) : val(i+j) {}
    Int_rep( int i, int j, int k) : val(i+j+k) {}
};

template < class Unify>
struct Int_t : public Handle_with_policy< Int_rep, Unify > {
    typedef Handle_with_policy< Int_rep, Unify > Base;
    Int_t( int i = 0) : Base( i) {}
    Int_t( int i, int j) : Base( i, j) {}     // template constructors
    Int_t( int i, int j, int k) : Base( Base::USE_WITH_INITIALIZE_WITH) {
        initialize_with( i, j + k);
    }
    int  value() const { return ptr()->val; }
    void set_value( int i) {
        copy_on_write();
        ptr()->val = i;
    }
    bool operator==( const Int_t<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            Base::unify(i);
        return equal;
    }
};
\endcode

\b Example

We use a class hierarchy of two representation classes: one base class
for representing one integer, and a derived class to represent an
additional integer. To also added virtual get and set functions to
make this example similar to the one above. 

We use the generic solution to pick the base class for \c Int_vrep
from the policy class. So all representations are class templates with
a policy and an allocator as parameter and the handle class
instantiates them. If this flexibility is not needed, one could derive
directly from the appropriate base class, i.e., \c
::CGAL::Reference_counted_hierarchy<Alloc> or \c
::CGAL::Reference_counted_hierarchy_with_union<Alloc>. \c Alloc is an
allocator of \c char's here.

\code
template <class Policy, class Alloc>
struct Int_vrep : public Policy::Hierarchy_base< Alloc>::Type {
    int val;
    virtual ::CGAL::Reference_counted_hierarchy<Alloc>* clone() { 
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
    virtual int get_val() const { return val + val2; }
    virtual void set_val( int i) { val = i - val2; }
    Int_vrep2( int i, int j) : Int_vrep<Policy,Alloc>(i), val2(j) {}
};

template < class Unify, class Alloc = CGAL_ALLOCATOR(char) >
struct Int_vt : public Handle_with_policy< Int_vrep<Unify,Alloc>, Unify > {
    typedef Handle_with_policy< Int_vrep<Unify,Alloc>, Unify > Base;
    Int_vt( int i = 0) : Base( new Int_vrep<Unify,Alloc>(i)) {}
    Int_vt( int i, int j) : Base( new Int_vrep2<Unify,Alloc>(i,j)) {}

    int  value() const { return ptr()->get_val(); }
    void set_value( int i) {
        copy_on_write();
        ptr()->set_val(i);
    }
    bool operator==( const Int_vt<Unify>& i) const {
        bool equal = (value() == i.value());
        if ( equal)
            Base::unify(i);
        return equal;
    }
};
\endcode

*/
//@{

// Forward declarations of HandlePolicy classes
class Handle_policy_in_place;
class Handle_policy_no_union;
class Handle_policy_union;
class Handle_policy_union_and_reset;

// Reference counted representation
// ================================

//! the base class for bodies of reference counted representations \c T.
template <class T_>
class Reference_counted {
public:
    typedef T_                          rep_type;
    typedef Reference_counted<rep_type> Self;
    typedef rep_type*                   Rep_pointer;
private:
    mutable unsigned int count;  // reference counter
    rep_type             rep;
public:
    Reference_counted() : count(1) {}
    Reference_counted( const rep_type& t) : count(1), rep(t) {}
    Reference_counted( const Self& r) : count(1), rep(r.rep) {}

    void clear() { rep = rep_type(); }
    Rep_pointer  base_ptr()  { return &rep; }
    void add_reference()     { ++count; }
    void remove_reference()  { --count; }
    bool is_shared() const   { return count > 1; }
    int  union_size() const  { return 1+count; }
    void add_union_size(int) {}
};

/*!\brief
 * Base class for bodies of reference counted representations \c T
 * with a forwarding pointer for identical representations.
 */
template <class T_>
class Reference_counted_with_forwarding {
public:
    typedef T_ rep_type;
    typedef Reference_counted_with_forwarding<rep_type> Self;
    typedef rep_type*  Rep_pointer;
    friend class Handle_policy_union;
    friend class Handle_policy_union_and_reset;
private:
    mutable unsigned int count;  // reference counter
    mutable Self*        next;   // forwarding pointer to valid rep or 0
    mutable int          u_size; // union set size incl this rep and its handle
    mutable rep_type     rep;
public:
    Reference_counted_with_forwarding()
        : count(1), next(0), u_size(2) {}
    Reference_counted_with_forwarding( const rep_type& t)
        : count(1), next(0), u_size(2), rep(t) {}
    Reference_counted_with_forwarding( const Self& r)
        : count(1), next(0), u_size(2), rep(r.rep) {}

    void clear() { rep = rep_type(); }
    Rep_pointer  base_ptr()    { return &rep; }
    void add_reference()       { ++count; }
    void remove_reference()    { --count; }
    bool is_shared() const     { return count > 1; }
    bool is_forwarding() const { return next != 0; }
    int  union_size() const    { return u_size; }
    void add_union_size(int a) {  
        CGAL_precondition( u_size + a > 0);
        u_size += a;
    }
};


struct Reference_counted_hierarchy_base {};


/*!\brief Base class for reference counted representations with a class 
 * hierarchy of different representations. Needs an allocator for \c char's
 * as parameter.
 */
template <class Allocator_  = CGAL_ALLOCATOR(char)>
class Reference_counted_hierarchy : public Reference_counted_hierarchy_base {
    // make sure it's always a char allocator
    typedef typename Allocator_::template rebind< char> Char_alloc_rebind;
    typedef typename Char_alloc_rebind::other   Char_allocator;

    static Char_allocator alloc;

public:
    void* operator new(size_t bytes) { return alloc.allocate( bytes); }
    void  operator delete(void* p, size_t bytes) { 
        alloc.deallocate((char*)p, bytes);
    }

public:
    typedef Allocator_ Allocator;
    typedef Reference_counted_hierarchy<Allocator> Self;
    typedef Self*  Rep_pointer;
private:
    mutable unsigned int count;  // reference counter
public:
    Reference_counted_hierarchy() : count(1) {}
    Reference_counted_hierarchy( const Self&) : count(1) {}

    Rep_pointer base_ptr()   { return this; }
    void add_reference()     { ++count; }
    void remove_reference()  { --count; }
    bool is_shared() const   { return count > 1; }
    int  union_size() const  { return 1+count; }
    void add_union_size(int) {}

    //! returns a copy of \c this. Can be implemented like
    //! <tt>return new Derived_type( *this);</tt>
    virtual Self* clone() = 0;
    //! the virtual destructor is essential for proper memory management here.
    virtual ~Reference_counted_hierarchy() {}
    //! can be used to minimize memory consumption once it is known that this
    //! representation is not used anymore and only needed to keep a fowarding
    //! pointer. One example would be cleaning up dynamically allocated
    //! data, or another example would be overwriting a \c leda::real with 
    //! a default constructed value to free its old expression tree. However,
    //! this function can also be savely ignored and kept empty.
    virtual void clear() {}
};

template <class Alloc>
typename Reference_counted_hierarchy<Alloc>::Char_allocator
    Reference_counted_hierarchy<Alloc>::alloc;

/*!\brief Base class for reference counted representations with a class 
 * hierarchy of different representations. Needs an allocator for \c char's
 * as parameter.
 */
template <class Allocator_  = CGAL_ALLOCATOR(char)>
class Reference_counted_hierarchy_with_union 
    : public Reference_counted_hierarchy<Allocator_> 
{
    friend class Handle_policy_union;
    friend class Handle_policy_union_and_reset;
public:
    typedef Allocator_ Allocator;
    typedef Reference_counted_hierarchy_with_union<Allocator> Self;
private:
    mutable Self*        next;   // forwarding pointer to valid rep or 0
    mutable int          u_size; // union set size incl this rep and its handle
public:
    Reference_counted_hierarchy_with_union() : 
        Reference_counted_hierarchy<Allocator_>(), next(0), u_size(2) {}
    bool is_forwarding() const { return next != 0; }
    int  union_size() const    { return u_size; }
    void add_union_size(int a) {  
        CGAL_precondition( u_size + a > 0);
        u_size += a;
    }
};


// Handle for reference counted representation
// ===========================================

namespace Intern {
    // Some helper classes to select representation between single class
    // representations and class hierarchy representations.

    // the representation type including a reference counter. 
    // The handle allocates objects of this type. This is the version
    // for the single representation type.
    template <class T, int HandleHierarchyPolicy>
    struct Rep_bind_reference_counted {
        typedef Reference_counted<T> Rep;
    };

    // the representation type including a reference counter. 
    // The handle allocates objects of this type. This is the version
    // for the class hierarchy of representation types.
    template <class T>
    struct Rep_bind_reference_counted<T, true> {
        typedef T Rep;
    };

    // the two versions for Reference_counted_with_forwarding
    template <class T, int HandleHierarchyPolicy>
    struct Rep_bind_reference_counted_with_forwarding {
        typedef Reference_counted_with_forwarding<T> Rep;
    };

    // the representation type including a reference counter. 
    // The handle allocates objects of this type. This is the version
    // for the class hierarchy of representation types.
    template <class T>
    struct Rep_bind_reference_counted_with_forwarding<T, true> {
        Rep_bind_reference_counted_with_forwarding() {
            // make sure we derived from the right type 
            typedef typename T::Allocator Alloc;
            typedef ::CGAL::Reference_counted_hierarchy_with_union<Alloc> 
                Reference_counted_hierarchy_with_union;
            CGAL_USE_TYPE(Reference_counted_hierarchy_with_union);
            CGAL_static_assertion((
              ::CGAL::is_same_or_derived< Reference_counted_hierarchy_with_union, T >::value ));
        }
        typedef T Rep;
    };

}

/*! \brief Policy class for \c Handle_with_policy that stores the
    representation directly without reference counting and without dynamic 
    memory allocation, is actually \e not a model of the \c HandlePolicy
    concept, but can be  used instead of one. It selects a different
    specialized implementation of \c Handle_with_policy. It works only with
    the single representation type, not with a class hierarchy of 
    representation types since they need the pointer in the handle
    for the polymorphy.
*/
class Handle_policy_in_place {};

/*!\brief
 * Policy class for \c Handle_with_policy<T> that ignores unifying of
 * identical representations \c T, is a model of the \c HandlePolicy concept.
 */
class Handle_policy_no_union {
public:
    /*!\brief
     * A rebind mechanism to create the representation type.
     */
    template <class T, int hierarchy>
    struct Rep_bind {
        //! the representation type including a reference counter. 
        //! The handle allocates objects of this type.
        typedef typename 
            Intern::Rep_bind_reference_counted<T,hierarchy>::Rep Rep;
    };

    /*!\brief
     * A rebind mechanism to access the base class for class hierarchies
     * of representations. 
     * 
     * The base classes can be used directly, but this
     * rebind mechamism allows the implementation of handle-rep classes
     * that are parameterized with the policy class only and adapt to
     * the necessary base class.
     */
    template <class Alloc>
    struct Hierarchy_base {
        //! type that can be used as base class for the representation type.
        typedef Reference_counted_hierarchy<Alloc> Type;
    };
    
    /*! \brief unifies the representations of the two handles \a h and \a g.
     *  The effect is void here.
     * 
     * \pre The representations represent the same value and one could be 
     *  replaced by the other.
    */
    template <class H>
    static void unify( const H& h, const H& g) {
        (void)h; // avoid warnings for unused parameters
        (void)g; // but keep the names in the definition for the doc.
    }

    //! finds the currently valid representation for the handle \a h
    //! and returns a pointer to its stored value of type \a T.
    template <class H>
    static typename H::Rep_pointer find( const H& h) {
        return h.ptr_->base_ptr();
    }
};

/*!\brief
 * Policy class for \c Handle_with_policy that implements unifying of
 * identical representations \c T with trees and path compression, is a
 * model of the \c HandlePolicy concept.
 */
class Handle_policy_union {
public:
    /*!\brief
     * A rebind mechanism to create the representation type.
     */
    template <class T, int hierarchy>
    struct Rep_bind {
        //! this default constructor contains some compile-time checks.
        Rep_bind() {
	  //Intern::Rep_bind_reference_counted_with_forwarding<T, hierarchy>
	  //     check;
          //  (void)check;
	  (void)Intern::Rep_bind_reference_counted_with_forwarding<T, hierarchy>();
        }
        //! the representation type including a reference counter. 
        //! The handle allocates objects of this type.
        typedef typename Intern::Rep_bind_reference_counted_with_forwarding<T,
            hierarchy>::Rep Rep;
    };

    /*!\brief
     * A rebind mechanism to access the base class for class hierarchies
     * of representations. 
     * 
     * The base classes can be used directly, but this
     * rebind mechamism allows the implementation of handle-rep classes
     * that are parameterized with the policy class only and adapt to
     * the necessary base class.
     */
    template <class Alloc>
    struct Hierarchy_base {
        //! type that can be used as base class for the representation type.
        typedef Reference_counted_hierarchy_with_union<Alloc> Type;
    };
    
    /*! \brief unifies the representations of the two handles \a h and \a g.
        Performs union.
        \pre The representations represent the same value and one can be 
        replaced by the other. The handles \a h and \a g are already
        the representatives found by the find operation and \a h is not
        equal to \a g. The tree representing the union of \a h has size
        not smaller than the corresponding tree size of \a g.
    */
    template <class H>
    static void unify_large_small( const H& h, const H& g) {
        typename H::Rep* hrep = h.ptr_;
        typename H::Rep* grep = g.ptr_;
        CGAL_precondition( ! grep->is_forwarding());
        CGAL_precondition( hrep->union_size() >= grep->union_size());
        grep->add_union_size(-1);
        // make g point to h's rep.
        if ( grep->is_shared()) {
            // grep survises the loss of one reference 
            // and hrep gets one more reference
            grep->remove_reference();
            hrep->add_reference();
            hrep->add_union_size( grep->union_size());
            grep->next = hrep;
        } else {
            g.delete_rep( grep); // did not survive loss of handle g
        }
        // redirect handle g and incr. hrep's counter
        g.ptr_ = hrep;
        hrep->add_reference();
        hrep->add_union_size(1);
    }

    /*! \brief unifies the representations of the two handles \a h and \a g.
        Performs union with path compression.
        \pre The representations represent the same value and one can be 
        replaced by the other.
    */
    template <class H>
    static void unify( const H& h, const H& g) {
        if ( find(h) !=  find(g)) {
            if ( h.ptr_->union_size() > g.ptr_->union_size())
                unify_large_small( h, g); // make g point to h's rep.
            else
                unify_large_small( g, h); // make h point to g's rep.
        }
    }

    /*! \brief finds the currently valid representation for the handle \a h
        and returns a pointer to its stored value of type \a T. Performs
        path-compression to speed-up later union operations.
    */
    template <class H>
    static typename H::Rep_pointer find( const H& h) {
        typedef typename H::Rep Rep;
        if ( h.ptr_->is_forwarding()) {
            // find new valid representation
            Rep* new_rep = h.ptr_;
            while ( new_rep->next != 0)
                new_rep = static_cast<Rep*>(new_rep->next);
            // path compression: assign new rep to all reps seen on the path
            // update reference count properly: all reps on the path loose
            // one reference, and the new_rep gains all of them unless
            // the rep on the path get actually deleted.
            Rep* rep = h.ptr_;
            while ( rep != new_rep) {
                Rep* tmp = static_cast<Rep*>(rep->next);
                if ( rep->is_shared()) {
                    // rep survives the loss of one reference 
                    // and new_rep gets one more reference
                    rep->remove_reference();
                    if ( tmp != new_rep) {
                        // re-link rep to the new_rep
                        rep->next = new_rep;
                        new_rep->add_reference();
                    }
                } else {
                    h.delete_rep( rep); // we have to delete the current rep
		    tmp->remove_reference();
                }                
                rep = tmp;
            }
            // hook h to new_rep
            h.ptr_ = new_rep;
            new_rep->add_reference();
        }
        return h.ptr_->base_ptr();
    }
};

/*!\brief Policy class for \c Handle_with_policy that implements unifying of
 * identical representations \c T with trees and path compression. 
 * 
 * It also 
 * sets the unused representation immediately to the default constructed
 * representation \c T(), which can help to free memory if the 
 * representation is dynamically allocated and potentially large, e.g.,
 * \c leda::real. This class is a model of the \c HandlePolicy concept.
 */
class Handle_policy_union_and_reset {
public:
    /*!\brief
     * A rebind mechanism to create the representation type.
     */
    template <class T, int hierarchy>
    struct Rep_bind {
        //! this default constructor contains some compile-time checks.
        Rep_bind() {
	  //Intern::Rep_bind_reference_counted_with_forwarding<T, hierarchy>
	  //     check;
	  // (void)check;
	  (void)Intern::Rep_bind_reference_counted_with_forwarding<T, hierarchy>();
        }
        //! the representation type including a reference counter. 
        //! The handle allocates objects of this type.
        typedef typename Intern::Rep_bind_reference_counted_with_forwarding<T,
            hierarchy>::Rep Rep;
    };

    /*!\brief
     * A rebind mechanism to access the base class for class hierarchies
     * of representations. 
     *
     * The base classes can be used directly, but this
     * rebind mechamism allows the implementation of handle-rep classes
     * that are parameterized with the policy class only and adapt to
     * the necessary base class.
     */
    template <class Alloc>
    struct Hierarchy_base {
        //! type that can be used as base class for the representation type.
        typedef Reference_counted_hierarchy_with_union<Alloc> Type;
    };

    // abbreviation to re-use its implementation below.
    typedef Handle_policy_union U;

    /*! \brief unifies the representations of the two handles \a h and \a g.
        Performs union with path compression and assigns a default
        constructed value of the representation type \c Rep to the
        superfluous representation.
        \pre The representations represent the same value and one can be 
        replaced by the other.
    */
    template <class H>
    static void unify( const H& h, const H& g) {
        if ( find(h) !=  find(g)) {
            if ( h.ptr_->union_size() > g.ptr_->union_size()) {
                // reset representation in g to default construction of T
                if ( g.ptr_->is_shared())
                    g.ptr_->clear();
                U::unify_large_small( h, g); // make g point to h's rep.
            } else {
                // reset representation in h to default construction of T
                if ( h.ptr_->is_shared())
                    h.ptr_->clear();
                U::unify_large_small( g, h); // make h point to g's rep.
            }
        }
    }

    /*! \brief finds the currently valid representation for the handle \a h
        and returns a pointer to its stored value of type \a T. Performs
        path-compression to speed-up later union operations.
    */
    template <class H>
    static typename H::Rep_pointer find( const H& h) { return U::find(h); }
};


/*! \brief the base class for handles of reference counted representations of
    \c T. 

    There are two fundamentally different usages of this base class:

        - with a single representation class. In this case the handle
          manages allocation and deallocation and the type \c T can
          be an arbitrary type---the handle adds the necessary reference
          counter internally.

        - with a hierarchy of representation classes. Type \c T will be 
          the common base class of this hierarchy and it has to be derived
          itself from either \c ::CGAL::Reference_counted_hierarchy or
          \c ::CGAL::Reference_counted_hierarchy_with_union, both parameterized 
          with an allocator. The allocator in the handle will not be used in 
          this scenario, since the handle class does not allocate any 
          representations. Instead, the handle class derived from this handle
          base class is allocating the different representations with the
          \c new operator. In this case, the allocator in the base class
          of \c T is used.

    The handle class distinguishes between these two alternative
    usages by checking if \c T is derived from one of the two base
    classes mentioned for the second alternative. If not, it picks the
    first alternative.

    In the second alternative, the correct base class, \c 
    ::CGAL::Reference_counted_hierarchy_with_union, has to be used
    if the policy class is one of \c class Handle_policy_union r \c 
    Handle_policy_union_and_reset. Otherwise, the other base class can 
    be used to save space.

    The policy class \c Handle_policy_in_place is incompatible with the class
    hierarchy for representation classes since the pointer in the
    handle class would be missing.

    The dependency of the base classes for \c T and the policy classes
    is also encoded in the policy classes and can be used to write
    generic handle-rep scheme classes. To do that one can derive \c T
    from the expressions \c Policy::Hierarchy_base<Alloc>::Type
    assuming that \c Policy is the handle policy and \c Alloc is the
    allocator. Btw, the allocator is used as an allocator of character
    arrays here.

    \see \link Handle Handle for Reference Counting\endlink for 
    an example for each of the two alternative usages.

    The template parameters are:
        - \b T: is one of the two following:
            - an arbitrary type but it must be a model of the
              \c DefaultConstructible concept if the default constructor
              of the handle is used.
            - a type derived from \c Reference_counted_hierarchy<Alloc> or
              \c Reference_counted_hierarchy_with_union<Alloc> implementing
              their virtual member function interface, namely a \c clone()
              function.

        - \b HandlePolicy: a model of the \c HandlePolicy concept or the
              \c Handle_policy_in_place class template that selects a specialized
              implementation without reference counting. Has the
              default \c Handle_policy_no_union.

        - \b Allocator_: a model of the \c Allocator concept,
          has the default \c CGAL_ALLOCATOR(T).

*/
template <class T_, 
          class HandlePolicy = Handle_policy_no_union, 
          class Allocator_ = CGAL_ALLOCATOR(T_)>
class Handle_with_policy {
public:

    //! first template parameter
    typedef T_ Handled_type;

    //! the handle type itself.
    typedef Handle_with_policy< Handled_type, HandlePolicy, Allocator_>    Self;

    //! the instantiated model of the \c HandlePolicy concept.
    typedef HandlePolicy                Handle_policy;

    //! the allocator type.
    typedef Allocator_                  Allocator;

    enum { is_class_hierarchy  = 
        ::CGAL::is_same_or_derived< Reference_counted_hierarchy_base, Handled_type>::value };
        
    typedef typename Handle_policy::template Rep_bind< Handled_type, is_class_hierarchy > Bind;
    
    // instantiate Rep_bind to activate compile time check in there
    static Bind bind;
    
    // Define type that is used for function matching 
    typedef typename ::boost::mpl::if_c< 
         is_class_hierarchy, 
           ::CGAL::Tag_true, 
           ::CGAL::Tag_false >::type 
         Class_hierarchy;

    //! the internal representation, i.e., \c T plus a reference count
    //! (if needed), or just \c T if we derived from the base class to
    //! support a class hierarchy for the representations.    
    typedef typename Bind::Rep  Rep;

    typedef typename Rep::Rep_pointer  Rep_pointer;

    typedef typename Allocator_::template rebind<Rep>::other  Rep_allocator;


    //! integer type for identifying a representation.
    typedef std::ptrdiff_t              Id_type;

    friend class Handle_policy_no_union;
    friend class Handle_policy_union;
    friend class Handle_policy_union_and_reset;
private:
    mutable Rep*  ptr_;

    // We have to distinguish between allocating single representations
    // and where we have a class hierarchy of representations, where the
    // user is responsible for allocating the first representations
    // and we can just \c clone and delete them.
    static Rep_allocator allocator;

    static Rep* new_rep( const Rep& rep) { 
        CGAL_static_assertion( !(
           ::CGAL::is_same_or_derived< Reference_counted_hierarchy_base, Handled_type >::value ));
        Rep* p = allocator.allocate(1);
        allocator.construct(p, rep);
        return p;
    }
    static void delete_rep( Rep* p, ::CGAL::Tag_false ) {
        allocator.destroy( p);
        allocator.deallocate( p, 1);
    }
    static void delete_rep( Rep* p, ::CGAL::Tag_true ) {
        delete p;
    }
    static void delete_rep( Rep* p) { delete_rep(p, Class_hierarchy()); }

    static Rep* clone_rep( Rep* p, ::CGAL::Tag_false ) {
        return new_rep( *p);
    }
    static Rep* clone_rep( Rep* p, ::CGAL::Tag_true ) {
        return static_cast<Rep*>(p->clone());
    }
    static Rep* clone_rep( Rep* p) { return clone_rep( p, Class_hierarchy()); }

    void remove_reference() {
        // cleans up the possible chain of forwarding reps
        Handle_policy::find( *this);
        if ( ! is_shared()) {
            delete_rep( ptr_);
        } else {
            ptr_->remove_reference();
            ptr_->add_union_size( -1);
        }
    }

    template <class TT>
    Rep* make_from_single_arg( const TT& t, ::CGAL::Tag_false ) {
        return new_rep( Rep( Handled_type(t)));
    }
    template <class TT>
    Rep* make_from_single_arg( TT t, ::CGAL::Tag_true ) {
      //Bind bind_; // trigger compile-time check
      // (void)bind_;
      (void)Bind(); // shouldn't this be enough to trigger?
        return t; // has to be a pointer convertible to Rep*
    }

protected:
    //! protected access to the stored representation
    Handled_type*       ptr()       { return static_cast<Handled_type*>(Handle_policy::find(*this));}
    //! protected access to the stored representation
    const Handled_type* ptr() const { 
        return static_cast<const Handled_type*>(Handle_policy::find( *this));
    }

    //! unify two representations. \pre The two representations describe
    //! the same value  and one can be replaced by the other, i.e., the 
    //! values are immutable, or protected from changes with \c copy_on_write()
    //! calls!
    void unify( const Self& h) const { Handle_policy::unify( *this, h); }

    //! can be called before modifying a shared representation
    //! to get an own copy of the representation which avoids effecting the
    //! other sharing handles. Does nothing if representation is actually
    //! not shared.
    void copy_on_write() {
        Handle_policy::find( *this);
        if ( is_shared() ) {
            Rep* tmp_ptr = clone_rep( ptr_);
            ptr_->remove_reference();
            ptr_->add_union_size( -1);
            ptr_ = tmp_ptr;
        }
    }

    //! used with special protected constructor
    enum Use_with_initialize_with {
        USE_WITH_INITIALIZE_WITH //!< used with special protected constructor
                                 //!< of \c Handle_with_policy.
    };

    //! special constructor, postpones the construction of the representation
    //! to one of the \c initialize_with() functions. An object is in an
    //! invalid state (and will report a failed precondition later) if
    //! it is not initialized with an \c initialize_with() function call
    //! after this constructor. Applicable for single representation but 
    //! also for a class hierarchy of representations.
    Handle_with_policy( Use_with_initialize_with) : ptr_( 0) {}

    //! constructor used for class hierarchies of representations, where
    //! the handle class derived from this handle creates the different 
    //! representations itself with the \c new operator. Except for this
    //! constructor, the the one with the \c Use_with_initialize_with
    //! argument, and the single argument template constructor no other 
    //! constructor will work for class hierarchies of representations.
    Handle_with_policy( Rep* p) : ptr_( p) {
        CGAL_static_assertion((
           ::CGAL::is_same_or_derived< Reference_counted_hierarchy_base, Handled_type >::value ));
        //Bind bind_; // trigger compile-time check
        //(void)bind_;
	(void)Bind();
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used. Applicable for a
    //! class hierarchy of representations only, where the derived handle class
    //! created the representation \c p with the \c new operator. No other
    //! version of \c initialize_with is applicable in this case except
    //! the template version with one argument.
    void initialize_with( Rep* p) {
        CGAL_static_assertion((
           ::CGAL::is_same_or_derived< Reference_counted_hierarchy_base, Handled_type >::value ));
        //Bind bind_; // trigger compile-time check
        //(void)bind_;
	(void)Bind();
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = p;
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    //! In case of the class hierarchy of representation classes,
    //! this function is also chosen for pointers to newly allocated
    //! representations that are types derived from \c T. In that case,
    //! the pointer is just assigned to the internal pointer.
    template <class T1>
    void initialize_with( const T1& t1) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = make_from_single_arg( t1, Class_hierarchy());
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2>
    void initialize_with( const T1& t1, const T2& t2) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4,t5)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4,t5,t6)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4,t5,t6,t7)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7, class T8>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7, const T8& t8) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4,t5,t6,t7,t8)));
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7, class T8, class T9>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7, const T8& t8, const T9& t9) {
        CGAL_precondition_msg( ptr_ == 0, "Handle_with_policy::initialize_with(): the "
                         "representation has already been initialized.");
        ptr_ = new_rep( Rep( Handled_type(t1,t2,t3,t4,t5,t6,t7,t8,t9)));
    }

public:
    //! default constructor.
    Handle_with_policy() : ptr_( new_rep( Rep())) {}

    //! copy constructor, increments reference count.
    Handle_with_policy(const Self& h) {
        CGAL_precondition_msg( h.ptr_ != 0, "Handle_with_policy::Handle_with_policy( Self): probably "
                         "used special protected constructor and not the "
                         "'initialize_with()' function.");
        Handle_policy::find( h);
        ptr_ = h.ptr_;
        ptr_->add_reference();
        ptr_->add_union_size( 1);
    }

    //! forwarding constructor passing its parameter to the representation
    //! constructor. In case of the class hierarchy of representation classes,
    //! this constructor is also chosen for pointers to newly allocated
    //! representations that are types derived from \c T. In that case,
    //! the pointer is just assigned to the internal pointer.
    template <class T1>
    explicit Handle_with_policy( const T1& t) 
        : ptr_( make_from_single_arg( t, Class_hierarchy())) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2>
    Handle_with_policy( const T1& t1, const T2& t2) : ptr_( new_rep( Rep( Handled_type( t1, t2)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4, t5)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4, t5, t6)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4, t5, t6, t7)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7, class T8>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7, const T8& t8) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4, t5, t6, t7, t8)))) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7, class T8, class T9>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7, const T8& t8,
            const T9& t9) 
        : ptr_( new_rep( Rep( Handled_type( t1, t2, t3, t4, t5, t6, t7, t8, t9)))) {}

    //! destructor, decrements reference count.
    ~Handle_with_policy() {
      //Bind bind_; // trigger compile-time check
      //(void)bind_;
      (void)Bind();
        CGAL_precondition_msg( ptr_ != 0, "Handle_with_policy::~Handle_with_policy(): probably used "
                         "special protected constructor and not the "
                         "'initialize_with()' function.");
        remove_reference();
    }

    //! assignment, updates reference count correspondingly.
    Self& operator=( const Self& h) {
        CGAL_precondition_msg( h.ptr_ != 0, "Handle_with_policy::operator=(): probably "
                         "used special protected constructor and not the "
                         "'initialize_with()' function.");
        Handle_policy::find( h);
        h.ptr_->add_reference();
        h.ptr_->add_union_size( 1);
        remove_reference();
        ptr_ = h.ptr_;
        return *this;
    }

    //! returns \c true if both share the same representation.
    bool is_identical( const Self& h) const { return ptr() == h.ptr(); }

    //! returns a unique id value. Two handles share their representation
    //! is their id values are identical.
    Id_type id() const { return reinterpret_cast<Id_type>(&*ptr()); }

    //! returns true if the representation is shared, i.e., the reference
    //! counter is greater than one.
    bool is_shared() const { return ptr_->is_shared(); }

    //! returns \c true if the representation is actually forwarding to
    //! another equivalent representation (happens only with the
    //! union-find policies).
    bool is_forwarding() const { return ptr_->is_forwarding(); }

    //! returns the size of the union set including all reference counts that
    //! have been accumulated so far for this representation.
    int  union_size() const { return ptr_->union_size(); }

    // backwards compatible
    bool identical( const Self& h) const { return is_identical(h); }

#ifdef CGAL_HANDLE_WITH_POLICY_INTERNAL_TEST
    // provide access to pointer for testing only!!
    const Rep* test_ptr() const { return ptr_; }
    // provide access to pointer for testing only!!
    bool test_identical_ptr( const Self& h) const { return ptr_ == h.ptr_; }
#endif // CGAL_HANDLE_WITH_POLICY_INTERNAL_TEST
};

// instantiate Rep_bind to activate compile time check in there
template <class T, class Policy, class Alloc>
typename Handle_with_policy<T,Policy,Alloc>::Bind Handle_with_policy<T,Policy,Alloc>::bind;


//! alternative syntax for \c h.id() to allow use with LEDA
/*! This is only provided for \c Handle_policy_no_union because
 *  ID numbers have to be fixed throughout an object's lifetime.
 */
template <class T, class A>
typename Handle_with_policy<T, Handle_policy_no_union, A>::Id_type
ID_Number(const Handle_with_policy<T, Handle_policy_no_union, A>& h)
    { return h.id(); }

template <class T, class Policy, class Alloc>
typename Handle_with_policy<T, Policy, Alloc>::Rep_allocator 
    Handle_with_policy<T, Policy, Alloc>::allocator;


/*! \brief specialization of the base class for handles for non-reference
    counted representations.
    Uses \c LEDA_MEMORY if available.
*/
template <class T_, class Allocator_>
class Handle_with_policy<T_, Handle_policy_in_place, Allocator_> {
public:

    //! first template paramter
    typedef T_ Handled_type;

    //! the handle type itself.
    typedef Handle_with_policy< Handled_type, Handle_policy_in_place, Allocator_>   Self;

    //! the model of the \c HandlePolicy concept.
    typedef Handle_policy_in_place                           Handle_policy;

    //! the allocator type.
    typedef Allocator_                                Allocator;

    //! identify \c T with the internal representation \c Rep.
    typedef Handled_type                              Rep;

    //! integer type for identifying a representation.
    typedef std::ptrdiff_t                            Id_type;
private:
    // store the rep in place
    Rep  rep;

protected:
    //! protected access to the stored representation
    Handled_type*       ptr()       { return &rep; }
    //! protected access to the stored representation
    const Handled_type* ptr() const { return &rep; }

    //! unify two representations, a null op here.
    void unify( const Self&) const {}

    //! can be called before modifying a shared representation
    //! to get an own copy of the representation, a null op here.
    void copy_on_write() {}

    //! used with special protected constructor
    enum Use_with_initialize_with {
        USE_WITH_INITIALIZE_WITH //!< used with special protected constructor
    };

    //! special constructor, postpones the construction of the representation
    //! to one of the \c initialize_with() functions. Requires default
    //! constructor for \c T.
    Handle_with_policy( Use_with_initialize_with) {}
    
    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1>
    void initialize_with( const T1& t1) { rep = Rep(t1); }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2>
    void initialize_with( const T1& t1, const T2& t2) { rep = Rep(t1,t2); }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3) {
        rep = Rep(t1,t2,t3);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4) {
        rep = Rep(t1,t2,t3,t4);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5) {
        rep = Rep(t1,t2,t3,t4,t5);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6) {
        rep = Rep(t1,t2,t3,t4,t5,t6);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7) {
        rep = Rep(t1,t2,t3,t4,t5,t6,t7);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7, class T8>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7, const T8& t8) {
        rep = Rep(t1,t2,t3,t4,t5,t6,t7,t8);
    }

    //! initializes the representation after the constructor from 
    //! \c USE_WITH_INITIALIZE_WITH has been used.
    template <class T1, class T2, class T3, class T4, class T5, class T6,
              class T7, class T8, class T9>
    void initialize_with( const T1& t1, const T2& t2, const T3& t3,
                          const T4& t4, const T5& t5, const T6& t6,
                          const T7& t7, const T8& t8, const T9& t9) {
        rep = Rep(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    }

public:
    //! default constructor.
    Handle_with_policy() {}

    //! copy constructor.
    Handle_with_policy(const Self& h) : rep( h.rep) {}

    //! forwarding constructor passing its parameter to the representation
    //! constructor.
    template <class T1>
    explicit Handle_with_policy( const T1& t) : rep( Rep(t)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2>
    Handle_with_policy( const T1& t1, const T2& t2) : rep( Rep(t1,t2)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3) : rep( Rep(t1,t2,t3)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4) 
        : rep( Rep( t1, t2, t3, t4)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5) 
        : rep( Rep( t1, t2, t3, t4, t5)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6) 
        : rep( Rep( t1, t2, t3, t4, t5, t6)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7) 
        : rep( Rep( t1, t2, t3, t4, t5, t6, t7)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7, class T8>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7, const T8& t8) 
        : rep( Rep( t1, t2, t3, t4, t5, t6, t7, t8)) {}

    //! forwarding constructor passing its parameters to the representation
    //! constructor.
    template <class T1, class T2, class T3, class T4, class T5, class T6, 
              class T7, class T8, class T9>
    Handle_with_policy( const T1& t1, const T2& t2, const T3& t3, const T4& t4,
            const T5& t5, const T6& t6, const T7& t7, const T8& t8,
            const T9& t9) 
        : rep( Rep( t1, t2, t3, t4, t5, t6, t7, t8, t9)) {}

    //! returns \c true if both share the same representation.
    bool is_identical( const Self& h) const { return this == &h; }

    //! returns a unique id value. Two handles share their representation
    //! is their id values are identical.
    Id_type id() const { return ptr() - static_cast<Handled_type const*>(0); }

    //! returns \c false since the representation is not shared for
    //! this specialization.
    bool is_shared() const { return false; }

    //! returns \c false since the representation is not forwarding for
    //! this specialization.
    bool is_forwarding() const { return false; }

    //! returns \c 1 as the union size for this specialization.
    int  union_size() const { return 1; }

    // backwards compatible
    bool identical( const Self& h) const { return is_identical(h); }

#ifdef CGAL_HANDLE_WITH_POLICY_INTERNAL_TEST
    // provide access to pointer for testing only!!
    const Rep* test_ptr() const { return *rep; }
    // provide access to pointer for testing only!!
    bool test_identical_ptr( const Self& h) const { return this == &h; }
#endif // CGAL_HANDLE_WITH_POLICY_INTERNAL_TEST

#ifdef CGAL_USE_LEDA
    LEDA_MEMORY( Self)
#endif
};

template <class T, class HandlePolicy, class Allocator>
inline bool identical(const Handle_with_policy<T,HandlePolicy,Allocator> &h1, const Handle_with_policy<T,HandlePolicy,Allocator> &h2) { return h1.is_identical(h2); }


/*\brief
 * This class' function call operator test whether one handle's \c id is
 * less than the \c id of the other handle.
 *
 * "Less" is defined in terms of the second template argument,
 * which defaults to \c std::less<Handle::Id_type>
 */
template <class Handle, class Less = std::less<typename Handle::Id_type> >
class Handle_id_less_than {
public:
    //! result_type
    typedef bool result_type;
    //! type of first argument
    typedef Handle first_argument_type;
    //! type of second argument
    typedef Handle second_argument_type;
    //! returns \c true iff \c h1.id() < \c h2.id()
    bool operator () (Handle h1, Handle h2) {
        Less is_less;
        return is_less(h1.id(), h2.id());
    }
    //! returns \c true iff \c h1.id() < \c h2.id()
    bool operator () (Handle h1, Handle h2) const {
        Less is_less;
        return is_less(h1.id(), h2.id());
    }
};


//@}

} //namespace CGAL

#endif // CGAL_HANDLE_WITH_POLICY_H
