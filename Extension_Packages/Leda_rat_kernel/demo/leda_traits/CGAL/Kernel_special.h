// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Kernel_special.h
// package       : 
// maintainer    : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : 
//
// ======================================================================

#ifndef CGAL_KERNEL_SPECIAL_H
#define CGAL_KERNEL_SPECIAL_H


#include <CGAL/basic.h>
#include <CGAL/Do_nothing.h>
#include <utility>
#include <iostream>



CGAL_BEGIN_NAMESPACE


template <class O1, 
          class PRE_OPERATION  = Do_nothing,
	  class POST_OPERATION = Do_nothing>
class Functor_enhancer
{
    O1 o1;
    
public:

    typedef Functor_enhancer<O1, PRE_OPERATION, POST_OPERATION>  Self;

    Functor_enhancer(const O1 &oo1 = O1())
	: o1(oo1) {}

    // the result type and arity must stay ...
    
    typedef typename O1::result_type result_type;
    typedef typename O1::Arity       Arity;

    template <class A1>
    result_type operator()(const A1 &a1) const
    { 
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1);
	    
	typename O1::result_type res = o1(a1);
	
	POST_OPERATION post_op;
	post_op(o1, a1, res);
	
	return res;
    }

    template <class A1, class A2>
    result_type
    operator()(const A1 &a1, const A2 &a2) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2);
	    
	typename O1::result_type res = o1(a1,a2);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, res);
	
	return res;
    }

    template <class A1, class A2, class A3>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3);
	    
	typename O1::result_type res = o1(a1,a2,a3);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, res);
	
	return res;
    }

    template <class A1, class A2, class A3, class A4>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3, a4);
	    
	typename O1::result_type res = o1(a1,a2,a3,a4);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, a4, res);
	
	return res; 
    }

    template <class A1, class A2, class A3, class A4, class A5>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3, a4, a5);
	    
	typename O1::result_type res = o1(a1,a2,a3,a4,a5);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, a4, a5, res);
	
	return res; 
    }

    template <class A1, class A2, class A3, class A4, class A5, class A6>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3, a4, a5, a6);
	    
	typename O1::result_type res = o1(a1,a2,a3,a4,a5,a6);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, a4, a5, a6, res);
	
	return res; 
    }

    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3, a4, a5, a6, a7);
	    
	typename O1::result_type res = o1(a1,a2,a3,a4,a5,a6,a7);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, a4, a5, a6, a7, res);
	
	return res; 
    }
    
    template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
    result_type
    operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	       const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
    {
        // perform an additional operation ...
	PRE_OPERATION  pre_op;
	pre_op(o1, a1, a2, a3, a4, a5, a6, a7, a8);
	    
	typename O1::result_type res = o1(a1,a2,a3,a4,a5,a6,a7,a8);
	
	POST_OPERATION post_op;
	post_op(o1, a1, a2, a3, a4, a5, a6, a7, a8, res);
	
	return res; 
    }    

};


// we inherit all geometric objects and constructions from K1

template <class K1, 
          class Pre_operation  = Do_nothing,
	  class Post_operation = Do_nothing>
class Kernel_special : public K1
{
    typedef K1     Kernel1;
    
public:

#define CGAL_special_predicate(X, Y) \
    typedef Functor_enhancer<typename K1::X, Pre_operation, Post_operation> X; \
    X Y() const { return X(K1::Y()); }
    
#define CGAL_special_construction(X, Y) \
    typedef Functor_enhancer<typename K1::X, Pre_operation, Post_operation> X; \
    X Y() const { return X(K1::Y()); }    

#define CGAL_Kernel_pred(Y,Z) CGAL_special_predicate(Y,Z)
#define CGAL_Kernel_cons(Y,Z) CGAL_special_construction(Y,Z)

public:

#include <CGAL/Kernel/interface_macros.h>
};

CGAL_END_NAMESPACE

#endif 
