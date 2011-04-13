// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 09
//
// file          : include/CGAL/Tree_base.h
// package       : SearchStructures (2.54)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// source        : include/CGAL/Tree_base.h 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#ifndef __CGAL_Tree_base_d__
#define __CGAL_Tree_base_d__

#include <iterator>
#include <iostream>
#include <functional>
#include <list>
#include <vector>
#include <CGAL/assertions.h>
#include <CGAL/Tree_assertions.h>

#ifndef  USE_ARGUMENT
#define  USE_ARGUMENT( X )  (void)(X);
#endif  
#ifndef  USE_VARIABLE
#define  USE_VARIABLE( X )  (void)(&(X));
#endif  

#ifndef TREE_BASE_NULL
#define TREE_BASE_NULL 0
#endif

#ifndef CGAL__static_cast
#define CGAL__static_cast(TYPE,EXPR) (TYPE)(EXPR)
#endif

#ifndef CGAL__const_cast
#define CGAL__const_cast(TYPE,EXPR) (TYPE)(EXPR)
#endif

#ifndef CGAL__reinterpret_cast
#define CGAL__reinterpret_cast(TYPE,EXPR) (TYPE)(EXPR)
#endif


#define stlvector

CGAL_BEGIN_NAMESPACE

template <class InIt>
int   count_elements__C( const InIt  first, const InIt  last )
{
  InIt z=first;
  int i=0;
  
  while ( z++ != last ) {  
    i++;
  }

  return  i;
}



//link type definition of an ordinary vertex of the tree
struct tree_node_base {
  void *parent_link;
  void *left_link;
  void *right_link;
  tree_node_base(){parent_link=0; left_link=0; right_link=0;}
};


// -------------------------------------------------------------------
// pure virtual abstract base class.
// Designed according to the Prototype Design Pattern 
// A tree class has to be derived from this class.

template <class C_Data, class C_Window>
class tree_base
{

protected:
  tree_base(tree_base const &); // prevent access
  void operator= (tree_base const &); // prevent access

public:
  typedef double vit;
  typedef int lit;
  typedef int lbit;
  typedef double vbit;
  typedef char oit;
  //  typedef std::vector<C_Data>::iterator vit;
  //typedef std::list<C_Data>::iterator lit;
  //typedef std::back_insert_iterator<lit>  lbit;
  //typedef std::back_insert_iterator<vit>  vbit;
  typedef tree_base<C_Data, C_Window> tree_base_type;
  tree_base() {}
  virtual ~tree_base() {}

  // 'clone()' returns an object which can be used as argument to 'delete'
  virtual tree_base<C_Data, C_Window>  *clone() const = 0;
  //virtual tree_base_type   *clone() const = 0;

  // 'make_tree()' returns an object which can be used as argument to 'delete'
  virtual bool make_tree(typename std::list<C_Data>::iterator& beg, 
			 typename std::list<C_Data>::iterator& end,
			 lit *dummy=0) =0;
#ifdef stlvector
  virtual bool make_tree(typename std::vector<C_Data>::iterator& beg, 
			 typename std::vector<C_Data>::iterator& end,
			 vit *dummy=0) =0;
#endif
#ifdef carray
  virtual bool make_tree(C_Data *beg, 
                         C_Data *end) =0;
#endif
  virtual std::back_insert_iterator< std::list<C_Data> > 
    window_query(C_Window const &win,  std::back_insert_iterator<
		 std::list<C_Data> > out,lbit *dummy=0 ) = 0; 
  virtual std::back_insert_iterator< std::vector<C_Data> >
    window_query(C_Window const &win,  std::back_insert_iterator<
		  std::vector<C_Data> > out,vbit *dummy=0) = 0; 
#ifdef carray
  virtual C_Data * window_query( C_Window const &win, 
			        C_Data * out) = 0; 
#endif
#ifdef ostreamiterator
  typedef std::ostream_iterator< C_Data> oit;
  virtual  std::ostream_iterator< C_Data> window_query( C_Window const &win, 
				     std::ostream_iterator< C_Data> out,
					oit *dummy=0	       ) = 0; 
#endif
  virtual  std::back_insert_iterator< std::list< C_Data> > 
    enclosing_query( C_Window const &win,  std::back_insert_iterator<
		     std::list< C_Data> > out, lbit *dummy=0 ) = 0; 
  virtual  std::back_insert_iterator< std::vector< C_Data> > 
    enclosing_query( C_Window const &win,  std::back_insert_iterator<
		     std::vector< C_Data> > out,vbit *dummy=0 ) = 0; 
#ifdef carray
  virtual   C_Data * enclosing_query( C_Window const &win, 
				    C_Data *out) = 0; 
#endif
#ifdef ostreamiterator
  virtual  std::ostream_iterator< C_Data> enclosing_query( C_Window const &win, 
				           std::ostream_iterator< C_Data> out,
					   oit *dummy=0) = 0; 
#endif
  virtual bool is_inside( C_Window const &win,
			  C_Data const& object)=0;  
  virtual bool is_anchor()=0;
  virtual bool is_valid()=0;
};


// -------------------------------------------------------------------
// Tree Anchor: this class is used as a recursion anchor.
// The derived tree classes can be nested. Use this class as the
// most inner class. This class is doing nothin exept stopping the recursion

template <class C_Data, class C_Window>
class tree_anchor: public tree_base< C_Data,  C_Window>
{
public:
  // Construct a factory with the given factory as sublayer
  tree_anchor() {}
  virtual ~tree_anchor(){}
  tree_base<C_Data, C_Window> *clone() const { return new tree_anchor(); }
  typedef tree_base<C_Data, C_Window> tbt;
//  tree_base_type *clone() const { return new tree_anchor(); }

  bool make_tree(typename std::list< C_Data>::iterator& beg, 
		 typename std::list< C_Data>::iterator& end, 
		 typename tbt::lit *dummy=0) 
  {
    USE_ARGUMENT(beg);
    USE_ARGUMENT(end);
    return true;
  }
#ifdef stlvector
  bool make_tree(typename std::vector< C_Data>::iterator& beg, 
		 typename std::vector< C_Data>::iterator& end, 
		 typename tbt::vit *dummy=0) 
  {
    USE_ARGUMENT(beg);
    USE_ARGUMENT(end);
    return true;
  }
#endif
#ifdef carray
  bool make_tree( C_Data *beg, 
                  C_Data *end) 
  {
    USE_ARGUMENT(beg);
    USE_ARGUMENT(end);
    return true;
  }
#endif
   std::back_insert_iterator< std::list< C_Data> > 
      window_query( 
       C_Window const &win, 
       std::back_insert_iterator< std::list< C_Data> > out,
       typename tbt::lbit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
   
  std::back_insert_iterator< std::vector< C_Data> >  
      window_query( C_Window const &win, 
		    std::back_insert_iterator< std::vector< C_Data> > out, 
                    typename tbt::vbit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
#ifdef carray
   C_Data * window_query( C_Window const &win, 
                     C_Data * out){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
#endif
#ifdef ostreamiterator
   std::ostream_iterator< C_Data> window_query( C_Window const &win, 
				        std::ostream_iterator< C_Data> out, 
					typename tbt::oit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
#endif
   std::back_insert_iterator< std::list< C_Data> > enclosing_query( C_Window const &win, 
                                   std::back_insert_iterator< std::list< C_Data> > out,
				   typename tbt::lbit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
   std::back_insert_iterator< std::vector< C_Data> > enclosing_query( C_Window const &win, 
                                   std::back_insert_iterator< std::vector< C_Data> > out,
				   typename tbt::vbit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out); 
    return out;
  }
#ifdef carray
   C_Data * enclosing_query( C_Window const &win, 
                        C_Data * out){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
#endif
#ifdef ostreamiterator
   std::ostream_iterator< C_Data> enclosing_query( C_Window const &win, 
					   std::ostream_iterator< C_Data> out,
                                           typename tbt::oit *dummy=0){
    USE_ARGUMENT(win);
    USE_ARGUMENT(out);
    return out;
  }
#endif
  bool is_valid(){ return true;}

protected:

  bool is_inside( C_Window const &win, 
		  C_Data const& object)
  {     
    USE_ARGUMENT(win);
    USE_ARGUMENT(object);
    return true;
  }
  bool is_anchor(){return true;}
};

CGAL_END_NAMESPACE
// -------------------------------------------------------------------
#endif
