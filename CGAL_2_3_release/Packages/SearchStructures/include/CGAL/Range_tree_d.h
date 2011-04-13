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
// file          : include/CGAL/Range_tree_d.h
// package       : SearchStructures (2.54)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// source        : include/CGAL/Range_tree_d.h
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


#ifndef CGAL_Range_tree_d__
#define CGAL_Range_tree_d__

#include <algorithm>
#include <iterator>
#include <functional>
#include <CGAL/Tree_base.h>
#include <list>
#include <vector>

// A d-dimensional Range Tree or a multilayer tree consisting of Range 
// and other trees that are derived public 
// Tree_base<C_Data, C_Window, C_Interface>
// can be construced within this class.
// C_Data: container class which contains the d-dimensional data the tree holds.
// C_Window: Query window -- a d-dimensional interval
// C_Interface: Interface for the class with functions that allow to 
// access the data.
// cf. file Tree_interface.h, class point_interface for the requirements.

CGAL_BEGIN_NAMESPACE

template <class C_Data, class C_Window, class C_Interface>
class Range_tree_d;

template <class C_Data, class C_Window, class C_Interface>
struct range_tree_node: public tree_node_base
{
  private	:
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef typename C_Interface::Key Key;
  typedef  C_Interface Interface;
  typedef typename tree_base< C_Data,  C_Window>::tree_base_type tree_base_type;
  protected:
  typedef Range_tree_d< C_Data,  C_Window,  C_Interface> rT_d;
public:
  friend rT_d;
  
  range_tree_node() 
  {
    sublayer = 0;// (tree_base_type *)0; 
  }
  
  range_tree_node( range_tree_node    * p_left,
		   range_tree_node    * p_right,
		   const  Data & v_obj,
		   const  Key  & v_key ) :
    object( v_obj ), key( v_key )
  {
    left_link = p_left;
    right_link = p_right;
    sublayer = 0;//(tree_base_type *)0; 
  }
  
  range_tree_node( range_tree_node    * p_left,
		   range_tree_node    * p_right,
		   const  Key  & v_key ) :
    key( v_key )
  {
    left_link = p_left;
    right_link = p_right;
    sublayer = 0;//(tree_base_type *)0; 
  }
  virtual ~range_tree_node()
  {
    if (sublayer != 0) //(tree_base_type *)
      delete sublayer;
  }
  
  Data object;
  Key key;
  tree_base_type *sublayer;
};

template <class C_Data, class C_Window, class C_Interface>
class Range_tree_d: public tree_base< C_Data,  C_Window>
{
 private:
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef typename C_Interface::Key Key;
  typedef  C_Interface Interface;
  typedef tree_base< C_Data,  C_Window>  tbt;
protected:
  typedef Range_tree_d< C_Data,  C_Window,  C_Interface> rT_d;
  tree_base<C_Data, C_Window> *sublayer_tree;
  //tree_base_type *sublayer_tree;
   C_Interface interface;
  int is_build;

 
  // A vertex is of this type:
  //  struct range_tree_node;

  friend range_tree_node<C_Data,C_Window,C_Interface>;

  typedef range_tree_node<C_Data,C_Window,C_Interface> range_tree_node2;
  typedef range_tree_node<C_Data,C_Window,C_Interface> *link_type;

  static link_type& left(link_type x) { 
  //  link_type left(link_type x) { 
    //    if ((*x).left_link!=0) 
      return CGAL__static_cast(link_type&, (*x).left_link);
    // return (*x).left_link;
      //return 0;
  }
  static link_type& right(link_type x) {
    // if ((*x).right_link!=0) 
      return CGAL__static_cast(link_type&, (*x).right_link);   
    //  return (*x).right_link;
      // return 0; 
  }

  static link_type& parent(link_type x) {
    //if ((*x).parent_link!=0) 
      return CGAL__static_cast(link_type&, (*x).parent_link);
    // return (*x).parent_link;
      //return 0;
  }

  link_type header;
  link_type node;
  link_type rightmost(){return right(header);}
  link_type leftmost(){return left(header);}
  link_type root(){
    if(header!=0) //TREE_BASE_NULL
      return CGAL__static_cast(link_type&, header->parent_link);
    // return parent(header);
    else 
      return 0; //TREE_BASE_NULL
  }

  bool is_less_equal(const Key  x, const Key  y)
  {
    return (!interface.comp(y,x));
  }  
  
  // this tree is not a recursion anchor
  bool is_anchor(){return false;}

  // returns true, if the object lies inside of win
  bool is_inside( C_Window const &win,  C_Data const& object)
  {
    if(is_less_equal(interface.get_left(win), interface.get_key(object)) 
       && interface.comp(interface.get_key(object),interface.get_right(win)))
   //half open
//       && is_less_equal(interface.get_key(object),interface.get_right(win)))
   //closed interval
    {
      return sublayer_tree->is_inside(win,object);
    }
    else
      return false;
  }


  // merge sort algorithms that takes O(n) time if the sequence to
  // be sorted consists of two sorted subsequences.
  template <class T>
  void dynamic_merge(T& first, T& last)
  {
    T prev, current=first;
    T current_first, current_middle, current_last;

    std::list<T> startpoints, tmp_startpoints;
    startpoints.push_back(current);
    prev = current++;

    while(current!=last)
    {
      if (interface.comp(interface.get_key(*current),interface.get_key(*prev)))
	startpoints.push_back(current);
      prev = current++;
    }
    while(startpoints.size()>1)
    {
      while(startpoints.size()>1)
      {
	current_first = startpoints.front();
	startpoints.erase(startpoints.begin());
	current_middle = startpoints.front();
	startpoints.erase(startpoints.begin());
	if(startpoints.size()>0)
	  current_last = startpoints.front();
	else 
	  current_last = last;
	tmp_startpoints.push_back(current_first);
	std::inplace_merge(current_first, current_middle, current_last, 
		      interface.key_comp);
      }
      if(startpoints.size()>0)
      {
	tmp_startpoints.push_back(startpoints.front());
	startpoints.erase(startpoints.begin());
      }
      startpoints.swap(tmp_startpoints);
    }
  }


  // recursive function 
  // (current,last) describe an interval of length n of sorted elements,
  // for this interval a tree is build containing these elements.
  // the most left child is returend in prevchild.

  template <class T>
  void build_range_tree(int n, link_type& leftchild, 
			link_type& rightchild,
			link_type& prevchild, 
			link_type& leftmostlink,
			T& current, 
			T& last,
			T& sublevel_first,
			T& sublevel_last)
  {
    // only two elements ==> two leaves and a parent is constructed
    if (n==2)
    {
      sublevel_first = current;

      link_type  vleft = new range_tree_node2( 0, 0,
                                  (*current), interface.get_key(*current) ); 
      //CGAL_NIL CGAL_NIL first two arguments
      CGAL_Tree_assertion( vleft != 0); //TREE_BASE_NULL

      ++current;
      link_type  vright = new range_tree_node2( 0,0,
                                  (*current), interface.get_key(*current) ); 
      //CGAL_NIL CGAL_NIL first two arguments
      CGAL_Tree_assertion( vright != 0); //TREE_BASE_NULL
      current++;
      sublevel_last = current;

      link_type  vparent = new range_tree_node2( vleft, vright, vleft->key );
      CGAL_Tree_assertion( vparent != 0); //TREE_BASE_NULL

      vleft->parent_link = vparent;
      vright->parent_link = vparent;
      leftchild = vleft;
      rightchild = vright;
      prevchild = vparent;
      if ( leftmostlink == 0) //TREE_BASE_NULL
	leftmostlink = leftchild;

      tree_base<C_Data, C_Window> *g = sublayer_tree->clone();
      
      T sub_first = sublevel_first;
      T sub_last = sublevel_last;
   
      g->make_tree(sub_first, sub_last);
      
      vparent->sublayer= g;
    }
    else
      // only one element ==> one leaf is constructed
      if(n==1)
      {
	sublevel_first = current;
	link_type vright = new range_tree_node2( 0, 0,
	                           (*current), interface.get_key(*current) );
	//CGAL_NIL CGAL_NIL first two arguments
        CGAL_Tree_assertion( vright != 0); //CGAL_NIL
	current++;
	sublevel_last = current;
	prevchild = vright;
	rightchild = vright;
      }
      else
      {
	// recursiv call for the construction. the interval is devided.
	T sublevel_left, sublevel_right;
	build_range_tree(n - (int)n/2, leftchild, rightchild, 
			 prevchild, leftmostlink, current, last, 
			 sublevel_first, sublevel_left);
	link_type vparent = new range_tree_node2( prevchild, 0,
                                        rightchild->key );
	//CGAL_NIL argument
        CGAL_Tree_assertion( vparent != 0); //TREE_BASE_NULL

	prevchild->parent_link = vparent;

	build_range_tree((int)n/2, leftchild, rightchild, 
			 prevchild, leftmostlink, current, 
			 last, sublevel_right, sublevel_last);
	vparent->right_link = prevchild;
	prevchild->parent_link = vparent;
	prevchild = vparent;
	tree_base<C_Data, C_Window> *g = sublayer_tree->clone();
	T sub_first = sublevel_first;
	T sub_last = sublevel_last;
	g->make_tree(sub_first, sub_last);
	vparent->sublayer = g;
      }
  }



  void delete_tree(link_type v)
  {
    if (v->left_link != 0) //TREE_BASE_NULL
    { 
       delete_tree(left(v));
       delete_tree(right(v));
    }
    delete v;
  }	    
		    
  
  // the vertex from that the way from root to the left interval bound 
  // and the right interval bound splits.
  link_type findSplitNode(Window const &key)
  {
    link_type v = root();

    while(v->left_link!=0) //TREE_BASE_NULL
    {
//      if(interface.comp(interface.get_right(key), v->key))
      if(is_less_equal(interface.get_right(key), v->key))
	v = left(v);
      else 
	if(interface.comp(v->key, interface.get_left(key)))
	  v = right(v);
	else
	  break;
    }

    return v;
  }

  template <class T>
  void report_subtree(link_type v, 
		      T result)
  {
    if(left(v)!=0) //TREE_BASE_NULL
    {
      report_subtree(left(v), result);
      report_subtree(right(v), result);
    }
    else
      (*result++)=v->object;
  }

  bool is_valid(link_type& v, link_type&  leftmost_child, 
		link_type& rightmost_child)
  {
    link_type leftmost_child_l, rightmost_child_l,  leftmost_child_r, 
      rightmost_child_r;
    if (v->sublayer != 0) //(tree_base_type *)
    {
      tree_base<C_Data, C_Window> *T= v->sublayer;
      if(!(*T).is_valid())
	return false;
    }
    if(left(v)!=0) //TREE_BASE_NULL
    {
      if(!is_valid(left(v), leftmost_child_l, rightmost_child_l))
	return false;
      if(!is_valid(right(v), leftmost_child_r, rightmost_child_r))
	return false;
      if(interface.comp((*v).key, (*rightmost_child_l).key) || 
	 interface.comp((*rightmost_child_l).key, (*v).key))
	return false;
      rightmost_child = rightmost_child_r;
      leftmost_child = leftmost_child_l;
    }
    else
    {
      rightmost_child = v;
      leftmost_child = v;      
    }
    return true;
  }




public:

  // construction of a tree
  Range_tree_d(Range_tree_d const &fact, bool):
    sublayer_tree(fact.sublayer_tree->clone()), is_build(false)
  {
    header = 0; //TREE_BASE_NULL
  }

  // construction of a tree
  Range_tree_d(tree_base<C_Data, C_Window> const &fact):
    sublayer_tree(fact.clone()), is_build(false)
  {
    header = 0; //TREE_BASE_NULL
  }

  // destruction
  virtual ~Range_tree_d()
  {
    link_type v=root();   

    if (v!=0)
      delete_tree(v);
      if (header!=0)
      	delete header;
      if (sublayer_tree!=0)
      	delete sublayer_tree;
  }


 // a prototype of the tree is returned
  tree_base<C_Data, C_Window> *clone() const 
  { 
    return new Range_tree_d(*this, true); 
  }
  
  bool make_tree(typename std::list< C_Data>::iterator& beg, 
		 typename std::list< C_Data>::iterator& end,
		 typename tbt::lit *dummy=0){ 
    return make_tree_impl(beg,end);
  }

#ifdef stlvector
  bool make_tree(typename std::vector< C_Data>::iterator& beg, 
		 typename std::vector< C_Data>::iterator& end,
		 typename tbt::vbit *dummy=0){ 
    return make_tree_impl(beg,end);
  }
#endif
#ifdef carray
  bool make_tree(C_Data *beg, 
		 C_Data *end){
    return make_tree_impl(beg,end);
  }
#endif

  // the tree is build according to the input elements in [first,last)
  template<class T>
  inline  
  bool make_tree_impl(T& first, 
		      T& last)
  {
    link_type leftchild, rightchild, prevchild, leftmostlink;

    if(!is_build)
      is_build = true;
    else
      return false;

    int n = count_elements__C( first, last );
    if(n==0)
    {
      is_build = false;
      return true;
    }

    dynamic_merge(first, last);
    
    leftmostlink = 0; //TREE_BASE_NULL
    T sublevel_first, sublevel_last;
    
    build_range_tree(n, leftchild, rightchild, prevchild, 
		     leftmostlink, first, last, 
		     sublevel_first, sublevel_last);
    
    header = new range_tree_node2();
    header->right_link = rightchild;
    header->parent_link = prevchild;
    header->left_link = leftmostlink;

    return true;
  }

  std::back_insert_iterator< std::list< C_Data> > window_query
          ( C_Window const &win, 
	    std::back_insert_iterator< std::list< C_Data> > out,
	    typename tbt::lbit *dummy=0){
    return window_query_impl(win,out);
  }


  std::back_insert_iterator< std::vector< C_Data> > window_query
          ( C_Window const &win, 
	    std::back_insert_iterator< std::vector< C_Data> > out,
	    typename tbt::vbit *dummy=0){
    return window_query_impl(win,out);
  }


#ifdef carray
  C_Data *window_query( C_Window const &win, C_Data *out){
    return window_query_impl(win,out);
  }
#endif

#ifdef ostreamiterator
  std::ostream_iterator< C_Data>  window_query( C_Window const &win, 
		     std::ostream_iterator< C_Data> out,
		     typename tbt::oit *dummy=0){
    return window_query_impl(win,out);
  }
#endif

  // all elements that ly in win are inserted in result
  template <class X>
  inline  
  X window_query_impl( C_Window const &win, X result)
  {
    if(is_less_equal(interface.get_right(win), interface.get_left(win)))
       return result;
    if(root()==0) //TREE_BASE_NULL
      return result;
    link_type split_node = findSplitNode(win);
    if(left(split_node)==0) //TREE_BASE_NULL
    {
      if(is_inside(win,split_node->object))
	(*result++)=split_node->object;
    }	  
    else
    {
      link_type v = (link_type)(*split_node).left_link;

      while(left(v)!=0) //TREE_BASE_NULL
      {
	if(is_less_equal(interface.get_left(win),v->key))
	{
	  link_type w = right(v);
	  if(left(w)!=0) //TREE_BASE_NULL
	  {
	    tree_base<C_Data, C_Window> *T= (w)->sublayer;
	    if(T->is_anchor())
	      report_subtree(w,result);
	    else
	      (*T).window_query(win, result);
	  }
	  else
	    if(is_inside(win,w->object))
	      (*result++)=(w)->object;
	  v = left(v);
	}
	else
	  v = right(v);
      }                 // end while
      if(is_inside(win,v->object))
	(*result++)=v->object;
      v = right(split_node);
      while(right(v)!=0) //TREE_BASE_NULL
      {
//	if(is_less_equal(v->key, interface.get_right(win))) closed interval
	if(interface.comp(v->key, interface.get_right(win))) 
	  //half open interval
	{
	  if(left(left(v))!=0) //TREE_BASE_NULL
	  {
	    tree_base<C_Data, C_Window> *T= (left(v))->sublayer;
	    if(T->is_anchor())
	      report_subtree(left(v),result);
	    else
	      (*T).window_query(win, result);
	  }
	  else
	  {
	    if(is_inside(win,left(v)->object))
	      (*result++)=left(v)->object; 
	  }
	  v = right(v);
	}
	else
	  v = left(v);
      }//end while
      if(is_inside(win,v->object))
      {
	(*result++)=v->object; 
      }
    }
    return result;
  }

  std::back_insert_iterator< std::list< C_Data> > enclosing_query( C_Window const &win, 
			     std::back_insert_iterator< std::list< C_Data> > out,
			     typename tbt::lbit *dummy=0){
    return enclosing_query_impl(win,out);
  }

  std::back_insert_iterator< std::vector< C_Data> > enclosing_query( C_Window const &win, 
			     std::back_insert_iterator< std::vector< C_Data> > out,
			     typename tbt::vbit *dummy=0){
    return enclosing_query_impl(win,out);
  }


#ifdef carray
  C_Data *enclosing_query( C_Window const &win, C_Data *out){
    return enclosing_query_impl(win,out);
  }
#endif

#ifdef ostreamiterator
  std::ostream_iterator< C_Data>  enclosing_query( C_Window const &win, 
			     std::ostream_iterator< C_Data> out,
			     typename tbt::oit *dummy=0){
    return enclosing_query_impl(win,out);
  }
#endif

  // a window query is performed 
  template <class T>
  inline
  T enclosing_query_impl(C_Window const &win, T result)
  {
    return window_query_impl(win, result);
  }

  bool is_valid()
  {
    link_type u,v,w;
    u=v=w=root();
    if(v!=0) //TREE_BASE_NULL
      return is_valid(u, v, w);
    return true;
  }
};

CGAL_END_NAMESPACE
#endif /* RANGE_TREE_H */




