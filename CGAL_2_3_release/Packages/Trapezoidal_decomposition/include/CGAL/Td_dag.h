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
// release_date  : 1999, October 13
//
// file          : include/CGAL/Td_dag.h
// package       : Trapezoidal decomposition 2
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

/* Directed acyclic binary graph template class */

#ifndef CGAL_TD_DAG_H
#define CGAL_TD_DAG_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <cstdlib>
#include <iostream>
#include <list>
#include <functional>
#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif
#ifndef CGAL_KERNEL_ASSERTIONS_H
#include <CGAL/kernel_assertions.h>
#endif
#ifndef CGAL_HANDLE_H
#include <CGAL/Handle.h>
#endif

CGAL_BEGIN_NAMESPACE

template<class T>
class Td_dag_base : public Leda_like_handle
{
public: //iddo (for CC-7.2) maybe protected?
  typedef T* pointer;
  typedef T& reference;
  typedef const T& const_reference;
protected:
  void init() {PTR=0;}
public:
  Td_dag_base() {init();}
  Td_dag_base(const Td_dag_base<T>& x) : Leda_like_handle(x) {}
  Td_dag_base& operator=(const Td_dag_base<T>& x) 
  {Leda_like_handle::operator=(x);return *this;}
  bool operator!() const {return PTR==0;}
};

template<class T>
class Td_dag : public Td_dag_base<T>
{
public:
  typedef T* pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef Td_dag_base<T> Td_dag_handle;
  typedef Td_dag<T> Self;
  typedef std::list<pointer> list_pointer;
protected:	
  class node : public Leda_like_rep
  {
#ifndef __BORLANDC__
    friend class Td_dag<T>;
#else
    typedef Td_dag<T> Td_dag;
    friend class Td_dag;
#endif
  public:
    node(const T& e,unsigned long _depth=0) : 
      data(e),leftPtr(),rightPtr(),depth_(_depth){}
    node(const T& e, const Td_dag_handle& left, 
         const Td_dag_handle& right,unsigned long _depth=0) : 
      data(e),leftPtr(left),rightPtr(right),depth_(_depth){}
    //		node(const T& e) : data(e),leftPtr(),rightPtr(){}
    //		node(const T& e, const Td_dag_handle& left, 
    // const Td_dag_handle& right) : data(e),leftPtr(left),rightPtr(right) {}
    ~node(){}
    bool is_inner_node() const {return !!leftPtr && !!rightPtr;}
    bool visited() const {return visited_;}
  protected:
    T data;			// information stored in node
    Td_dag_handle leftPtr,rightPtr;
    mutable unsigned long depth_;
    mutable bool visited_;
  };
  
public:
  /* -------constructors destructors -----*/
  Td_dag(){}
  Td_dag(const Td_dag_handle& dag):Td_dag_handle(dag){}
  Td_dag(const Self& dag):Td_dag_handle(dag){}
  Td_dag(const T& rootValue){PTR = new node(rootValue);}
  Td_dag(const T& rootValue, const Self& left, const Self& right)
  {PTR= new node(rootValue, left, right);rebalance_depth();}
  ~Td_dag(){}
  /* --------information retrieval -------*/
  const Self& left() const
  {
    CGAL_precondition(!operator!());
    return *(const Self*)&ptr()->leftPtr;
    
  }
  const Self& right() const
  {
    CGAL_precondition(!operator!());
    return *(const Self*)&ptr()->rightPtr;
  }
  reference operator*() const
  {
    CGAL_precondition(!operator!());
    return ptr()->data;
  }
  pointer operator->() const
  {
    CGAL_precondition(!operator!());
    return &operator*();
  }
  bool is_inner_node() const 
  {return !operator!() && ptr()->is_inner_node();}
  unsigned long size() const
  {
    visit_none();
    return recursive_size();
  }
  unsigned long depth() const
  {
    visit_none();
    return recursive_depth();
  }
  bool operator==(const Self& b) const
  {
    return PTR==b.PTR;
  }
  bool operator!=(const Self& b) const
  {
    return !operator==(b);
  }
  /* dynamic management ---------*/
  
  /* description:
     Shallow copy	*/
  Self& operator=(const Self& b)
  {
    Leda_like_handle::operator=(b);
    return *this;
  }
  /*
  Self& deep_copy(const Self& b)
  {
    if (this != &b)
      {
        clear();
        operator=(b);
      }
    return *this;
  }
  void clear()
  {
    if (!operator!())
      {
        left().clear();
        right().clear();
        operator=(Self());
      }
  }
  void detach_left()
  {
    if (!operator!())
      {
        // create dummy Td_dag
        T tmp;
        Self dummy(tmp);
        // detach left son,redirect to dummy
        set_left(dummy);
        // set left son pointer to 0
        ptr()->leftPtr.PTR=0;
        // delete dummy Td_dag
        delete dummy.ptr();
      }
  }
  void detach_right()
  {
    if (!operator!())
      {
				// create dummy Td_dag
        T tmp;
        Self dummy(tmp);
				// detach right son,redirect to dummy
        set_right(dummy);
				// set right son pointer to 0
        ptr()->rightPtr.PTR=0;
				// delete dummy Td_dag
        delete dummy.ptr();
      }
  }
  void detach()
  {
    detach_left();
    detach_right();
  }
  */
  void set_data(const T& data)
  {
    if (!operator!()) ptr()->data=data;
    else operator=(Self(data));
  }
  void set_depth(unsigned long _depth) const
  {
    ptr()->depth_=_depth;
  }
  void set_left(const Self& left)
  {
    CGAL_precondition(!operator!());
    ptr()->leftPtr=left;
    if (left.depth()<depth()+1) left.ptr()->depth_=depth()+1;
    left.rebalance_depth(); 
    // does nothing if right is a leaf
  }
  void set_right(const Self& right)
  {
    CGAL_precondition(!operator!());
    ptr()->rightPtr=right;
    if (right.depth()<depth()+1) right.ptr()->depth_=depth()+1;
    right.rebalance_depth(); 
    // does nothing if right is a leaf
  }
  void replace(const T& data,const Self& left,const Self& right)
  {
    set_data(data);
    set_left(left);
    set_right(right);
  }
  void visit_none() const
  {
    if (!operator!())
      {
        ptr()->visited_=false;
        left().visit_none();
        right().visit_none();
      }
  }
  void visit_one() const
  {
    if (!operator!())
      ptr()->visited_=true;
  }
  /* -----output ---------------*/
#ifdef CGAL_PRE_IN_POST_ORDER
  void preorder() const
  {
    if (!operator!())
      {
        std::cout << operator*() << '\t';
        left().preorder();
        right().preorder();
      }
  }
  void inorder() const
  {
    if (!operator!())
      {
        left().inorder();
        std::cout << operator*() << '\t';
        right().inorder();
      }
  }
  void postorder() const
  {
    if (!operator!())
      {
        left().postorder();
        right().postorder();
        std::cout << operator*() << '\t';
      }
  }
#endif
  
#if _MSC_VER>=1100
  friend std::ostream& operator<<(std::ostream&  out, const Self& t);
#endif
  
  template <class Container,class Predicate>
  Container& filter(Container& c,const Predicate& pr) const
  {
    visit_none();
    return recursive_filter(c,pr);
  }
protected:
  void rebalance_depth() const
  {
    if (is_inner_node()) 
      {  
        unsigned long depth_=depth();
        if (left().depth()<depth_+1)
          {
            left().set_depth(depth_+1);
            left().rebalance_depth();
          }
        if (right().depth()<depth_+1)
          {
            right().set_depth(depth_+1);
            right().rebalance_depth();
          }
      }
  }
  
  unsigned long recursive_size() const
  {
    if (!operator!() && !ptr()->visited())
      return 1+left().recursive_size()+right().recursive_size();
    else
      return 0;
  }
  unsigned long recursive_depth() const
  {
    if (!operator!() && !ptr()->visited())
      return 1+ std::max(left().recursive_depth(),right().recursive_depth());
    else
      return 0;
  }
  template <class Container,class Predicate>
  Container& recursive_filter(Container& c,const Predicate& pr) const
  {
    if (!operator!() && !ptr()->visited())
      {
        if (pr(operator*())) 
          c.insert(c.end(),operator*());
        visit_one();
        left().recursive_filter(c,pr);
        right().recursive_filter(c,pr);
      }
    return c;
  }
private:
  node* ptr() const {return (node*)PTR;}
};

template<class T,class Traits> 
std::ostream& write(
	std::ostream&  out, 
	const Td_dag<T>& t,
	const Traits& traits)
{
  static int depth;
  int i;
  if (!!t)
	{
		out << "\n";
		for(i=0;i<depth;i++)out << ">";
		out << "Data=";
		write(out,*t,traits);
		{
			depth++;
			out << "\n";
			for(i=0;i<depth;i++)out << ">";
			out << "left=";
			write(out,t.left(),traits);
			out << "\n";
			for(i=0;i<depth;i++)out << ">";
			out << "right=";
			write(out,t.right(),traits);
			depth--;
		}
	}
	else
	{
		out << "Empty";
	}
	return out ;
}

template<class T> std::ostream& operator<<(std::ostream&  out, 
                                           const Td_dag<T>& t)
{
  static int depth;
  int i;
  if (!!t)
	{
		out << "\n";
		for(i=0;i<depth;i++)out << ">";
		out << "Data=" << *t;
		{
			depth++;
			out << "\n";
			for(i=0;i<depth;i++)out << ">";
			out << "left=" << t.left();
			out << "\n";
			for(i=0;i<depth;i++)out << ">";
			out << "right=" <<t.right();
			depth--;
		}
	}
	else
	{
		out << "Empty";
	}
	return out ;
}

CGAL_END_NAMESPACE

#endif
/* 
   tech notes:
   The code is Handle designed.
   left(),right() are designed to cope with Handle(Handle& x) 
     precondition x.PTR!=0
   operator=() performs shallow copy
   operator*() returns data type
   output is done as a binary tree.
*/














