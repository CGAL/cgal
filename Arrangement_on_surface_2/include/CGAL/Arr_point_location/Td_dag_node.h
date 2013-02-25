// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Michal Balas <balasmic@post.tau.ac.il>
//                 (based on old version by Iddo Hanniel & Oren Nechushtan)

/* Directed acyclic binary graph template class */

#ifndef CGAL_TD_DAG_NODE_H
#define CGAL_TD_DAG_NODE_H

#include <CGAL/basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle.h>

#include <cstdlib>
#include <iostream>
#include <list>
#include <functional>

#include <boost/variant.hpp>

namespace CGAL {

/*! \class
 * Trapezoidal decomposition DAG node base class derived from Handle
 */
template<class Traits>
class Td_dag_node_base : public Handle
{
protected:
  void init() { PTR = 0; } //MICHAL: I think it is not used - so need to be removed

public:
  //c'tors
  
  Td_dag_node_base() {   init();   }

  Td_dag_node_base(const Td_dag_node_base<Traits>& x) : Handle(x) { }

  
  //operators overloading

  Td_dag_node_base& operator=(const Td_dag_node_base<Traits> & x) 
  {
    Handle::operator=(x); 
    return *this; 
  }

  //bool operator!() const {  return PTR == 0;  } //MICHAL: maybe use ptr(), and also can change to is_null or something similar
  bool is_null() const { return PTR == 0; }
protected:
  Rep * ptr() const { return (Rep*) PTR; }
  //Rep *& ptr() { return (Rep*) PTR; }
  void set_ptr(Rep* rep) { PTR = rep; } 

};




/*! \class
 * Trapezoidal decomposition DAG node class derived from 
 * Td_dag_node_base class
 */
template<class Traits>
class Td_dag_node : public Td_dag_node_base<Traits>
{

public:
  //type of base class
  typedef Td_dag_node_base<Traits>   Td_dag_node_handle;
  typedef Td_dag_node_base<Traits>   Base;

  //type of Td_dag_node (Self)
  typedef Td_dag_node<Traits>        Self;
  
  //type of td_map_item
  typedef typename Traits::Td_map_item Td_map_item;

  //type of Td_active_trapezoid
  typedef typename Traits::Td_active_trapezoid Td_active_trapezoid;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2

public:
  //using Td_dag_node_handle::PTR;
  //using Td_dag_node_handle::operator!;

#endif //CGAL_CFG_USING_BASE_MEMBER_BUG_2



protected:

  /*! \class
   * Inner class Node derived from Rep class
   */
  class Node : public Rep
  {
    friend class Td_dag_node<Traits>;

  public:
  
    class clear_neighbors_visitor : public boost::static_visitor< void  >
    {
    public:
      void operator()(Td_active_trapezoid& t) const
      {
        t.clear_neighbors();
      }

      template < typename Tp >
      void operator()(Tp& t) const
      {
      }

    };

    //c'tors
    Node(const Td_map_item& e,unsigned long depth=0) : 
      m_data(e), m_left_child(), m_right_child(), 
      m_depth(depth), m_visited(false)
    {}

    Node(const Td_map_item& e, const Td_dag_node_handle& left, 
         const Td_dag_node_handle& right, unsigned long depth=0) : 
      m_data(e), m_left_child(left), m_right_child(right), 
      m_depth(depth), m_visited(false)
    {}

    //d'tor
    ~Node() 
    {  
      boost::apply_visitor(clear_neighbors_visitor(), m_data);
    }

    bool is_inner_node() const //MICHAL: a node with only left child (like removed node) will be concidered as a leaf
    {
      //return !!m_left_child && !!m_right_child;
      return (!m_left_child.is_null() && !m_right_child.is_null());
    }

    bool visited() const {  return m_visited; }

  protected:
    //protected data members    
    
    //information stored in node
    Td_map_item m_data;			

    //left & right child nodes
    Td_dag_node_handle m_left_child;
    Td_dag_node_handle m_right_child; 
    
    //depth of the node in the DAG
    mutable unsigned long m_depth; 
    
    //visited flag for traversing the DAG
    mutable bool m_visited;
  };



public:
  //c'tors

  Td_dag_node() { }
  
  Td_dag_node(const Td_dag_node_handle& dag) : Td_dag_node_handle(dag) { }
  
  Td_dag_node(const Self& dag) : Td_dag_node_handle(dag) { }

  Td_dag_node(const Td_map_item& rootValue){  this->set_ptr(new Node(rootValue));  }

  Td_dag_node(const Td_map_item& rootValue, unsigned long depth)
  {  this->set_ptr(new Node(rootValue, depth));  }

  Td_dag_node(const Td_map_item& rootValue, const Self& left, const Self& right)
  {
    this->set_ptr(new Node( rootValue, left, right)); 
    depth_propagation();
  }

  Td_dag_node(const Td_map_item& rootValue, const Self& left, const Self& right,
              unsigned long depth)
  {
    this->set_ptr(new Node( rootValue, left, right, depth)); 
    depth_propagation();
  }

  //d'tor
  ~Td_dag_node() { }

  //information retrieval

  const Self& left_child() const
  {
    CGAL_precondition(!this->is_null());
    return *(const Self*)&this->node()->m_left_child;  
  }

  Self& left_child()
  {
    CGAL_precondition(!this->is_null());
    return (Self &)this->node()->m_left_child;  
  }

  const Self& right_child() const
  {
    CGAL_precondition(!this->is_null());
    return *(const Self*)&this->node()->m_right_child;
  }

  Self& right_child()
  {
    CGAL_precondition(!this->is_null());
    return (Self &)this->node()->m_right_child;
  }

  Td_map_item& get_data() const
  {
    CGAL_precondition(!this->is_null());
    return node()->m_data;
  }

  Td_map_item& operator*() const
  {
    return this->get_data();
  }

  Td_map_item* data_ptr() const
  {
    CGAL_precondition(!this->is_null());
    return &operator*();
  }

  Td_map_item* operator->() const
  {
    return this->data_ptr();
  }

  bool is_inner_node() const 
  {
    return !this->is_null() && this->node()->is_inner_node();
  }

  unsigned long size_inaccurate() const //exponential
  {
    init_visited();
    unsigned long res = recursive_size_inaccurate();
    init_visited();
    return res;
  }

  unsigned long size() const
  {
    init_visited();
    unsigned long res = recursive_size();
    init_visited();
    return res;
  }

  unsigned long rec_depth() const //exponential
  {
    init_visited();
    unsigned long res = recursive_depth();
    init_visited();
    return res;
  }
#if 0
  unsigned long rec_check(unsigned long rec_bound) const //exponential
  {
    unsigned long res = recursive_check(1,rec_bound);
    return res;
  }
#endif //0
  unsigned long max_depth() const
  {
    init_visited();
    unsigned long res = rec_max_depth();
    init_visited();
    return res;
  }

  const unsigned long& depth() const
  {
    return node()->m_depth;
  }

  unsigned long& depth()
  {
    return node()->m_depth;
  }

  bool operator==(const Self& b) const
  {
    return this->ptr() == b.ptr();
  }

  bool operator!=(const Self& b) const
  {
    return !operator==(b);
  }

  //dynamic management:
  
  //Shallow copy
  Self& operator=(const Self& b)
  {
    Handle::operator=(b);
    return *this;
  }

  void set_data(const Td_map_item& data)
  {
    if (!this->is_null()) 
      node()->m_data = data;
    else 
      operator=(Self(data));
  }

  void set_left_child(Self& left)
  {
    CGAL_precondition(!this->is_null());
    node()->m_left_child = left;
    if (left.depth() < depth()+1) 
      left.depth() = depth()+1;
    left.depth_propagation(); 
    // does nothing if left is a leaf
  }

  void set_right_child(Self& right)
  {
    CGAL_precondition(!this->is_null());
    node()->m_right_child = right;
    if (right.depth() < depth()+1) 
      right.depth() = depth()+1;
    right.depth_propagation(); 
    // does nothing if right is a leaf
  }

  void replace(const Td_map_item& data, Self& left, Self& right)
  {
    set_data(data);
    set_left_child(left);
    set_right_child(right);
  }

  // Td_dag implementation not thread safe!
  void init_visited() const
  {
    if (this->is_null() || node()->m_visited == false)
      return;
    node()->m_visited = false;
    left_child().init_visited();
    right_child().init_visited();
  }

  void visit_node() const
  {
    if (!this->is_null())
      node()->m_visited = true;
  }

  /* -----output ---------------*/
#ifdef CGAL_PRE_IN_POST_ORDER

  void preorder() const  //exponential
  {
    if (!is_null())
    {
      std::cout << operator*() << '\t';
      left_child().preorder();
      right_child().preorder();
    }
  }

  void inorder() const //exponential
  {
    if (!is_null())
    {
      left_child().inorder();
      std::cout << operator*() << '\t';
      right_child().inorder();
    }
  }

  void postorder() const //exponential
  {
    if (!is_null())
    {
      left_child().postorder();
      right_child().postorder();
      std::cout << operator*() << '\t';
    }
  }

#endif //CGAL_PRE_IN_POST_ORDER
  
  template <class Container,class Predicate>
  Container& filter(Container& c,const Predicate& pr) const
  {
    init_visited();
    Container& res = recursive_filter(c,pr);
    init_visited();
    return res;
  }


  

protected:
  
  //
  //Propagating depth for left child & right child if they exist
  //
  void depth_propagation() //exponential
  {
    if (!is_inner_node()) 
      return;
    
    if (left_child().depth() < depth() + 1)
    {
      left_child().depth() = depth() + 1;
      left_child().depth_propagation();
    }
    if (right_child().depth() < depth() + 1)
    {
      right_child().depth() = depth() + 1;
      right_child().depth_propagation();
    }
  }
  
  unsigned long recursive_depth() const
  {
    if (this->is_null() || node()->visited())
      return 0;
    return 1 + (std::max)(left_child().recursive_depth(), 
                            right_child().recursive_depth());
  }
#if 0
  unsigned long recursive_check(unsigned long curr_rec_depth, unsigned long rec_bound) const
  {
    if (is_null())
      return 0;
   //std::cout << curr_rec_depth << "," << std::flush;
    if ( curr_rec_depth > rec_bound + 30)
    {
      std::cout << "passed " << rec_bound + 30 << ", stopping\n";
      return 0;
    }
    return 1 + (std::max)(left_child().recursive_check(curr_rec_depth + 1, rec_bound), 
                          right_child().recursive_check(curr_rec_depth + 1, rec_bound));
  }
#endif //0

  unsigned long rec_max_depth() const
  {
    if (this->is_null() || node()->visited())
      return 0;
    visit_node();
    if (is_inner_node())
      return (std::max)(left_child().rec_max_depth(), 
                        right_child().rec_max_depth());
    else
      return depth();
  }

  unsigned long recursive_size_inaccurate() const
  {
    if (!this->is_null() && !node()->visited())
      return 1+ left_child().recursive_size_inaccurate() + right_child().recursive_size_inaccurate();
    return 0;
  }
  
  unsigned long recursive_size() const
  {
    if (this->is_null() || node()->visited())
      return 0;
    visit_node();
    return (1 + left_child().recursive_size() + right_child().recursive_size());  
  }

  template <class Container,class Predicate>
  Container& recursive_filter(Container& c,const Predicate& pr) const
  {
    if (this->is_null() || node()->visited())
      return c;
    if (pr(operator*())) 
      c.insert(c.end(),operator*());
    visit_node();
    left_child().recursive_filter(c, pr);
    right_child().recursive_filter(c, pr);
    return c;
  }

private:
  
  Node* node() const {   return (Node*)Base::PTR;  }

};

/*
//io methods

template<class T,class Traits> 
std::ostream& write (std::ostream&  out, 
                     const Td_dag_node<T>& t,
                     const Traits& traits)
{
  static int depth;
  int i;
  if (!t.is_null()) 
  {
    out << "\n";
    for(i=0; i<depth; i++)  out << ">";
    out << "Data=";
    write(out,*t,traits);
    {
      depth++;
      out << "\n";
      for(i=0; i<depth; i++)  out << ">";
      out << "left_child=";
      write(out,t.left_child(),traits);
      out << "\n";
      for(i=0; i<depth; i++)  out << ">";
      out << "right_child=";
      write(out,t.right_child(),traits);
      depth--;
    }
  }
  else
  {
    out << "Empty";
  }
  return out ;
}

template<class T> 
std::ostream& operator<< (std::ostream&  out, 
                          const Td_dag_node<T>& t)
{
  static int depth;
  int i;
  if (!t.is_null()) 
  {
    out << "\n";
    for(i=0; i<depth; i++)  out << ">";
    out << "Data=" << *t;
    {
      depth++;
      out << "\n";
      for(i=0; i<depth; i++)  out << ">";
      out << "left_child=" << t.left_child();
      out << "\n";
      for(i=0; i<depth; i++)  out << ">";
      out << "right_child=" <<t.right_child();
      depth--;
    }
  }
  else
  {
    out << "Empty";
  }
  return out ;
}
*/


} //namespace CGAL

#endif //CGAL_TD_DAG_NODE_H


/* 
   tech notes:
   The code is Handle designed.
   left_child(),right_child() are designed to cope with Handle(Handle& x) 
     precondition x.PTR!=0
   operator=() performs shallow copy
   operator*() returns data type
   output is done as a binary tree.
*/
