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

namespace CGAL {

/*! \class
 * Trapezoidal decomposition DAG node base class derived from Handle
 */
template<class T>
class Td_dag_node_base : public Handle
{
protected:
  void init() { PTR = 0; }

public:
  //c'tors
  
  Td_dag_node_base() {   init();   }

  Td_dag_node_base(const Td_dag_node_base<T>& x) : Handle(x) { }

  
  //operators overloading

  Td_dag_node_base& operator=(const Td_dag_node_base<T> & x) 
  {
    Handle::operator=(x); 
    return *this; 
  }

  bool operator!() const {  return PTR == 0;  }

};

/*! \class
 * Trapezoidal decomposition DAG node class derived from 
 * Td_dag_node_base class
 */
template<class T>
class Td_dag_node : public Td_dag_node_base<T>
{
public:
  //type of base class
  typedef Td_dag_node_base<T>   Td_dag_node_handle;

  //type of Td_dag_node (Self)
  typedef Td_dag_node<T>        Self;
  

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2

public:
  using Td_dag_node_handle::PTR;
  using Td_dag_node_handle::operator!;

#endif //CGAL_CFG_USING_BASE_MEMBER_BUG_2

protected:	
  /*! \class
   * Inner class node derived from Rep class
   */
  class node : public Rep
  {
    friend class Td_dag_node<T>;

  public:
    //c'tors
    node(const T& e,unsigned long depth=0) : 
      m_data(e),m_left_child(),m_right_child(),m_depth(depth)
    {}

    node(const T& e, const Td_dag_node_handle& left, 
         const Td_dag_node_handle& right,unsigned long depth=0) : 
      m_data(e),m_left_child(left),m_right_child(right),m_depth(depth)
    {}
    
    //d'tor
    ~node() {  }

    bool is_inner_node() const 
    {
      return !!m_left_child && !!m_right_child;
    }

    bool visited() const {  return m_visited; }

  protected:
    T m_data;			// information stored in node
    Td_dag_node_handle m_left_child, m_right_child; //left & right child nodes
    mutable unsigned long m_depth; 
    mutable bool m_visited;
  };

public:
  //c'tors

  Td_dag_node() { }
  
  Td_dag_node(const Td_dag_node_handle& dag):Td_dag_node_handle(dag) { }
  
  Td_dag_node(const Self& dag):Td_dag_node_handle(dag) { }

  Td_dag_node(const T& rootValue){  PTR = new node(rootValue);  }

  Td_dag_node(const T& rootValue, const Self& left, const Self& right)
  {
    PTR = new node(rootValue, left, right); 
    rebalance_depth();
  }

  //d'tor
  ~Td_dag_node() { }


  //information retrieval

  const Self& left_child() const
  {
    CGAL_precondition(!operator!());
    return *(const Self*)&ptr()->m_left_child;  
  }

  const Self& right_child() const
  {
    CGAL_precondition(!operator!());
    return *(const Self*)&ptr()->m_right_child;
  }

  T& get_data() const
  {
    CGAL_precondition(!operator!());
    return ptr()->m_data;
  }

  T& operator*() const
  {
    return get_data();
  }

  T* data_ptr() const
  {
    CGAL_precondition(!operator!());
    return &operator*();
  }

  T* operator->() const
  {
    return data_ptr();
  }

  bool is_inner_node() const 
  {
    return !operator!() && ptr()->is_inner_node();
  }

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

  //dynamic management:
  
  //Shallow copy
  Self& operator=(const Self& b)
  {
    Handle::operator=(b);
    return *this;
  }

  void set_data(const T& data)
  {
    if (!operator!()) 
      ptr()->m_data = data;
    else 
      operator=(Self(data));
  }

  void set_depth(unsigned long depth) const
  {
    ptr()->m_depth = depth;
  }

  void set_left_child(const Self& left)
  {
    CGAL_precondition(!operator!());
    ptr()->m_left_child = left;
    if (left.depth() < depth()+1) 
      left.ptr()->m_depth = depth()+1;
    left.rebalance_depth(); 
    // does nothing if left is a leaf
  }

  void set_right_child(const Self& right)
  {
    CGAL_precondition(!operator!());
    ptr()->m_right_child = right;
    if (right.depth() < depth()+1) 
      right.ptr()->m_depth = depth()+1;
    right.rebalance_depth(); 
    // does nothing if right is a leaf
  }

  void replace(const T& data,const Self& left,const Self& right)
  {
    set_data(data);
    set_left_child(left);
    set_right_child(right);
  }

  // Td_dag implementation not thread safe!
  void visit_none() const
  {
    if (!operator!())
    {
      ptr()->m_visited = false;
      left_child().visit_none();
      right_child().visit_none();
    }
  }

  void visit_one() const
  {
    if (!operator!())
      ptr()->m_visited = true;
  }

  /* -----output ---------------*/
#ifdef CGAL_PRE_IN_POST_ORDER

  void preorder() const
  {
    if (!operator!())
    {
      std::cout << operator*() << '\t';
      left_child().preorder();
      right_child().preorder();
    }
  }

  void inorder() const
  {
    if (!operator!())
    {
      left_child().inorder();
      std::cout << operator*() << '\t';
      right_child().inorder();
    }
  }

  void postorder() const
  {
    if (!operator!())
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
    visit_none();
    return recursive_filter(c,pr);
  }


protected:
  
  void rebalance_depth() const
  {
    if (is_inner_node()) 
    {  
      unsigned long depth_ = depth();
      if (left_child().depth() < depth_+1)
      {
        left_child().set_depth(depth_+1);
        left_child().rebalance_depth();
      }
      if (right_child().depth() < depth_+1)
      {
        right_child().set_depth(depth_+1);
        right_child().rebalance_depth();
      }
    }
  }
  
  unsigned long recursive_size() const
  {
    if (!operator!() && !ptr()->visited())
      return 1+ left_child().recursive_size() + right_child().recursive_size();
    return 0;
  }

  unsigned long recursive_depth() const
  {
    if (!operator!() && !ptr()->visited())
      return 1 + (std::max)(left_child().recursive_depth(), 
                            right_child().recursive_depth());
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
      left_child().recursive_filter(c, pr);
      right_child().recursive_filter(c, pr);
    }
    return c;
  }

private:
  
  node* ptr() const {   return (node*)PTR;  }

};


//io methods

template<class T,class Traits> 
std::ostream& write (std::ostream&  out, 
                     const Td_dag_node<T>& t,
                     const Traits& traits)
{
  static int depth;
  int i;
  if (!!t) 
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
  if (!!t) 
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
