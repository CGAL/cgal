/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 *
 * The internal implementation of red-black trees is based on that of SGI STL
 * stl_tree.h file: 
 *
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_ORD_INDEX_NODE_HPP
#define BOOST_MULTI_INDEX_DETAIL_ORD_INDEX_NODE_HPP

#include <algorithm>
#include <cstddef>

namespace boost{

namespace multi_index{

namespace detail{

/* definition of red-black nodes for ordered_index */

enum ordered_index_color{red=false,black=true};

struct ordered_index_node_impl
{
  ordered_index_color&            color(){return color_;}
  const ordered_index_color&      color()const{return color_;}
  ordered_index_node_impl*&       parent(){return parent_;}
  ordered_index_node_impl*const & parent()const{return parent_;}
  ordered_index_node_impl*&       left(){return left_;}
  ordered_index_node_impl*const & left()const{return left_;}
  ordered_index_node_impl*&       right(){return right_;}
  ordered_index_node_impl*const & right()const{return right_;}

  /* interoperability with index_iterator */

  static void increment(ordered_index_node_impl*& x)
  {
    if(x->right()){
      x=x->right();
      while(x->left())x=x->left();
    }
    else{
      ordered_index_node_impl* y=x->parent();
      while(x==y->right()){
        x=y;
        y=y->parent();
      }
      if(x->right()!=y)x=y;
    }
  }

  static void decrement(ordered_index_node_impl*& x)
  {
    if(x->color()==red&&x->parent()->parent()==x){
      x=x->right();
    }
    else if(x->left()){
      ordered_index_node_impl* y=x->left();
      while(y->right())y=y->right();
      x=y;
    }else{
      ordered_index_node_impl* y=x->parent();
      while(x==y->left()){
        x=y;
        y=y->parent();
      }
      x=y;
    }
  }

  /* interoperability with index_proxy */

  static ordered_index_node_impl* begin(ordered_index_node_impl* header)
  {
    return header->left();
  }

  static ordered_index_node_impl* end(ordered_index_node_impl* header)
  { 
    return header;
  }

  /* algorithmic stuff */

  static void rotate_left(
    ordered_index_node_impl* x,ordered_index_node_impl*& root)
  {
    ordered_index_node_impl* y=x->right();
    x->right()=y->left();
    if(y->left())y->left()->parent()=x;
    y->parent()=x->parent();
    
    if(x==root)                    root=y;
    else if(x==x->parent()->left())x->parent()->left()=y;
    else                           x->parent()->right()=y;
    y->left()=x;
    x->parent()=y;
  }

  static ordered_index_node_impl* minimum(ordered_index_node_impl* x)
  {
    while(x->left())x=x->left();
    return x;
  }

  static ordered_index_node_impl* maximum(ordered_index_node_impl* x)
  {
    while(x->right())x=x->right();
    return x;
  }

  static void rotate_right(
    ordered_index_node_impl* x,ordered_index_node_impl*& root)
  {
    ordered_index_node_impl* y=x->left();
    x->left()=y->right();
    if(y->right())y->right()->parent()=x;
    y->parent()=x->parent();

    if(x==root)                     root=y;
    else if(x==x->parent()->right())x->parent()->right()=y;
    else                            x->parent()->left()=y;
    y->right()=x;
    x->parent()=y;
  }

  static void rebalance(
    ordered_index_node_impl* x,ordered_index_node_impl*& root)
  {
    x->color()=red;
    while(x!=root&&x->parent()->color()==red){
      if(x->parent()==x->parent()->parent()->left()){
        ordered_index_node_impl* y=x->parent()->parent()->right();
        if(y&&y->color()==red){
          x->parent()->color()=black;
          y->color()=black;
          x->parent()->parent()->color()=red;
          x=x->parent()->parent();
        }
        else{
          if(x==x->parent()->right()){
            x=x->parent();
            rotate_left(x,root);
          }
          x->parent()->color()=black;
          x->parent()->parent()->color()=red;
          rotate_right(x->parent()->parent(),root);
        }
      }
      else{
        ordered_index_node_impl* y=x->parent()->parent()->left();
        if(y&&y->color()==red){
          x->parent()->color()=black;
          y->color()=black;
          x->parent()->parent()->color()=red;
          x=x->parent()->parent();
        }
        else{
          if(x==x->parent()->left()){
            x=x->parent();
            rotate_right(x,root);
          }
          x->parent()->color()=black;
          x->parent()->parent()->color()=red;
          rotate_left(x->parent()->parent(),root);
        }
      }
    }
    root->color()=black;
  }

  static ordered_index_node_impl* rebalance_for_erase(
    ordered_index_node_impl* z,ordered_index_node_impl*& root,
    ordered_index_node_impl*& leftmost,ordered_index_node_impl*& rightmost)
  {
    ordered_index_node_impl* y=z;
    ordered_index_node_impl* x=0;
    ordered_index_node_impl* x_parent=0;
    if(y->left()==0){        /* z has at most one non-null child. y==z. */
      x=y->right();          /* x might be null */
    }
    else{
      if(y->right()==0)  {     /* z has exactly one non-null child. y==z. */
        x=y->left();           /* x is not null */
      }
      else{                    /* z has two non-null children.  Set y to */
        y=y->right();          /*   z's successor.  x might be null.     */
        while(y->left())y=y->left();
        x=y->right();
      }
    }
    if(y!=z){
      z->left()->parent()=y;   /* relink y in place of z. y is z's successor */
      y->left()=z->left();
      if(y!=z->right()){
        x_parent=y->parent();
        if(x) x->parent()=y->parent();
        y->parent()->left()=x; /* y must be a child of left */
        y->right()=z->right();
        z->right()->parent()=y;
      }
      else{
        x_parent=y;
      }

      if(root==z)                    root=y;
      else if(z->parent()->left()==z)z->parent()->left()=y;
      else                           z->parent()->right()=y;
      y->parent()=z->parent();
      std::swap(y->color(),z->color());
      y=z;                    /* y now points to node to be actually deleted */
    }
    else{                     /* y==z */
      x_parent=y->parent();
      if(x)x->parent()=y->parent();   
      if(root==z){
        root=x;
      }
      else{
        if(z->parent()->left()==z)z->parent()->left()=x;
        else                      z->parent()->right()=x;
      }
      if(leftmost==z){
        if(z->right()==0){      /* z->left() must be null also */
          leftmost=z->parent();
        }
        else{              
          leftmost=minimum(x);  /* makes leftmost==header if z==root */
        }
      }
      if(rightmost==z){
        if(z->left()==0){       /* z->right() must be null also */
          rightmost=z->parent();
        }
        else{                   /* x==z->left() */
          rightmost=maximum(x); /* makes rightmost==header if z==root */
        }
      }
    }
    if(y->color()!=red){
      while(x!=root&&(x==0 || x->color()==black)){
        if(x==x_parent->left()){
          ordered_index_node_impl* w=x_parent->right();
          if(w->color()==red){
            w->color()=black;
            x_parent->color()=red;
            rotate_left(x_parent,root);
            w=x_parent->right();
          }
          if((w->left()==0||w->left()->color()==black) &&
             (w->right()==0||w->right()->color()==black)){
            w->color()=red;
            x=x_parent;
            x_parent=x_parent->parent();
          } 
          else{
            if(w->right()==0 
                || w->right()->color()==black){
              if(w->left()) w->left()->color()=black;
              w->color()=red;
              rotate_right(w,root);
              w=x_parent->right();
            }
            w->color()=x_parent->color();
            x_parent->color()=black;
            if(w->right())w->right()->color()=black;
            rotate_left(x_parent,root);
            break;
          }
        } 
        else{                   /* same as above,with right <-> left */
          ordered_index_node_impl* w=x_parent->left();
          if(w->color()==red){
            w->color()=black;
            x_parent->color()=red;
            rotate_right(x_parent,root);
            w=x_parent->left();
          }
          if((w->right()==0||w->right()->color()==black) &&
             (w->left()==0||w->left()->color()==black)){
            w->color()=red;
            x=x_parent;
            x_parent=x_parent->parent();
          }
          else{
            if(w->left()==0||w->left()->color()==black){
              if(w->right())w->right()->color()=black;
              w->color()=red;
              rotate_left(w,root);
              w=x_parent->left();
            }
            w->color()=x_parent->color();
            x_parent->color()=black;
            if(w->left())w->left()->color()=black;
            rotate_right(x_parent,root);
            break;
          }
        }
      }
      if(x)x->color()=black;
    }
    return y;
  }

  static void restore(
    ordered_index_node_impl* x,ordered_index_node_impl* prior,
    ordered_index_node_impl* next,ordered_index_node_impl* header)
  {
    if(next->left()==0){
      next->left()=x;
      x->parent()=next;
      if(next==header->left()){
        header->left()=x;         /* maintain leftmost pointing to min node */
      }
    }
    else if(x!=prior){            /* prior->right() must be null */
      prior->right()=x;
      x->parent()=prior;
      if(prior==header->right()){
        header->right()=x;        /* maintain rightmost pointing to max node */
      }
    }
    else{
      header->parent()=x;
      header->left()=x;
      header->right()=x;
      x->parent()=header;
    }
    x->left()=0;
    x->right()=0;
    rebalance(x,header->parent());
  }

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)
  /* invariant stuff */

  static std::size_t black_count(
    ordered_index_node_impl* node,ordered_index_node_impl* root)
  {
    if(!node)return 0;
    std::size_t sum=0;
    for(;;){
      if(node->color()==black)++sum;
      if(node==root)break;
      node=node->parent();
    } 
    return sum;
  }
#endif

private:
  ordered_index_node_impl();

  ordered_index_color      color_; 
  ordered_index_node_impl* parent_;
  ordered_index_node_impl* left_;
  ordered_index_node_impl* right_;
};

template<typename Super>
struct ordered_index_node_trampoline:ordered_index_node_impl{};

template<typename Super>
struct ordered_index_node:Super,ordered_index_node_trampoline<Super>
{
  ordered_index_color&            color(){return impl_type::color();}
  const ordered_index_color&      color()const{return impl_type::color();}
  ordered_index_node_impl*&       parent(){return impl_type::parent();}
  ordered_index_node_impl*const & parent()const{return impl_type::parent();}
  ordered_index_node_impl*&       left(){return impl_type::left();}
  ordered_index_node_impl*const & left()const{return impl_type::left();}
  ordered_index_node_impl*&       right(){return impl_type::right();}
  ordered_index_node_impl*const & right()const{return impl_type::right();}

  ordered_index_node_impl* impl(){return static_cast<impl_type*>(this);}
  const ordered_index_node_impl* impl()const
    {return static_cast<const impl_type*>(this);}

  static ordered_index_node* from_impl(ordered_index_node_impl *x)
  {
    return static_cast<ordered_index_node*>(static_cast<impl_type*>(x));
  }
  
  static const ordered_index_node* from_impl(const ordered_index_node_impl* x)
  {
    return static_cast<const ordered_index_node*>(
      static_cast<const impl_type*>(x));
  }

  /* interoperability with index_iterator */

  static void increment(ordered_index_node*& x)
  {
    ordered_index_node_impl* xi=x->impl();
    impl_type::increment(xi);
    x=from_impl(xi);
  }

  static void decrement(ordered_index_node*& x)
  {
    ordered_index_node_impl* xi=x->impl();
    impl_type::decrement(xi);
    x=from_impl(xi);
  }

  /* interoperability with index_proxy */

  static ordered_index_node* begin(ordered_index_node* header)
  {
    return from_impl(impl_type::begin(header->impl()));
  }

  static ordered_index_node* end(ordered_index_node* header)
  {
    return from_impl(impl_type::end(header->impl()));
  }

private:
  typedef ordered_index_node_trampoline<Super> impl_type;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
