/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_SEQ_INDEX_NODE_HPP
#define BOOST_MULTI_INDEX_DETAIL_SEQ_INDEX_NODE_HPP

#include <algorithm>

namespace boost{

namespace multi_index{

namespace detail{

/* doubly-linked node for use by sequenced_index */

struct sequenced_index_node_impl
{
  sequenced_index_node_impl*&       prior(){return prior_;}
  sequenced_index_node_impl*const & prior()const{return prior_;}
  sequenced_index_node_impl*&       next(){return next_;}
  sequenced_index_node_impl*const & next()const{return next_;}

  /* interoperability with index_iterator */

  static void increment(sequenced_index_node_impl*& x){x=x->next();}
  static void decrement(sequenced_index_node_impl*& x){x=x->prior();}

  /* interoperability with index_proxy */

  static sequenced_index_node_impl* begin(sequenced_index_node_impl* header)
  {
    return header->next();
  }

  static sequenced_index_node_impl* end(sequenced_index_node_impl* header)
  {
    return header;
  }

  /* algorithmic stuff */

  static void link(
    sequenced_index_node_impl* x,sequenced_index_node_impl* header)
  {
    x->prior()=header->prior();
    x->next()=header;
    x->prior()->next()=x->next()->prior()=x;
  };

  static void unlink(sequenced_index_node_impl* x)
  {
    x->prior()->next()=x->next();
    x->next()->prior()=x->prior();
  }

  static void relink(
    sequenced_index_node_impl* position,sequenced_index_node_impl* x)
  {
    unlink(x);
    x->prior()=position->prior();
    x->next()=position;
    x->prior()->next()=x->next()->prior()=x;
  }

  static void relink(
    sequenced_index_node_impl* position,
    sequenced_index_node_impl* x,sequenced_index_node_impl* y)
  {
    /* position is assumed not to be in [x,y) */

    if(x!=y){
      sequenced_index_node_impl* z=y->prior();
      x->prior()->next()=y;
      y->prior()=x->prior();
      x->prior()=position->prior();
      z->next()=position;
      x->prior()->next()=x;
      z->next()->prior()=z;
    }
  }

  static void reverse(sequenced_index_node_impl* header)
  {
    sequenced_index_node_impl* x=header;
    do{
      sequenced_index_node_impl* y=x->next();
      std::swap(x->prior(),x->next());
      x=y;
    }while(x!=header);
  }

  static void swap(sequenced_index_node_impl* x,sequenced_index_node_impl* y)
  {
    /* This swap function does not exchange the header nodes,
     * but rather their pointers. This is *not* used for implementing
     * sequenced_index::swap.
     */

    if(x->next()!=x){
      if(y->next()!=y){
        std::swap(x->next(),y->next());
        std::swap(x->prior(),y->prior());
        x->next()->prior()=x->prior()->next()=x;
        y->next()->prior()=y->prior()->next()=y;
      }
      else{
        y->next()=x->next();
        y->prior()=x->prior();
        x->next()=x->prior()=x;
        y->next()->prior()=y->prior()->next()=y;
      }
    }
    else if(y->next()!=y){
      x->next()=y->next();
      x->prior()=y->prior();
      y->next()=y->prior()=y;
      x->next()->prior()=x->prior()->next()=x;
    }
  }

private:
  sequenced_index_node_impl();

  sequenced_index_node_impl* prior_;
  sequenced_index_node_impl* next_;
};

template<typename Super>
struct sequenced_index_node_trampoline:sequenced_index_node_impl{};

template<typename Super>
struct sequenced_index_node:Super,sequenced_index_node_trampoline<Super>
{
  sequenced_index_node_impl*&       prior(){return impl_type::prior();}
  sequenced_index_node_impl*const & prior()const{return impl_type::prior();}
  sequenced_index_node_impl*&       next(){return impl_type::next();}
  sequenced_index_node_impl*const & next()const{return impl_type::next();}

  sequenced_index_node_impl*       impl()
    {return static_cast<impl_type*>(this);}
  const sequenced_index_node_impl* impl()const
    {return static_cast<const impl_type*>(this);}

  static sequenced_index_node* from_impl(sequenced_index_node_impl *x)
    {return static_cast<sequenced_index_node*>(static_cast<impl_type*>(x));}
  static const sequenced_index_node* from_impl(
    const sequenced_index_node_impl* x)
  {
    return static_cast<const sequenced_index_node*>(
      static_cast<const impl_type*>(x));
  }

  /* interoperability with index_iterator */

  static void increment(sequenced_index_node*& x)
  {
    sequenced_index_node_impl* xi=x->impl();
    impl_type::increment(xi);
    x=from_impl(xi);
  }

  static void decrement(sequenced_index_node*& x)
  {
    sequenced_index_node_impl* xi=x->impl();
    impl_type::decrement(xi);
    x=from_impl(xi);
  }

  /* interoperability with index_proxy */

  static sequenced_index_node* begin(sequenced_index_node* header)
  {
    return from_impl(impl_type::begin(header->impl()));
  }

  static sequenced_index_node* end(sequenced_index_node* header)
  {
    return from_impl(impl_type::end(header->impl()));
  }

private:
  typedef sequenced_index_node_trampoline<Super> impl_type;
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
