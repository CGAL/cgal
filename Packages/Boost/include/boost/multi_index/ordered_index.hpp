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

#ifndef BOOST_MULTI_INDEX_ORDERED_INDEX_HPP
#define BOOST_MULTI_INDEX_ORDERED_INDEX_HPP

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <algorithm>
#include <boost/call_traits.hpp>
#include <boost/detail/no_exceptions_support.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/multi_index/detail/access_specifier.hpp>
#include <boost/multi_index/detail/index_iterator.hpp>
#include <boost/multi_index/detail/modify_key_adaptor.hpp>
#include <boost/multi_index/detail/ord_index_node.hpp>
#include <boost/multi_index/detail/ord_index_ops.hpp>
#include <boost/multi_index/detail/prevent_eti.hpp>
#include <boost/multi_index/detail/safe_mode.hpp>
#include <boost/multi_index/detail/scope_guard.hpp>
#include <boost/multi_index/detail/unbounded.hpp>
#include <boost/multi_index/detail/value_compare.hpp>
#include <boost/multi_index/ordered_index_fwd.hpp>
#include <boost/ref.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)
#define BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT                          \
  detail::scope_guard BOOST_JOIN(check_invariant_,__LINE__)=                 \
    detail::make_obj_guard(*this,&ordered_index::check_invariant_);          \
  BOOST_JOIN(check_invariant_,__LINE__).touch();
#else
#define BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT
#endif

namespace boost{

namespace multi_index{

namespace detail{

/* ordered_index adds a layer of indexing to a given Super */

/* Most of the implementation of unique and non-unique indices is
 * shared. We tell from one another on instantiation time by using
 * these tags.
 */

struct ordered_unique_tag{};
struct ordered_non_unique_tag{};

template<
  typename KeyFromValue,typename Compare,
  typename Super,typename TagList,typename Category
>

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
class ordered_index:
  BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS Super,
  public index_proxy<ordered_index_node<typename Super::node_type> >
#else
class ordered_index:
  BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS Super,
  public safe_container<
    ordered_index<KeyFromValue,Compare,Super,TagList,Category> >
#endif
#else
class ordered_index:
  BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS Super
#endif

{ 
#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)&&\
    BOOST_WORKAROUND(__MWERKS__,<=0x3003)
/* The "ISO C++ Template Parser" option in CW8.3 has a problem with the
 * lifetime of const references bound to temporaries --precisely what
 * scopeguards are.
 */

#pragma parse_mfunc_templ off
#endif

protected:
  typedef ordered_index_node<
      typename Super::node_type>                     node_type;

public:
  /* types */

  typedef typename KeyFromValue::result_type         key_type;
  typedef typename node_type::value_type             value_type;
  typedef KeyFromValue                               key_from_value;
  typedef Compare                                    key_compare;
  typedef value_comparison<
    value_type,KeyFromValue,Compare>                 value_compare;
  typedef tuple<key_from_value,key_compare>          ctor_args;
  typedef typename Super::final_allocator_type       allocator_type;
  typedef typename allocator_type::reference         reference;
  typedef typename allocator_type::const_reference   const_reference;

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
  typedef index_iterator<node_type>                  iterator;
  typedef index_iterator<node_type>                  const_iterator;
#else
  typedef index_iterator<node_type,ordered_index>    iterator;
  typedef index_iterator<node_type,ordered_index>    const_iterator;
#endif
#else
  typedef index_iterator<node_type>                  iterator;
  typedef index_iterator<node_type>                  const_iterator;
#endif

  typedef std::size_t                                size_type;      
  typedef std::ptrdiff_t                             difference_type;
  typedef typename allocator_type::pointer           pointer;
  typedef typename allocator_type::const_pointer     const_pointer;
  typedef typename
    boost::reverse_iterator<iterator>                reverse_iterator;
  typedef typename
    boost::reverse_iterator<const_iterator>          const_reverse_iterator;
  typedef typename TagList::type                     tag_list;

protected:
  typedef typename Super::final_node_type            final_node_type;
  typedef tuples::cons<
    ctor_args, 
    typename Super::ctor_args_list>                  ctor_args_list;
  typedef typename mpl::push_front<
    typename Super::index_type_list,
    ordered_index>::type                             index_type_list;
  typedef typename mpl::push_front<
    typename Super::iterator_type_list,
    iterator>::type    iterator_type_list;
  typedef typename mpl::push_front<
    typename Super::const_iterator_type_list,
    const_iterator>::type                            const_iterator_type_list;
  typedef typename Super::copy_map_type              copy_map_type;

private:
#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
#if BOOST_WORKAROUND(BOOST_MSVC,<1300)
  typedef index_proxy<
      ordered_index_node<
        typename Super::node_type> >                 safe_super;
#else
  typedef safe_container<ordered_index>              safe_super;
#endif
#endif

  typedef typename call_traits<
    value_type>::param_type                          value_param_type;
  typedef typename call_traits<
    key_type>::param_type                            key_param_type;

public:

  /* construct/copy/destroy
   * Default and copy ctors are in the protected section as indices are
   * not supposed to be created on their own. No range ctor either.
   */

  ordered_index<KeyFromValue,Compare,Super,TagList,Category>& operator=(
    const ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x)
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    this->final()=x.final();
    return *this;
  }

  allocator_type get_allocator()const
  {
    return this->final().get_allocator();
  }

  /* iterators */

  iterator               begin(){return make_iterator(leftmost());}
  const_iterator         begin()const{return make_iterator(leftmost());}
  iterator               end(){return make_iterator(header());}
  const_iterator         end()const{return make_iterator(header());}
  reverse_iterator       rbegin(){return make_reverse_iterator(end());}
  const_reverse_iterator rbegin()const{return make_reverse_iterator(end());}
  reverse_iterator       rend(){return make_reverse_iterator(begin());}
  const_reverse_iterator rend()const{return make_reverse_iterator(begin());}
 
  /* capacity */

  bool      empty()const{return this->final_empty_();}
  size_type size()const{return this->final_size_();}
  size_type max_size()const{return this->final_max_size_();}

  /* modifiers */

  std::pair<iterator,bool> insert(value_param_type x)
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    std::pair<final_node_type*,bool> p=this->final_insert_(x);
    return std::pair<iterator,bool>(make_iterator(p.first),p.second);
  }

  iterator insert(iterator position,value_param_type x)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(position,*this);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    std::pair<final_node_type*,bool> p=this->final_insert_(
      x,static_cast<final_node_type*>(position.get_node()));
    return make_iterator(p.first);
  }
    
  template<typename InputIterator>
  void insert(InputIterator first,InputIterator last)
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    iterator hint=end();
    for(;first!=last;++first)hint=insert(hint,*first);
  }

  void erase(iterator position)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(position,*this);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
    /* MSVC++ 6.0 optimizer on safe mode code chokes if this
     * this is not added. Left it for all compilers as it does no
     * harm.
     */

    position.detach();
#endif

    this->final_erase_(static_cast<final_node_type*>(position.get_node()));
  }
  
  size_type erase(key_param_type x)
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    std::pair<iterator,iterator> p=equal_range(x);
    size_type s=0;
    while(p.first!=p.second){
      erase(p.first++);
      ++s;
    }
    return s;
  }

  void erase(iterator first,iterator last)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(first);
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(last);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(first,*this);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(last,*this);
    BOOST_MULTI_INDEX_CHECK_VALID_RANGE(first,last);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    while(first!=last){
      erase(first++);
    }
  }

  bool replace(iterator position,value_param_type x)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(position,*this);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    return this->final_replace_(
      x,static_cast<final_node_type*>(position.get_node()));
  }

  template<typename Modifier>
  bool modify(iterator position,Modifier mod)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(position,*this);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
    /* MSVC++ 6.0 optimizer on safe mode code chokes if this
     * this is not added. Left it for all compilers as it does no
     * harm.
     */

    position.detach();
#endif

    return this->final_modify_(
      mod,static_cast<final_node_type*>(position.get_node()));
  }

  template<typename Modifier>
  bool modify_key(iterator position,Modifier mod)
  {
    BOOST_MULTI_INDEX_CHECK_VALID_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_DEREFERENCEABLE_ITERATOR(position);
    BOOST_MULTI_INDEX_CHECK_IS_OWNER(position,*this);
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    return modify(
      position,modify_key_adaptor<Modifier,value_type,KeyFromValue>(mod,key));
  }

  void swap(ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x)
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    this->final_swap_(x.final());
  }

  void clear()
  {
    BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT;
    erase(begin(),end());
  }

  /* observers */

  key_from_value key_extractor()const{return key;}
  key_compare    key_comp()const{return comp;}
  value_compare  value_comp()const{return value_compare(key,comp);}

  /* set operations */

  /* no need to provide non-const versions as
   * ordered_index::iterator==ordered_index::const_iterator
   */

  template<typename CompatibleKey>
  const_iterator find(const CompatibleKey& x)const
  {
    return make_iterator(ordered_index_find(header(),key,x,comp));
  }

  template<typename CompatibleKey,typename CompatibleCompare>
  const_iterator find(
    const CompatibleKey& x,const CompatibleCompare& comp)const
  {
    return make_iterator(ordered_index_find(header(),key,x,comp));
  }

  template<typename CompatibleKey>
  size_type count(const CompatibleKey& x)const
  {
    return count(x,comp);
  }

  template<typename CompatibleKey,typename CompatibleCompare>
  size_type count(const CompatibleKey& x,const CompatibleCompare& comp)const
  {
    std::pair<const_iterator, const_iterator> p=equal_range(x,comp);
    size_type n=std::distance(p.first,p.second);
    return n;
  }

  template<typename CompatibleKey>
  const_iterator lower_bound(const CompatibleKey& x)const
  {
    return make_iterator(ordered_index_lower_bound(header(),key,x,comp));
  }

  template<typename CompatibleKey,typename CompatibleCompare>
  const_iterator lower_bound(
    const CompatibleKey& x,const CompatibleCompare& comp)const
  {
    return make_iterator(ordered_index_lower_bound(header(),key,x,comp));
  }

  template<typename CompatibleKey>
  const_iterator upper_bound(const CompatibleKey& x)const
  {
    return make_iterator(ordered_index_upper_bound(header(),key,x,comp));
  }

  template<typename CompatibleKey,typename CompatibleCompare>
  const_iterator upper_bound(
    const CompatibleKey& x,const CompatibleCompare& comp)const
  {
    return make_iterator(ordered_index_upper_bound(header(),key,x,comp));
  }

  template<typename CompatibleKey>
  std::pair<const_iterator,const_iterator> equal_range(
    const CompatibleKey& x)const
  {
    return equal_range(x,comp);
  }

  template<typename CompatibleKey,typename CompatibleCompare>
  std::pair<const_iterator,const_iterator> equal_range(
    const CompatibleKey& x,const CompatibleCompare& comp)const
  {
    return std::pair<const_iterator,const_iterator>(
      lower_bound(x,comp),upper_bound(x,comp));
  }

  /* range */

  /* no need to provide non-const version as
   * ordered_index::iterator==ordered_index::const_iterator
   */

  template<typename LowerBounder,typename UpperBounder>
  std::pair<const_iterator,const_iterator>
  range(LowerBounder lower,UpperBounder upper)const
  {
    std::pair<const_iterator,const_iterator> p(
      lower_range(lower),upper_range(upper));
    if(p.second!=end()&&(p.first==end()||comp(key(*p.second),key(*p.first)))){
      p.second=p.first;
    }
    return p;
  }

BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS:
  ordered_index(const ctor_args_list& args_list,const allocator_type& al):
    Super(args_list.get_tail(),al),

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)&&\
    BOOST_WORKAROUND(BOOST_MSVC,<1300)
    safe_super(final_header()),
#endif

    key(tuples::get<0>(args_list.get_head())),
    comp(tuples::get<1>(args_list.get_head()))
  {
    empty_initialize();
  }

  ordered_index(
    const ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x):
    Super(x),

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)&&\
    BOOST_WORKAROUND(BOOST_MSVC,<1300)
    safe_super(final_header()),
#endif

    key(x.key),
    comp(x.comp)
  {
    /* Copy ctor just takes the key and compare objects from x. The rest is
     * done in subsequent call to copy_().
     */
  }

  ~ordered_index()
  {
    /* the container is guaranteed to be empty by now */
  }

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
  iterator       make_iterator(node_type* node){return iterator(node,this);}
  const_iterator make_iterator(node_type* node)const
    {return const_iterator(node,const_cast<ordered_index*>(this));}
#else
  iterator       make_iterator(node_type* node){return iterator(node);}
  const_iterator make_iterator(node_type* node)const
                   {return const_iterator(node);}
#endif

  void copy_(
    const ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x,
    const copy_map_type& map)
  {
    if(!x.root()){
      empty_initialize();
    }
    else{
      header()->color()=x.header()->color();

      node_type* root_cpy=map.find(static_cast<final_node_type*>(x.root()));
      header()->parent()=root_cpy->impl();

      node_type* leftmost_cpy=map.find(
        static_cast<final_node_type*>(x.leftmost()));
      header()->left()=leftmost_cpy->impl();

      node_type* rightmost_cpy=map.find(
        static_cast<final_node_type*>(x.rightmost()));
      header()->right()=rightmost_cpy->impl();

      typedef typename copy_map_type::const_iterator copy_map_iterator;
      for(copy_map_iterator it=map.begin(),it_end=map.end();it!=it_end;++it){
        node_type* org=it->first;
        node_type* cpy=it->second;

        cpy->color()=org->color();

        ordered_index_node_impl* parent_org=org->parent();
        if(!parent_org)cpy->parent()=0;
        else{
          node_type* parent_cpy=map.find(
            static_cast<final_node_type*>(node_type::from_impl(parent_org)));
          cpy->parent()=parent_cpy->impl();
          if(parent_org->left()==org->impl()){
            parent_cpy->left()=cpy->impl();
          }
          else if(parent_org->right()==org->impl()){
            /* header() does not satisfy this nor the previous check */
            parent_cpy->right()=cpy->impl();
          }
        }

        if(!org->left())cpy->left()=0;
        if(!org->right())cpy->right()=0;
      }
    }
    
    Super::copy_(x,map);
  }

  node_type* insert_(value_param_type v,node_type* x)
  {
    node_type* res=link2(key(v),x,Category());
    if(res!=x)return res;
    else{
      BOOST_TRY{
        res=static_cast<node_type*>(Super::insert_(v,x));
        if(res!=x){
          ordered_index_node_impl::rebalance_for_erase(
            x->impl(),header()->parent(),header()->left(),header()->right());
        }
        return res;
      }
      BOOST_CATCH(...){
        ordered_index_node_impl::rebalance_for_erase(
          x->impl(),header()->parent(),header()->left(),header()->right());
        BOOST_RETHROW;
      }
      BOOST_CATCH_END
    }
  }

  node_type* insert_(value_param_type v,node_type* position,node_type* x)
  {
    node_type* res=link3(key(v),position,x,Category());
    if(res!=x)return res;
    else{
      BOOST_TRY{
        res=static_cast<node_type*>(Super::insert_(v,position,x));
        if(res!=x){
          ordered_index_node_impl::rebalance_for_erase(
            x->impl(),header()->parent(),header()->left(),header()->right());
        }
        return res;
      }
      BOOST_CATCH(...){
        ordered_index_node_impl::rebalance_for_erase(
          x->impl(),header()->parent(),header()->left(),header()->right());
        BOOST_RETHROW;
      }
      BOOST_CATCH_END
    }
  }

  void erase_(node_type* x)
  {
    Super::erase_(x);
    ordered_index_node_impl::rebalance_for_erase(
      x->impl(),header()->parent(),header()->left(),header()->right());

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
    detach_iterators(x);
#endif
  }

  void swap_(ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x)
  {
    std::swap(key,x.key);
    std::swap(comp,x.comp);

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
    safe_super::swap(x);
#endif

    Super::swap_(x);
  }

  bool replace_(value_param_type v,node_type* x)
  {
    if(!comp(key(v),key(x->value))&&
       !comp(key(x->value),key(v))){
      return Super::replace_(v,x);
    }

    node_type* prior=x;
    node_type::decrement(prior);
    node_type* next=x;
    node_type::increment(next);

    ordered_index_node_impl::rebalance_for_erase(
      x->impl(),header()->parent(),header()->left(),header()->right());

    BOOST_TRY{
      if(link2(key(v),x,Category())!=x){
        ordered_index_node_impl::restore(
          x->impl(),prior->impl(),next->impl(),header()->impl());
        return false;
      }

      BOOST_TRY{
        if(!Super::replace_(v,x)){
          ordered_index_node_impl::rebalance_for_erase(
            x->impl(),header()->parent(),header()->left(),header()->right());
          ordered_index_node_impl::restore(
            x->impl(),prior->impl(),next->impl(),header()->impl());
          return false;
        }
        else return true;
      }
      BOOST_CATCH(...){
        ordered_index_node_impl::rebalance_for_erase(
          x->impl(),header()->parent(),header()->left(),header()->right());
        BOOST_RETHROW;
      }
      BOOST_CATCH_END
    }
    BOOST_CATCH(...){
      ordered_index_node_impl::restore(
        x->impl(),prior->impl(),next->impl(),header()->impl());
      BOOST_RETHROW;
    }
    BOOST_CATCH_END
  }

  bool modify_(node_type* x)
  {
    bool b;
    BOOST_TRY{
      b=in_place(x,Category());
    }
    BOOST_CATCH(...){
      erase_(x);
      BOOST_RETHROW;
    }
    BOOST_CATCH_END
    if(!b){
      ordered_index_node_impl::rebalance_for_erase(
        x->impl(),header()->parent(),header()->left(),header()->right());
      BOOST_TRY{
        if(link2(key(x->value),x,Category())!=x){
          Super::erase_(x);

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
          detach_iterators(x);
#endif
          return false;
        }
      }
      BOOST_CATCH(...){
        Super::erase_(x);

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
        detach_iterators(x);
#endif

        BOOST_RETHROW;
      }
      BOOST_CATCH_END
    }

    BOOST_TRY{
      if(!Super::modify_(x)){
        ordered_index_node_impl::rebalance_for_erase(
          x->impl(),header()->parent(),header()->left(),header()->right());

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
        detach_iterators(x);
#endif

        return false;
      }
      else return true;
    }
    BOOST_CATCH(...){
      ordered_index_node_impl::rebalance_for_erase(
        x->impl(),header()->parent(),header()->left(),header()->right());

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
      detach_iterators(x);
#endif

      BOOST_RETHROW;
    }
    BOOST_CATCH_END
  }

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)
  /* invariant stuff */

  bool invariant_()const
  {
    if(size()==0||begin()==end()){
      if(size()!=0||begin()!=end()||
         header()->left()!=header()->impl()||
         header()->right()!=header()->impl())return false;
    }
    else{
      if((size_type)std::distance(begin(),end())!=size())return false;

      std::size_t len=ordered_index_node_impl::black_count(
        leftmost()->impl(),root()->impl());
      for(const_iterator it=begin(),it_end=end();it!=it_end;++it){
        node_type* x=it.get_node();
        node_type* left_x=node_type::from_impl(x->left());
        node_type* right_x=node_type::from_impl(x->right());

        if(x->color()==red){
          if((left_x&&left_x->color()==red)||
             (right_x&&right_x->color()==red))return false;
        }
        if(left_x&&comp(key(x->value),key(left_x->value)))return false;
        if(right_x&&comp(key(right_x->value),key(x->value)))return false;
        if(!left_x&&!right_x&&
           ordered_index_node_impl::black_count(
             x->impl(),root()->impl())!=len)
          return false;
      }
    
      if(leftmost()->impl()!=
         ordered_index_node_impl::minimum(root()->impl()))
        return false;
      if(rightmost()->impl()!=
         ordered_index_node_impl::maximum(root()->impl()))
        return false;
    }

    return Super::invariant_();
  }

  
  /* This forwarding function eases things for the boost::mem_fn construct
   * in BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT. Actually,
   * final_check_invariant is already an inherited member function of
   * ordered_index.
   */
  void check_invariant_()const{this->final_check_invariant_();}
#endif

private:
  node_type* header()const{return this->final_header();}
  node_type* root()const{return node_type::from_impl(header()->parent());}
  node_type* leftmost()const{return node_type::from_impl(header()->left());}
  node_type* rightmost()const{return node_type::from_impl(header()->right());}

  void empty_initialize()
  {
    header()->color()=red;
    /* used to distinguish header() from root, in iterator.operator++ */
    
    header()->parent()=0;
    header()->left()=header()->impl();
    header()->right()=header()->impl();
  }

  node_type* link4(key_param_type k,node_type* x,node_type* y,node_type* z)
  {
    if(x!=0||y==header()||comp(k,key(y->value))){
      y->left()=z->impl(); /* also makes leftmost()=z when y==header() */
      if (y==header()){
        header()->parent()=z->impl();
        header()->right()=z->impl();
      }
      else if(y==leftmost()){
        header()->left()=z->impl();
        /* maintain leftmost() pointing to min node */
      }
    }
    else{
      y->right()=z->impl();
      if(y==rightmost()){
        header()->right()=z->impl();
        /* maintain rightmost() pointing to max node */
      }
    }
    z->parent()=y->impl();
    z->left()=0;
    z->right()=0;
    ordered_index_node_impl::rebalance(z->impl(),header()->parent());
    return z;
  }

  node_type* link2(key_param_type k,node_type* z,ordered_unique_tag)
  {
    node_type* y=header();
    node_type* x=root();
    bool c=true;
    while(x){
      y=x;
      c=comp(k,key(x->value));
      x=node_type::from_impl(c?x->left():x->right());
    }
    iterator j=make_iterator(y);   
    if(c){
      if(j==begin())return link4(k,x,y,z);
      else --j;
    }

    if(comp(key(*j),k))return link4(k,x,y,z);
    else return j.get_node();
  }

  node_type* link2(key_param_type k,node_type* z,ordered_non_unique_tag)
  {
    node_type* y=header();
    node_type* x=root();
    while (x){
     y=x;
     x=node_type::from_impl(comp(k,key(x->value))?x->left():x->right());
    }
    return link4(k,x,y,z);
  }

  node_type* link3(
    key_param_type k,node_type* position,node_type* z,ordered_unique_tag)
  {
    if(position->impl()==header()->left()){ 
      if(size()>0&&comp(k,key(position->value))){
        return link4(k,position,position,z);
      }
      else return link2(k,z,ordered_unique_tag());
    } 
    else if(position==header()){ 
      if(comp(key(rightmost()->value),k)){
        return link4(k,0,rightmost(),z);
      }
      else return link2(k,z,ordered_unique_tag());
    } 
    else{
      node_type* before=position;
      node_type::decrement(before);
      if(comp(key(before->value),k)&&comp(k,key(position->value))){
        if(before->right()==0)return link4(k,0,before,z); 
        else return link4(k,position,position,z);
      } 
      else return link2(k,z,ordered_unique_tag());
    }
  }

  node_type* link3(
    key_param_type k,node_type* position,node_type* z,ordered_non_unique_tag)
  {
    if(position->impl()==header()->left()){ 
      if(size()>0&&!comp(key(position->value),k)){
        return link4(k,position,position,z);
      }
      else return link2(k,z,ordered_non_unique_tag());
    } 
    else if(position==header()){
      if(!comp(k,key(rightmost()->value))){
        return link4(k,0,rightmost(),z);
      }
      else return link2(k,z,ordered_non_unique_tag());
    } 
    else{
      node_type* before=position;
      node_type::decrement(before);
      if (!comp(k,key(before->value))&&!comp(key(position->value),k)){
        if(before->right()==0)return link4(k,0,before,z); 
        else return link4(k,position,position,z);
      } 
      else return link2(k,z,ordered_non_unique_tag());
    }
  }

  bool in_place(node_type* x,ordered_unique_tag)
  {
    node_type* y=x;
    node_type::decrement(y);
    if(y!=header()&&!comp(key(y->value),key(x->value)))return false;

    y=x;
    node_type::increment(y);
    return y==header()||comp(key(x->value),key(y->value));
  }

  bool in_place(node_type* x,ordered_non_unique_tag)
  {
    node_type* y=x;
    node_type::decrement(y);
    if(y!=header()&&comp(key(x->value),key(y->value)))return false;

    y=x;
    node_type::increment(y);
    return y==header()||!comp(key(y->value),key(x->value));
  }

#if defined(BOOST_MULTI_INDEX_ENABLE_SAFE_MODE)
  void detach_iterators(node_type* x)
  {
    iterator it=make_iterator(x);
    safe_mode::detach_equivalent_iterators(it);
  }
#endif

  template<typename LowerBounder>
  const_iterator lower_range(LowerBounder lower)const
  {
    node_type* y=header();
    node_type* z=root();

    while(z){
      if(lower(key(z->value))){
        y=z;
        z=node_type::from_impl(z->left());
      }
      else z=node_type::from_impl(z->right());
    }

    return make_iterator(y);
  }

  const_iterator lower_range(unbounded_type)const
  {
    return begin();
  }

  template<typename UpperBounder>
  const_iterator upper_range(UpperBounder upper)const
  {
    node_type* y=header();
    node_type* z=root();

    while(z){
      if(!upper(key(z->value))){
        y=z;
        z=node_type::from_impl(z->left());
      }
      else z=node_type::from_impl(z->right());
    }

    return make_iterator(y);
  }

  const_iterator upper_range(unbounded_type)const
  {
    return end();
  }

  key_from_value key;
  key_compare    comp;

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)&&\
    BOOST_WORKAROUND(__MWERKS__,<=0x3003)
#pragma parse_mfunc_templ reset
#endif
};

/* comparison */

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator==(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return x.size()==y.size()&&std::equal(x.begin(),x.end(),y.begin());
}

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator<(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return std::lexicographical_compare(x.begin(),x.end(),y.begin(),y.end());
}

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator!=(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return !(x==y);
}

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator>(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return y<x;
}

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator>=(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return !(x<y);
}

template<
  typename KeyFromValue1,typename Compare1,
  typename Super1,typename TagList1,typename Category1,
  typename KeyFromValue2,typename Compare2,
  typename Super2,typename TagList2,typename Category2
>
bool operator<=(
  const ordered_index<KeyFromValue1,Compare1,Super1,TagList1,Category1>& x,
  const ordered_index<KeyFromValue2,Compare2,Super2,TagList2,Category2>& y)
{
  return !(x>y);
}

/*  specialized algorithms */

template<
  typename KeyFromValue,typename Compare,
  typename Super,typename TagList,typename Category
>
void swap(
  ordered_index<KeyFromValue,Compare,Super,TagList,Category>& x,
  ordered_index<KeyFromValue,Compare,Super,TagList,Category>& y)
{
  x.swap(y);
}

} /* namespace multi_index::detail */

/* ordered_index specifiers */

template<typename Arg1,typename Arg2,typename Arg3>
struct ordered_unique
{
  typedef typename detail::ordered_index_args<
    Arg1,Arg2,Arg3>                                index_args;
  typedef typename index_args::tag_list_type       tag_list_type;
  typedef typename index_args::key_from_value_type key_from_value_type;
  typedef typename index_args::compare_type        compare_type;

  template<typename Super>
  struct node_class
  {
    typedef detail::ordered_index_node<Super> type;
  };

  template<typename Super>
  struct index_class
  {
    typedef detail::ordered_index<
      key_from_value_type,compare_type,
      Super,tag_list_type,detail::ordered_unique_tag> type;
  };
};

template<typename Arg1,typename Arg2,typename Arg3>
struct ordered_non_unique
{
  typedef detail::ordered_index_args<
    Arg1,Arg2,Arg3>                                index_args;
  typedef typename index_args::tag_list_type       tag_list_type;
  typedef typename index_args::key_from_value_type key_from_value_type;
  typedef typename index_args::compare_type        compare_type;

  template<typename Super>
  struct node_class
  {
    typedef detail::ordered_index_node<Super> type;
  };

  template<typename Super>
  struct index_class
  {
    typedef detail::ordered_index<
      key_from_value_type,compare_type,
      Super,tag_list_type,detail::ordered_non_unique_tag> type;
  };
};

} /* namespace multi_index */

} /* namespace boost */

#undef BOOST_MULTI_INDEX_ORD_INDEX_CHECK_INVARIANT

#endif
