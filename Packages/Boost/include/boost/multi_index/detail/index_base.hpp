/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_INDEX_BASE_HPP
#define BOOST_MULTI_INDEX_DETAIL_INDEX_BASE_HPP

#include <boost/config.hpp> /* keep it first to prevent nasty warns in MSVC */
#include <boost/call_traits.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/multi_index/detail/copy_map.hpp>
#include <boost/multi_index/detail/node_type.hpp>
#include <boost/multi_index_container_fwd.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>

namespace boost{

namespace multi_index{

namespace detail{

/* The role of this class is threefold:
 *   - tops the linear hierarchy of indices.
 *   - terminates some cascading backbone function calls (insert_, etc.),
 *   - grants access to the backbone functions of the final
 *     multi_index_container class (for access restriction reasons, these
 *     cannot be called directly from the index classes.)
 */

template<typename Value,typename IndexSpecifierList,typename Allocator>
class index_base
{
protected:
  typedef index_node_base<Value>              node_type;
  typedef typename multi_index_node_type<
    Value,IndexSpecifierList,Allocator>::type final_node_type;
  typedef multi_index_container<
    Value,IndexSpecifierList,Allocator>       final_type;
  typedef tuples::null_type                   ctor_args_list;
  typedef typename 
    boost::detail::allocator::rebind_to<
      Allocator,
      typename Allocator::value_type>::type   final_allocator_type;
  typedef mpl::vector0<>                      index_type_list;
  typedef mpl::vector0<>                      iterator_type_list;
  typedef mpl::vector0<>                      const_iterator_type_list;
  typedef copy_map<
    final_node_type,
    final_allocator_type>                     copy_map_type;

private:
  typedef typename call_traits<Value>::param_type value_param_type;

protected:
  explicit index_base(const ctor_args_list&,const Allocator&){}

  void copy_(
    const index_base<Value,IndexSpecifierList,Allocator>&,const copy_map_type&)
  {}

  node_type* insert_(value_param_type v,node_type* x)
  {
    boost::detail::allocator::construct(&x->value,v);
    return x;
  }

  node_type* insert_(value_param_type v,node_type*,node_type* x)
  {
    boost::detail::allocator::construct(&x->value,v);
    return x;
  }

  void erase_(node_type* x)
  {
    boost::detail::allocator::destroy(&x->value);
  }

  void swap_(index_base<Value,IndexSpecifierList,Allocator>&){}

  bool replace_(value_param_type v,node_type* x)
  {
    x->value=v;
    return true;
  }

  bool modify_(node_type*){return true;}

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)
  /* invariant stuff */

  bool invariant_()const{return true;}
#endif

  /* access to backbone memfuns of Final class */

  final_type&       final(){return *static_cast<final_type*>(this);}
  const final_type& final()const{return *static_cast<const final_type*>(this);}

  final_node_type* final_header()const{return final().header();}

  bool        final_empty_()const{return final().empty_();}
  std::size_t final_size_()const{return final().size_();}
  std::size_t final_max_size_()const{return final().max_size_();}

  std::pair<final_node_type*,bool> final_insert_(value_param_type x)
    {return final().insert_(x);}
  std::pair<final_node_type*,bool> final_insert_(
    value_param_type x,final_node_type* position)
    {return final().insert_(x,position);}

  void final_erase_(final_node_type* x){final().erase_(x);}
  void final_swap_(final_type& x){final().swap_(x);}
  bool final_replace_(
    value_param_type k,final_node_type* x)
    {return final().replace_(k,x);}

  template<typename Modifier>
  bool final_modify_(Modifier mod,final_node_type* x)
    {return final().modify_(mod,x);}

#if defined(BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING)
  void final_check_invariant_()const{final().check_invariant_();}
#endif
};

} /* namespace multi_index::detail */

} /* namespace multi_index */

} /* namespace boost */

#endif
