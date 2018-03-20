// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DOUBLE_MAP_H
#define CGAL_DOUBLE_MAP_H

#include <CGAL/disable_warnings.h>

#include <map>
#include <utility>
#include <algorithm>
#include <CGAL/assertions.h>
#include <CGAL/function_objects.h> // for CGAL::Identity

#include <boost/version.hpp>
#if BOOST_VERSION >= 103500
#  define CGAL_USE_BOOST_BIMAP
#endif

#if defined(CGAL_USE_BOOST_BIMAP) && BOOST_VERSION == 104100
#include <CGAL/internal/container_fwd_fixed.hpp>
#endif

#ifdef CGAL_USE_BOOST_BIMAP
#include <CGAL/boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#endif

namespace CGAL {

template <class _Key, class _Data, class _Direct_compare = std::less<_Key>, 
          class _Reverse_compare = std::less<_Data> >
class Double_map
{
public:
  typedef _Key Key;
  typedef _Data Data;
  typedef _Direct_compare Direct_compare;
  typedef _Reverse_compare Reverse_compare;

  typedef Double_map<Key, Data, Direct_compare, Reverse_compare> Self;
  
#ifdef CGAL_USE_BOOST_BIMAP
  typedef ::boost::bimap< ::boost::bimaps::set_of<Key, Direct_compare>,
				 ::boost::bimaps::multiset_of<Data, Reverse_compare> > Boost_bimap;

  typedef typename Boost_bimap::left_map Direct_func;
  typedef typename Boost_bimap::right_map Reverse_func;
  typedef typename Reverse_func::iterator reverse_iterator;
  typedef typename Boost_bimap::relation relation;
  typedef typename Boost_bimap::left_key_type left_key_type;
  typedef typename Boost_bimap::left_value_type left_value_type;
  typedef typename Boost_bimap::right_key_type right_key_type;
  typedef typename Boost_bimap::right_value_type right_value_type;
#else
  typedef std::multimap <Data, Key, Reverse_compare> Reverse_func;
  typedef typename Reverse_func::iterator reverse_iterator;
  typedef std::map <Key, reverse_iterator, Direct_compare> Direct_func;
  typedef typename Self::Direct_func::key_type left_key_type;
  typedef typename Self::Direct_func::value_type left_value_type;
  typedef typename Self::Reverse_func::key_type right_key_type;
  typedef typename Self::Reverse_func::value_type right_value_type;
#endif
  typedef typename Direct_func::value_type Direct_entry;
               // std::pair<Key, reverse_iterator> 
  typedef typename Reverse_func::value_type Reverse_entry;
               // std::pair<Data, Key> ;

  typedef typename Direct_func::size_type size_type;


  typedef typename Direct_func::iterator direct_iterator;
  typedef reverse_iterator iterator;
  typedef typename Direct_func::const_iterator direct_const_iterator;

  typedef direct_const_iterator const_iterator;

private:
  // Private members
#ifdef CGAL_USE_BOOST_BIMAP
  Boost_bimap boost_bimap;

  const Direct_func& direct_func() const { return boost_bimap.left;}
  const Reverse_func& reverse_func() const { return boost_bimap.right;}
  Direct_func& direct_func() { return boost_bimap.left;}
  Reverse_func& reverse_func() { return boost_bimap.right;}
#else
  Direct_func direct_func_;
  Reverse_func reverse_func_;

  const Direct_func& direct_func() const { return direct_func_;}
  const Reverse_func& reverse_func() const { return reverse_func_;}
  Direct_func& direct_func() { return direct_func_;}
  Reverse_func& reverse_func() { return reverse_func_;}
#endif

public :
  // The default constructor
  Double_map () {}

  // The copy constructor
  Double_map(const Self& dm)
  {
    copy(dm);
  }

  // Queries
  bool empty() const
  {
    CGAL_assertion(is_valid());
    return(direct_func().empty());
  }

  size_type size() const
  {
    CGAL_assertion(is_valid());
    return direct_func().size();
  }

  bool is_valid() const
  {
    return(direct_func().size()==reverse_func().size());
  }
  
  void clear()
  {
#ifdef CGAL_USE_BOOST_BIMAP
    boost_bimap.clear();
#else
    direct_func().clear();
    reverse_func().clear();
#endif
    
  }

  Self& operator=(const Self& dm)
  {
    if(&dm != this)
      copy(dm);
    return *this;
  }

  void copy(const Self& dm)
  {
    clear();

#ifdef CGAL_USE_BOOST_BIMAP
    for(direct_const_iterator rit = dm.direct_func().begin();
	rit != dm.direct_func().end();
        ++rit)
    {
      direct_func().insert(*rit);
    }
#else
    reverse_func() = dm.reverse_func();
    
    for(reverse_iterator rit = reverse_func().begin();
        rit != reverse_func().end();
        ++rit)
    {
      // Fix an error with -D_GLIBCXX_DEBUG -D_GLIBCPP_CONCEPT_CHECKS
      // The following (commented) line triggered the error:
      //     direct_func()[rit->second] = rit;
      // The new following line fixes the bug. Actually, it is even more
      // efficient.
      direct_func().insert(std::make_pair(rit->second, rit));
    }
#endif
  }

  // Assignation
  bool insert(const Key& k, const Data& d)
    {
#ifdef CGAL_USE_BOOST_BIMAP
      direct_iterator hint = boost_bimap.left.lower_bound(k);

      if(hint != boost_bimap.left.end() && hint->first == k)
	return false;

      boost_bimap.left.insert(hint, Direct_entry(k, d));
      return true;
#else
      direct_iterator direct_hint = direct_func().lower_bound(k);

      if(direct_hint != direct_func().end() &&
	 direct_hint->first == k)
        return false;
      
      reverse_iterator reverse_it = reverse_func().insert(Reverse_entry(d, k));
      
      direct_func().insert(direct_hint, Direct_entry(k, reverse_it));

      CGAL_assertion(is_valid());
      return(true);
#endif
    }

  bool erase(const Key& k);

  // Access
  reverse_iterator front()
    {
      CGAL_assertion(is_valid() && !empty());
      return(reverse_func().begin());
    }

  void pop_front()
    {
#ifdef CGAL_USE_BOOST_BIMAP
      boost_bimap.right.erase(boost_bimap.right.begin());
#else
      CGAL_assertion(is_valid());
      reverse_iterator rit = reverse_func().begin();
      direct_iterator pos = direct_func().find(rit->second);
      CGAL_assertion(pos != direct_func().end());
      CGAL_assertion(pos->second == rit);
      
      direct_func().erase(pos);
      reverse_func().erase(rit);
      CGAL_assertion(is_valid());
#endif
    }

  const_iterator begin() const
  {
    return direct_func().begin();
  }
    
  const_iterator end() const
  {
    return direct_func().end();
  }

  template <typename Func_key, typename Func_data>
  void dump_direct_func(std::ostream& out,
			Func_key func_key,
			Func_data func_data)
  {
    for(typename Direct_func::iterator it = direct_func().begin();
	it != direct_func().end();
	++it)
    {
      out << func_key(it->first) << " -> "
#ifdef CGAL_USE_BOOST_BIMAP
	  << func_data(it->second)
#else
	  << func_data(it->second->first)
#endif
	  << std::endl;
    }
  }

  void dump_direct_func(std::ostream& out)
  {
    dump_direct_func(out, CGAL::Identity<Key>(), CGAL::Identity<Data>() );
  }

  template <typename Func_key, typename Func_data>
  void dump_reverse_func(std::ostream& out,
			 Func_key func_key,
			 Func_data func_data)
  {
    for(typename Reverse_func::iterator it = reverse_func().begin();
	it != reverse_func().end();
	++it)
    {
      out << func_data(it->first) << " "
	  << func_key(it->second) << std::endl;
    }
  }

  void dump_reverse_func(std::ostream& out) 
  {
    dump_reverse_func(out, CGAL::Identity<Data>(), CGAL::Identity<Key>());
  }
};

template <class _Key, class _Data, class _Direct_compare, 
  class _Reverse_compare>
bool
Double_map<_Key, _Data, _Direct_compare, _Reverse_compare>::
erase(const Key& k)
{
  CGAL_assertion(is_valid());
  direct_iterator pos = direct_func().find(k);
  if (pos == direct_func().end())
    return false;
  else
    {
#ifdef CGAL_USE_BOOST_BIMAP
      boost_bimap.left.erase(pos);
#else
      reverse_func().erase(pos->second);
      direct_func().erase(pos);
#endif
    }

  CGAL_assertion(is_valid());
  return true;
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_DOUBLE_MAP_H
