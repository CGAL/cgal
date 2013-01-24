// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DOUBLE_MAP_H
#define CGAL_DOUBLE_MAP_H

#include <map>
#include <utility>
#include <algorithm>
#include <CGAL/assertions.h>
#include <CGAL/function_objects.h> // for CGAL::Identity

#include <CGAL/internal/container_fwd_fixed.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

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
  Boost_bimap boost_bimap;

  const Direct_func& direct_func() const { return boost_bimap.left;}
  const Reverse_func& reverse_func() const { return boost_bimap.right;}
  Direct_func& direct_func() { return boost_bimap.left;}
  Reverse_func& reverse_func() { return boost_bimap.right;}

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
    boost_bimap.clear();   
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

    for(direct_const_iterator rit = dm.direct_func().begin();
      	rit != dm.direct_func().end();
        ++rit)
    {
      direct_func().insert(*rit);
    }
  }

  // Assignation
  bool insert(const Key& k, const Data& d)
    {
      direct_iterator hint = boost_bimap.left.lower_bound(k);

      if(hint != boost_bimap.left.end() && hint->first == k)
	return false;

      boost_bimap.left.insert(hint, Direct_entry(k, d));
      return true;
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
      boost_bimap.right.erase(boost_bimap.right.begin());
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
	  << func_data(it->second)
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
      boost_bimap.left.erase(pos);
    }

  CGAL_assertion(is_valid());
  return true;
}

} // end namespace CGAL

#endif // CGAL_DOUBLE_MAP_H
