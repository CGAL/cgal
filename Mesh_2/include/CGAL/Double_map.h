// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
  
  typedef std::multimap <Data, Key, Reverse_compare> Reverse_func;
  typedef typename Reverse_func::iterator reverse_iterator;
  typedef std::map <Key, reverse_iterator, Direct_compare> Direct_func;

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
  Direct_func direct_func;
  Reverse_func reverse_func;

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
    return(direct_func.empty());
  }

  size_type size() const
  {
    CGAL_assertion(is_valid());
    return direct_func.size();
  }

  bool is_valid() const
  {
    return(direct_func.size()==reverse_func.size());
  }
  
  void clear()
  {
    direct_func.clear();
    reverse_func.clear();
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

    reverse_func = dm.reverse_func;
    
    for(reverse_iterator rit = reverse_func.begin();
        rit != reverse_func.end();
        ++rit)
    {
      // Fix an error with -D_GLIBCXX_DEBUG -D_GLIBCPP_CONCEPT_CHECKS
      // The following (commented) line triggered the error:
      //     direct_func[rit->second] = rit;
      // The new following line fixes the bug. Actually, it is even more
      // efficient.
      direct_func.insert(std::make_pair(rit->second, rit));
    }
  }

  // Assignation
  bool insert(const Key& k, const Data& d)
    {
      direct_iterator direct_hint = direct_func.lower_bound(k);

      if(direct_hint != direct_func.end() &&
	 direct_hint->first == k)
        return false;
      
      reverse_iterator reverse_it = reverse_func.insert(Reverse_entry(d, k));
      
      direct_func.insert(direct_hint, Direct_entry(k, reverse_it));

      CGAL_assertion(is_valid());
      return(true);
    }

  bool erase(const Key& k);

  // Access
  reverse_iterator front()
    {
      CGAL_assertion(is_valid() && !empty());
      return(reverse_func.begin());
    }

  void pop_front()
    {
      CGAL_assertion(is_valid());
      reverse_iterator rit = reverse_func.begin();
      direct_iterator pos = direct_func.find(rit->second);
      CGAL_assertion(pos != direct_func.end());
      CGAL_assertion(pos->second == rit);
      
      direct_func.erase(pos);
      reverse_func.erase(rit);
      CGAL_assertion(is_valid());
    }

  const_iterator begin() const
  {
    return direct_func.begin();
  }
    
  const_iterator end() const
  {
    return direct_func.end();
  }

  void dump_direct_func(std::ostream& out)
    {
      for(typename Direct_func::iterator it = direct_func.begin();
	  it != direct_func.end();
	  ++it)
	out << it->second->first << " " 
	    << "("
	    << it->first->vertex(0)->point() << ", "
	    << it->first->vertex(1)->point() << ", "
	    << it->first->vertex(2)->point()
	    << ")" << std::endl;
    }

  void dump_reverse_func(std::ostream& out)
    {
      for(typename Reverse_func::iterator it = reverse_func.begin();
	  it != reverse_func.end();
	  ++it)
        out << it->first << " " 
            << "("
            << it->second->vertex(0)->point() << ", "
            << it->second->vertex(1)->point() << ", "
            << it->second->vertex(2)->point()
	    << ")" << std::endl;
    }
};

template <class _Key, class _Data, class _Direct_compare, 
  class _Reverse_compare>
bool
Double_map<_Key, _Data, _Direct_compare, _Reverse_compare>::
erase(const Key& k)
{
  CGAL_assertion(is_valid());
  direct_iterator pos = direct_func.find(k);
  if (pos == direct_func.end())
    return false;
  else
    {
      reverse_func.erase(pos->second);
      direct_func.erase(pos);
    }

  CGAL_assertion(is_valid());
  return true;
}

} // end namespace CGAL

#endif // CGAL_DOUBLE_MAP_H
