// Copyright (c) 2001-2004  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_FILTRED_CONTAINER_H
#define CGAL_FILTRED_CONTAINER_H

CGAL_BEGIN_NAMESPACE

template <class Cont, class Pred>
class Filtered_container
{
  Cont cont;
  Pred test;
public:
  Filtered_container(Pred p=Pred()) : cont(), test(p) {};
  Filtered_container(Cont& c, Pred p=Pred()) : cont(c), test(p) {};

  typedef typename Cont::reference reference;
  typedef typename Cont::iterator iterator;
  typedef typename Cont::const_iterator const_iterator;
  typedef typename Cont::value_type value_type;

  inline
  reference front()
    {
      iterator r=cont.begin();
      while(!test(*r))
	{
	  cont.erase(r);
	  r=cont.begin();
	}
      return *r;
    }

  inline
  bool empty()
    {
      if(cont.empty())
	return true;
      else
	{
	  while(!cont.empty())
	    {
	      iterator r=cont.begin();
	    if(!test(*r))
	      pop_front();
	    else
	      return false;
	    }
	  return true;
	}
    }

  inline
  void pop_front() { cont.erase(cont.begin()); }

  inline
  void push_back(const value_type& e) { cont.push_back(e); }

  inline
  void insert(const value_type& e) { cont.insert(e); }

  inline
  void clear() { cont.clear(); }

};

CGAL_END_NAMESPACE

#endif
