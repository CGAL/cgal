// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_SPARSE_VECTOR_ITERATOR_H
#define CGAL_SPARSE_VECTOR_ITERATOR_H

#include <CGAL/basic.h>
#include <iterator>
#include <utility>
#include <vector>

CGAL_BEGIN_NAMESPACE

template < class NT >
class Sparse_vector_iterator {
  public:
    typedef  NT                                value_type;
    typedef  ptrdiff_t                         difference_type;
    typedef  value_type*                       pointer;
    typedef  value_type&                       reference;
    typedef  std::random_access_iterator_tag   iterator_category;

    typedef  Sparse_vector_iterator<NT>        Self;
    typedef  value_type                        Val;
    typedef  difference_type                   Dist;

    typedef  std::pair<Dist,Val>               Entry;
    typedef  std::vector<Entry>                Entries;
    typedef  typename Entries::const_iterator  EntryIterator;

    // forward operations
    Sparse_vector_iterator( ) : curr( 0), zero( 0) { }

    bool   operator == ( const Self& it) const { return ( curr == it.curr); }
    bool   operator != ( const Self& it) const { return ( curr != it.curr); }

    Val    operator *  ( ) const
	{
	    EntryIterator it;
	    for ( it = entries.begin(); it != entries.end(); ++it)
		if ( curr == it->first) return it->second;
	    return zero;
	}

    Self&  operator ++ (    ) {                   ++curr; return *this; }
    Self   operator ++ ( int) { Self tmp = *this; ++curr; return tmp;   }

    // bidirectional operations
    Self&  operator -- (    ) {                   --curr; return *this; }
    Self   operator -- ( int) { Self tmp = *this; --curr; return tmp;   }

    // random access operations
    Self&  operator += ( Dist n) { curr += n; return *this; }
    Self&  operator -= ( Dist n) { curr -= n; return *this; }

    Self   operator +  ( Dist n) const { Self tmp=*this; return tmp+=n; }
    Self   operator -  ( Dist n) const { Self tmp=*this; return tmp-=n; }

    Dist   operator - ( const Self& it) const { return curr - it.curr; }

    Val    operator [] ( Dist i) const
	{
	    Dist index = curr+i;
	    EntryIterator it;
	    for ( it = entries.begin(); it != entries.end(); ++it)
		if ( index == it->first) return it->second;
	    return zero;
	}

    bool   operator <  ( const Self& it) const { return ( curr <  it.curr); }
    bool   operator >  ( const Self& it) const { return ( curr >  it.curr); }
    bool   operator <= ( const Self& it) const { return ( curr <= it.curr); }
    bool   operator >= ( const Self& it) const { return ( curr >= it.curr); }

    // additional operations
    void  push_back( const Entry& entry) { entries.push_back( entry); }

  private:
    Dist                curr;
    NT                  zero;
    std::vector<Entry>  entries;
};

CGAL_END_NAMESPACE
  
#endif

// ===== EOF ==================================================================
