//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
// If not, contact Office of Research, University of Notre Dame, Notre
// Dame, IN 46556.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
#ifndef BOOST_GRAPH_DETAIL_IS_SAME_HPP
#define BOOST_GRAPH_DETAIL_IS_SAME_HPP

#include <boost/pending/ct_if.hpp>

namespace boost {
  struct false_tag;
  struct true_tag;

  namespace graph_detail {
    
#if !defined BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    template <class U, class V>
    struct is_same {
      typedef boost::false_tag is_same_tag; 
    };
    template <class U>
    struct is_same<U, U> {
      typedef boost::true_tag is_same_tag;
    };
#else
    template <class U, class V>
    struct is_same {
      enum { Unum = U::num, Vnum = V::num };
      typedef typename boost::ct_if< (Unum == Vnum),
               boost::true_tag, boost::false_tag>::type is_same_tag;
    };
#endif
  } // namespace graph_detail
} // namespace boost

#endif
