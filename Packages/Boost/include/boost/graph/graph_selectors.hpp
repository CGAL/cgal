//=======================================================================
// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// This file is part of the Boost Graph Library
//
// You should have received a copy of the License Agreement for the
// Boost Graph Library along with the software; see the file LICENSE.
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

#ifndef BOOST_GRAPH_SELECTORS_HPP
#define BOOST_GRAPH_SELECTORS_HPP

namespace boost {

  //===========================================================================
  // Selectors for the Directed template parameter of adjacency_list
  // and adjacency_matrix.

  struct directedS { enum { is_directed = true, is_bidir = false }; 
    typedef true_type is_directed_t; 
    typedef false_type is_bidir_t;
  };
  struct undirectedS { 
    enum { is_directed = false, is_bidir = false }; 
    typedef false_type is_directed_t;
    typedef false_type is_bidir_t;
  };
  struct bidirectionalS { 
    enum { is_directed = true, is_bidir = true }; 
    typedef true_type is_directed_t;
    typedef true_type is_bidir_t;
  };

} // namespace boost

#endif // BOOST_GRAPH_SELECTORS_HPP
