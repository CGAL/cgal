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

#ifndef BOOST_GRAPH_EXCEPTION_HPP
#define BOOST_GRAPH_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace boost {

  struct bad_graph : public std::invalid_argument {
    bad_graph(const std::string& what_arg)
      : std::invalid_argument(what_arg) { }
  };

  struct not_a_dag : public bad_graph {
    not_a_dag()
        : bad_graph("The graph must be a DAG.") { } 
  };

  struct negative_edge : public bad_graph {
    negative_edge()
      : bad_graph("The graph may not contain an edge with negative weight."){ }
  };

  struct negative_cycle : public bad_graph {
    negative_cycle()
      : bad_graph("The graph may not contain negative cycles.") { }
  };
  struct not_connected : public bad_graph {
    not_connected()
      : bad_graph("The graph must be connected.") { }
  };

} // namespace boost

#endif // BOOST_GRAPH_EXCEPTION_HPP
