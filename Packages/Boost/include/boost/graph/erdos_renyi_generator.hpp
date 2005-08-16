// Copyright 2004 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine
#ifndef BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP
#define BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP

#include <iterator>
#include <utility>
#include <boost/random/uniform_int.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost {

  template<typename RandomGenerator, typename Graph>
  class erdos_renyi_iterator
  {
    typedef typename graph_traits<Graph>::directed_category directed_category;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;

    BOOST_STATIC_CONSTANT
      (bool,
       is_undirected = (is_base_and_derived<undirected_tag,
                                            directed_category>::value
                        || is_same<undirected_tag, directed_category>::value));

  public:
    typedef std::input_iterator_tag iterator_category;
    typedef std::pair<vertices_size_type, vertices_size_type> value_type;
    typedef const value_type& reference;
    typedef const value_type* pointer;
    typedef void difference_type;

    erdos_renyi_iterator() : gen(0), n(0), edges(0), allow_self_loops(false) {}
    erdos_renyi_iterator(RandomGenerator& gen, vertices_size_type n, 
                         double prob = 0.0, bool allow_self_loops = false)
      : gen(&gen), n(n), edges(edges_size_type(prob * n * n)),
        allow_self_loops(allow_self_loops)
    { 
      if (is_undirected) edges = edges / 2;
      next(); 
    }

    reference operator*() const { return current; }
    pointer operator->() const { return &current; }
    
    erdos_renyi_iterator& operator++()
    { 
      --edges;
      next();
      return *this;
    }

    erdos_renyi_iterator operator++(int)
    {
      erdos_renyi_iterator temp(*this);
      ++(*this);
      return temp;
    }

    bool operator==(const erdos_renyi_iterator& other) const
    { return edges == other.edges; }

    bool operator!=(const erdos_renyi_iterator& other) const
    { return !(*this == other); }

  private:
    void next()
    {
      uniform_int<vertices_size_type> rand_vertex(0, n-1);
      current.first = rand_vertex(*gen);
      do {
        current.second = rand_vertex(*gen);
      } while (current.first == current.second && !allow_self_loops);
    }

    RandomGenerator* gen;
    vertices_size_type n;
    edges_size_type edges;
    bool allow_self_loops;
    value_type current;
  };

} // end namespace boost

#endif // BOOST_GRAPH_ERDOS_RENYI_GENERATOR_HPP
