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

#ifndef BOOST_BITSET_ADAPTOR_HPP
#define BOOST_BITSET_ADAPTOR_HPP

    template <class T, class Derived>
    struct bitset_adaptor {
      Derived& derived() { return static_cast<Derived&>(*this); }
      const Derived& derived() const { 
        return static_cast<const Derived&>(*this); 
      }
    };

    template <class T, class D, class V>
    bool set_contains(const bitset_adaptor<T,D>& s, const V& x) {
      return s.derived().test(x);
    }
    
    template <class T, class D>
    bool set_equal(const bitset_adaptor<T,D>& x,
                   const bitset_adaptor<T,D>& y) {
      return x.derived() == y.derived();
    }

    template <class T, class D>
    int set_lex_order(const bitset_adaptor<T,D>& x,
                      const bitset_adaptor<T,D>& y) {
      return compare_3way(x.derived(), y.derived());
    }

    template <class T, class D>
    void set_clear(bitset_adaptor<T,D>& x) {
      x.derived().reset();
    }

    template <class T, class D>
    bool set_empty(const bitset_adaptor<T,D>& x) {
      return x.derived().none();
    }

    template <class T, class D, class V>
    void set_insert(bitset_adaptor<T,D>& x, const V& a) {
      x.derived().set(a);
    }

    template <class T, class D, class V>
    void set_remove(bitset_adaptor<T,D>& x, const V& a) {
      x.derived().set(a, false);
    }
    
    template <class T, class D>    
    void set_intersect(const bitset_adaptor<T,D>& x,
                       const bitset_adaptor<T,D>& y,
                       bitset_adaptor<T,D>& z)
    {
      z.derived() = x.derived() & y.derived();
    }

    template <class T, class D>    
    void set_union(const bitset_adaptor<T,D>& x,
                   const bitset_adaptor<T,D>& y,
                   bitset_adaptor<T,D>& z)
    {
      z.derived() = x.derived() | y.derived();
    }

    template <class T, class D>    
    void set_difference(const bitset_adaptor<T,D>& x,
                        const bitset_adaptor<T,D>& y,
                        bitset_adaptor<T,D>& z)
    {
      z.derived() = x.derived() - y.derived();
    }

    template <class T, class D>    
    void set_compliment(const bitset_adaptor<T,D>& x,
                        bitset_adaptor<T,D>& z)
    {
      z.derived() = x.derived();
      z.derived().flip();
    }
    
#endif // BOOST_BITSET_ADAPTOR_HPP
