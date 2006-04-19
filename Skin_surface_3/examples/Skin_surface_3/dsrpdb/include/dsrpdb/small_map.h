/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */


#ifndef DSR_PDB_SMALL_MAP_H
#define DSR_PDB_SMALL_MAP_H

#include <vector>
#include <algorithm>

namespace dsrpdb {
  template <class Key, class Data>
  class small_map {
  public:
    typedef std::pair<Key, Data> value_type;
    typedef Key key_type;
    typedef Data data_type;
    typedef std::vector<value_type> container;
    typedef typename container::iterator iterator;
    typedef typename container::const_iterator const_iterator;


    small_map(std::size_t sz=0){c_.reserve(sz);}

    iterator find(key_type k) {
      for (iterator it= c_.begin(); it != c_.end(); ++it){
	if (it->first==k) return it;
      }
      return end();
    }
    const_iterator find(key_type k) const {
      for (const_iterator it= c_.begin(); it != c_.end(); ++it){
	if (it->first==k) return it;
      }
      return end();
    }

    iterator begin() {
      return c_.begin();
    }

    iterator end() {
      return c_.end();
    }

    const_iterator begin() const {
      return c_.begin();
    }

    const_iterator end() const {
      return c_.end();
    }

    data_type& operator[](key_type k){
      iterator it= find(k);
      if (it != end()) return it->second;
      else {
	c_.push_back(value_type(k,data_type()));
	return c_.back().second;
      }
    }

    void insert(const value_type &v) {
      c_.push_back(v);
    }

    std::size_t size() const {
      return c_.size();
    }

  protected:
	container c_;
  };
}
#endif
