// Copyright (c) 2021  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Giles Bathgate <giles.bathgate@gmail.com>

#ifndef CGAL_HANDLE_HASH_MAP
#define CGAL_HANDLE_HASH_MAP

#include <CGAL/config.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/Tools/robin_hood.h>
#include <cstddef>
#include <functional>

namespace CGAL
{
namespace internal
{

template <class Key,
          class T,
          class Hash = Handle_hash_function,
          class KeyEqual = std::equal_to<Key>,
          class Allocator = void>
class Handle_hash_map
  : private robin_hood::unordered_node_map<Key,T,Hash,KeyEqual>
{
    typedef robin_hood::unordered_node_map<Key,T,Hash,KeyEqual> Base;

public:
    // STL compliance
    typedef Key                           key_type;
    typedef T                             mapped_type;
    typedef std::size_t                   size_type;
    typedef Hash                          hasher;
    typedef KeyEqual                      key_equal;
    typedef typename Base::const_iterator const_iterator;

    // Backwards compatibility
    typedef mapped_type  Data;

private:
    mapped_type _default;

public:
    Handle_hash_map()
        : _default(mapped_type())
    {
    }

    Handle_hash_map(const mapped_type& def)
        : _default(def)
    {
    }

    Handle_hash_map(const mapped_type& def,size_type table_size)
        : _default(def)
    {
        Base::reserve(table_size);
    }

    Handle_hash_map(key_type key_begin,key_type key_end,mapped_type values_begin,
                    const mapped_type& def,size_type table_size)
        : _default(def)
    {
        Base::reserve(table_size);
        insert(key_begin,key_end,values_begin);
    }

    void clear()
    {
        Base::clear();
    }

    void clear(const mapped_type& def)
    {
        clear();
        _default = def;
    }

    bool is_defined(const key_type& key) const
    {
        return Base::contains(key);
    }

    const mapped_type& operator[](const key_type& key) const
    {
        const_iterator f = Base::find(key);
        if(f != Base::end())
            return f->second;
        return _default;
    }

    mapped_type& operator[](const key_type& key)
    {
        return Base::try_emplace(key,_default).first->second;
    }

    mapped_type insert(key_type key_begin,key_type key_end,
                       mapped_type values_begin)
    {
        for(; key_begin != key_end; (++key_begin,++values_begin)) {
            Base::insert_or_assign(key_begin,values_begin);
        }
        return values_begin;
    }
};

}
}

#endif // CGAL_HANDLE_HASH_MAP

