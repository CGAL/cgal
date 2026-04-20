// Copyright (c) 2026 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Giles Bathgate

#ifndef CGAL_INTERNAL_FLAT_MAP_H
#define CGAL_INTERNAL_FLAT_MAP_H

#include <CGAL/config.h>
#include <CGAL/memory.h>
#include <CGAL/unordered_flat_map.h>

#include <cstddef>
#include <memory>
#include <utility>

namespace CGAL {
namespace internal {

/**
 * @brief A Unique_hash_map backend that wraps CGAL::unordered_flat_map.
 */
template <typename Value, typename Allocator = CGAL_ALLOCATOR(Value)>
class unordered_flat_map_adaptor {

    struct identity_hash {
        std::size_t operator()(std::size_t x) const noexcept { return x; }
    };

public:
    typedef std::size_t Key;
    typedef std::pair<const Key, Value> Entry;
    typedef typename std::allocator_traits<Allocator>::template rebind_alloc<Entry> Entry_allocator;
    typedef CGAL::unordered_flat_map<Key, Value, identity_hash, std::equal_to<Key>, Entry_allocator> Map;
    typedef const Entry* Item;
    typedef unordered_flat_map_adaptor<Value, Allocator> Self;

    static constexpr std::size_t default_size = 8;

    unordered_flat_map_adaptor(std::size_t reserve = default_size,
             const Value& default_val = Value(),
             const Allocator& allocator = Allocator())
        : map_(reserve, identity_hash(), std::equal_to<Key>(), allocator),
          default_value_(default_val) {
    }

    unordered_flat_map_adaptor(const Self& other)
        : map_(other.map_),
          default_value_(other.default_value_) {
    }

    Self& operator=(const Self& other) {
        if (this != &other) {
            map_ = other.map_;
            default_value_ = other.default_value_;
        }
        return *this;
    }

    unordered_flat_map_adaptor(Self&& other) noexcept
        : map_(std::move(other.map_)),
          default_value_(std::move(other.default_value_)) {
    }

    Self& operator=(Self&& other) noexcept {
        if (this != &other) {
            map_ = std::move(other.map_);
            default_value_ = std::move(other.default_value_);
        }
        return *this;
    }

    Value& xdef() { return default_value_; }
    const Value& cxdef() const { return default_value_; }

    Key index(Item it) const { return it->first; }
    const Value& inf(Item it) const { return it->second; }

    void reserve(std::size_t count) { map_.reserve(count); }
    void clear() { map_.clear(); }

    Value& access(Key key) {
        auto [it, inserted] = map_.try_emplace(key, default_value_);
        return it->second;
    }

    Item lookup(Key key) const {
        auto it = map_.find(key);
        if (it != map_.end()) return &*it;
        return nullptr;
    }

    void statistics() const {
        std::cout << "flat_map: size=" << map_.size() << "\n";
    }

private:
    Map map_;
    Value default_value_;
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_FLAT_MAP_H
