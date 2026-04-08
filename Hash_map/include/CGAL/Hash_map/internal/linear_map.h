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

#ifndef CGAL_INTERNAL_LINEAR_MAP_H
#define CGAL_INTERNAL_LINEAR_MAP_H

#include <CGAL/config.h>
#include <CGAL/memory.h>

namespace CGAL {
namespace internal {

template <typename T, typename Allocator = CGAL_ALLOCATOR(T)>
class linear_map {
public:
    typedef std::size_t Key;
    typedef T Value;

    /**
     * @brief Internal storage entry for the hash map.
     */
    struct Entry {
        Key key_;
        Value value_;

        static constexpr Key null_key() {
            return (std::numeric_limits<std::size_t>::max)();
        }

        Entry() : key_(null_key()) {}

        bool is_empty() const { return key_ == null_key(); }
        void set_empty() { key_ = null_key(); }
    };

    typedef Entry* Item;

private:
    typedef std::allocator_traits<Allocator> Alloc_traits;
    typedef typename Alloc_traits::template rebind_alloc<Entry> Entry_allocator;
    typedef std::allocator_traits<Entry_allocator> Entry_alloc_traits;

    Entry_allocator allocator_;
    Entry* table_;
    std::size_t capacity_;
    std::size_t size_;
    std::size_t mask_;
    int hash_shift_;
    Value default_value_;
    std::size_t reserve_;

    /**
     * @brief Fibonacci Hashing (Multiplicative Hashing)
     *
     * Uses the golden ratio constant to map keys into the table range.
     */
    inline std::size_t hash_key(Key key) const {
        return (key * 11400714819323198485ULL) >> hash_shift_;
    }

    /**
     * @brief Computes hash shift and mask for the current capacity.
     */
    void compute_hash_params() {
        mask_ = capacity_ - 1;
        hash_shift_ = 64;
        std::size_t temp_capacity = capacity_;
        while (temp_capacity > 1) {
            temp_capacity >>= 1;
            hash_shift_--;
        }
    }

    /**
     * @brief Finds the slot for a key, or the first available empty slot.
     */
    inline Entry* find_slot(Key key) const {
        std::size_t h = hash_key(key);
        while (!table_[h].is_empty()) {
            if (table_[h].key_ == key) {
                return &table_[h];
            }
            h = (h + 1) & mask_;
        }
        return &table_[h];
    }

    Entry* allocate_table(std::size_t capacity) {
        Entry* table = Entry_alloc_traits::allocate(allocator_, capacity);
        for (std::size_t i = 0; i < capacity; ++i) {
            Entry_alloc_traits::construct(allocator_, table + i);
        }
        return table;
    }

    /**
     * @brief Deallocates the hash table.
     * @pre table must not be nullptr.
     */
    void deallocate_table(Entry* table, std::size_t capacity) {
        for (std::size_t i = 0; i < capacity; ++i) {
            Entry_alloc_traits::destroy(allocator_, table + i);
        }
        Entry_alloc_traits::deallocate(allocator_, table, capacity);
    }

    /**
     * @brief Initializes the internal hash table.
     */
    void init_table(std::size_t reserve = default_size) {
        capacity_ = default_size;
        while (capacity_ < reserve * 2) {
            capacity_ <<= 1;
        }
        compute_hash_params();
        table_ = allocate_table(capacity_);
    }

    /**
     * @brief Resizes the table and re-hashes existing entries.
     */
    void rehash(std::size_t count) {
        Entry* old_table = table_;
        std::size_t old_capacity = capacity_;

        capacity_ = count;
        compute_hash_params();

        table_ = allocate_table(capacity_);

        for (std::size_t i = 0; i < old_capacity; ++i) {
            if (!old_table[i].is_empty()) {
                Entry* slot = find_slot(old_table[i].key_);
                *slot = std::move(old_table[i]);
            }
        }

        deallocate_table(old_table, old_capacity);
    }

public:
    static constexpr std::size_t default_size = 8;

    explicit linear_map(std::size_t reserve = default_size,
                        const Value& default_val = Value(),
                        const Allocator& allocator = Allocator())
        : allocator_(allocator),
          table_(nullptr),
          capacity_(0),
          size_(0),
          mask_(0),
          hash_shift_(0),
          default_value_(default_val),
          reserve_(reserve) {
    }

    ~linear_map() {
        if (table_) deallocate_table(table_, capacity_);
    }

    linear_map(const linear_map&) = delete;
    linear_map& operator=(const linear_map&) = delete;
    linear_map(linear_map&&) = delete;
    linear_map& operator=(linear_map&&) = delete;

    Value& xdef() { return default_value_; }
    const Value& cxdef() const { return default_value_; }

    std::size_t index(Item it) const { return it->key_; }
    Value& inf(Item it) const { return it->value_; }

    /**
     * @brief Reserves space for at least count elements.
     */
    void reserve(std::size_t count) {
        if (!table_) {
            reserve_ = count;
        } else if (count * 2 > capacity_) {
            std::size_t new_capacity = capacity_;
            while (new_capacity < count * 2) {
                new_capacity <<= 1;
            }
            rehash(new_capacity);
        }
    }

    /**
     * @brief Clears all entries from the map.
     */
    void clear() {
        if (!table_) return;
        for (std::size_t i = 0; i < capacity_; ++i) {
            table_[i].set_empty();
        }
        size_ = 0;
    }

    /**
     * @brief Returns a reference to the value associated with the key,
     * inserting the default value if it doesn't exist.
     */
    Value& access(Key key) {
        if (!table_) {
            init_table(reserve_);
        } else if (size_ * 2 >= capacity_) {
            rehash(capacity_ * 2);
        }

        Entry* slot = find_slot(key);
        if (slot->is_empty()) {
            slot->key_ = key;
            slot->value_ = default_value_;
            size_++;
        }

        return slot->value_;
    }

    /**
     * @brief Finds the item associated with the key, or returns nullptr.
     */
    Item lookup(Key key) const {
        if (!table_) return nullptr;

        Entry* slot = find_slot(key);
        if (slot->is_empty()) return nullptr;

        return slot;
    }

    void statistics() const {
        if (!table_) {
            std::cout << "linear_map: not initialized (reserve=" << reserve_ << ")\n";
        } else {
            std::cout << "linear_map: capacity=" << capacity_
                      << ", size=" << size_
                      << ", load_factor=" << (double)size_ / capacity_ << "\n";
        }
    }
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_LINEAR_MAP_H
