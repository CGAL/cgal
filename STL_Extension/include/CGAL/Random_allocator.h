// Copyright (c) 2025 GeometryFactory Sarl (France).
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_RANDOM_ALLOCATOR_H
#define CGAL_RANDOM_ALLOCATOR_H

#include <boost/container/static_vector.hpp>
#include <boost/dynamic_bitset.hpp>
#include <memory>
#include <vector>

// This header requires C++20 or later.
#include <format>
#include <random>
#include <ranges>
#include <source_location>

namespace CGAL {

// A random allocator that allocates blocks of memory and returns random
// locations. To use only for debugging purposes. That allocator is:
//   - not efficient,
//   - not thread-safe,
//   - and increases memory-fragmentation and non-locality.
template <typename T, typename Upstream_allocator = std::allocator<T>> class Random_allocator
{
  constexpr static std::size_t minimal_block_size = 1024;
  constexpr static std::size_t random_size = 64;

  struct Blocks
  {
    struct Block
    {
      T* data = nullptr;
      boost::dynamic_bitset<> available;
      std::size_t maximal_continuous_free_space = 0;

      auto size() const { return available.size(); }
    };

    Upstream_allocator alloc;
    std::vector<Block> blocks; // List of allocated blocks

    void allocate_new_block(std::size_t block_size = minimal_block_size)
    {
      auto& block = blocks.emplace_back(nullptr, boost::dynamic_bitset<>(block_size));
      block.data = alloc.allocate(block_size);
      block.available.set();
      block.maximal_continuous_free_space = block_size;
    }

    Blocks(const Upstream_allocator& alloc) : alloc(alloc) { allocate_new_block(); }

    ~Blocks()
    {
      for(auto& block : blocks) {
        alloc.deallocate(block.data, block.size());
      }
    }
  };
  using Block = typename Blocks::Block;
  using Ptr = std::shared_ptr<Blocks>;

  Ptr ptr_;
  std::mt19937 gen;

public:
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using allocator_type = Upstream_allocator;

  Random_allocator(const Upstream_allocator& alloc = {}) noexcept;

  pointer allocate(size_type n, const void* hint = 0);
  void deallocate(pointer p, size_type n);

  template <typename U, typename... Args> void construct(U* p, Args&&... args);
  template <typename U> void destroy(U* p);

  size_type max_size() const noexcept;

  template <typename U> struct rebind
  {
    using other_upstream_allocator = std::allocator_traits<Upstream_allocator>::template rebind_alloc<U>;
    using other = Random_allocator<U, other_upstream_allocator>;
  };

  bool operator==(const Random_allocator&) const noexcept;
  bool operator!=(const Random_allocator&) const noexcept;
};

// Implementation of Random_allocator methods
template <typename T, typename Upstream_allocator>
Random_allocator<T, Upstream_allocator>::Random_allocator(const Upstream_allocator& alloc) noexcept
    : ptr_(std::make_shared<Blocks>(alloc)), gen(std::random_device()())
{
}

template <typename T, typename Upstream_allocator>
typename Random_allocator<T, Upstream_allocator>::pointer
Random_allocator<T, Upstream_allocator>::allocate(size_type n, const void* hint)
{
  boost::container::static_vector<std::pair<Block*, size_type>, random_size> found_spaces;
  for(auto& block : ptr_->blocks | std::views::reverse) {
    if(block.maximal_continuous_free_space < n)
      continue;
    auto& available = block.available;
    const auto block_size = block.size();
    size_type found_max_free_space = 0;
    auto index = available.find_first();
    while(index != boost::dynamic_bitset<>::npos) {
      available.flip();
      const auto end_of_free_block = (std::min)(available.find_next(index), block_size);
      available.flip();
      auto free_space = end_of_free_block - index;
      found_max_free_space = (std::max)(found_max_free_space, free_space);
      while(free_space > n && found_spaces.size() < found_spaces.capacity()) {
        found_spaces.push_back({std::addressof(block), index});
        free_space -= n;
        index += n;
      }
      index = block.available.find_next(end_of_free_block);
    }
    block.maximal_continuous_free_space = found_max_free_space;
    if(found_spaces.size() == found_spaces.capacity())
      break;
  }
  if(!found_spaces.empty()) {
    std::uniform_int_distribution<> dis(0, found_spaces.size() - 1);
    auto i = dis(gen);
    auto [block, index] = found_spaces[i];
    block->available.set(index, n, false);
#if CGAL_DEBUG_RANDOM_ALLOCATOR
    std::clog << std::format("CGAL::Random_allocator debug info: n = {}, found_spaces.size() = {}, i = {},"
                             "block nb = {}, block size = {}, index = {}\n",
                             n, found_spaces.size(), i, block - ptr_->blocks.data(), block->size(), index);
#endif // CGAL_DEBUG_RANDOM_ALLOCATOR
    return block->data + index;
  }
  size_type block_size = std::max(n * random_size, minimal_block_size);
  ptr_->allocate_new_block(block_size);
  return allocate(n, hint);
}

template <typename T, typename Upstream_allocator>
void Random_allocator<T, Upstream_allocator>::deallocate(pointer p, size_type n)
{
  for(auto& block : ptr_->blocks) {
    if(block.data <= p && p < block.data + block.size()) {
      block.available.set(p - block.data, n, true);
      return;
    }
  }
}

template <typename T, typename Upstream_allocator>
template <typename U, typename... Args>
void Random_allocator<T, Upstream_allocator>::construct(U* p, Args&&... args)
{
  ::new((void*)p) U(std::forward<Args>(args)...);
}

template <typename T, typename Upstream_allocator>
template <typename U>
void Random_allocator<T, Upstream_allocator>::destroy(U* p)
{
  p->~U();
}

template <typename T, typename Upstream_allocator>
typename Random_allocator<T, Upstream_allocator>::size_type
Random_allocator<T, Upstream_allocator>::max_size() const noexcept
{
  return std::numeric_limits<size_type>::max() / sizeof(T);
}

} // namespace CGAL

#endif // CGAL_RANDOM_ALLOCATOR_H
