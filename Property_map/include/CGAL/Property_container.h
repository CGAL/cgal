// Copyright (c) 2023 INRIA
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro

#ifndef CGAL_PROPERTY_CONTAINTER_H
#define CGAL_PROPERTY_CONTAINTER_H

#include <CGAL/assertions.h>

#include <map>
#include <optional>
#include <iterator>

#include <boost/property_map/property_map.hpp>

#ifndef DOXYGEN_RUNNING

namespace CGAL::Properties::Experimental {

template <typename Index>
class Property_array_base {
public:

  Property_array_base() = default;

  Property_array_base(const Property_array_base<Index>& rhs) = delete;

  virtual ~Property_array_base() = default;

  // Declare virtual functions here, for things which need to be done within the Property container
  // todo: these should mostly be private, and made available using friend

  virtual std::shared_ptr<Property_array_base<Index>> empty_clone(const std::vector<bool>& active_indices) = 0;

  virtual std::shared_ptr<Property_array_base<Index>> clone(const std::vector<bool>& active_indices) = 0;

  virtual void copy(const Property_array_base<Index>& other) = 0;

  // desactived as MSVC 2017 as an issue with that but it is not currently used.
#if 0
  virtual void move(Property_array_base<Index>&& other) = 0;
#endif

  virtual void append(const Property_array_base<Index>& other) = 0;

  virtual void reserve(std::size_t n) = 0;

  virtual void shrink_to_fit() = 0;

  virtual void swap(Index a, Index b) = 0;

  virtual void reset(Index i) = 0;

  virtual const std::type_info& type() const = 0;

  virtual void transfer_from(const Property_array_base<Index>& other_base, Index other_index, Index this_index) = 0;

};

/*!
 * \brief Indexed storage for arbitrary types
 *
 * todo: make this effectively private, prioritize the use of Property_array_handle
 *
 * @tparam T
 */
template <typename Index, typename T>
class Property_array : public Property_array_base<Index> {

  std::vector<T> m_data;
  const std::vector<bool>& m_active_indices;
  T m_default_value;

public:

  using value_type = T;
  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Property_array(const std::vector<bool>& active_indices, const T& default_value) :
    m_data(), m_active_indices(active_indices), m_default_value(default_value) {

    m_data.reserve(active_indices.capacity());
    m_data.resize(active_indices.size(), m_default_value);
  }

  virtual std::shared_ptr<Property_array_base<Index>> empty_clone(const std::vector<bool>& active_indices) override {
    return std::make_shared<Property_array<Index, T>>(active_indices, m_default_value);
  }

  virtual std::shared_ptr<Property_array_base<Index>> clone(const std::vector<bool>& active_indices) override {
    auto new_array = std::make_shared<Property_array<Index, T>>(active_indices, m_default_value);
    new_array->m_data = m_data;
    return new_array;
  }

  virtual void copy(const Property_array_base<Index>& other_base) override {
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    m_data = other.m_data;
    CGAL_precondition(m_active_indices.size() == m_data.size());
  }

// deactived as MSVC 2017 as an issue with that but it is not currently used.
#if 0
  virtual void move(Property_array_base<Index>&& other_base) override {
    auto&& other = static_cast<Property_array<Index, T>&&>(other_base);
    m_data = std::move(other.m_data);
    CGAL_precondition(m_active_indices.size() == m_data.size());
  }
#endif

  virtual void append(const Property_array_base<Index>& other_base) override {
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    CGAL_precondition(m_data.size() + other.m_data.size() == m_active_indices.size());
    m_data.insert(m_data.end(), other.m_data.begin(), other.m_data.end());
  }

  virtual void reserve(std::size_t n) override {
    CGAL_precondition(m_active_indices.size() == n);
    m_data.resize(n, m_default_value);
  };

  virtual void shrink_to_fit() override {
    m_data.shrink_to_fit();
  }

  virtual void swap(Index a, Index b) override {
    // todo: maybe cast to index, instead of casting index to size?
    CGAL_precondition(std::size_t(a) < m_data.size() && std::size_t(b) < m_data.size());
    std::iter_swap(m_data.begin() + a, m_data.begin() + b);
  };

  virtual void reset(Index i) override {
    CGAL_precondition(std::size_t(i) < m_data.size());
    m_data[std::size_t(i)] = m_default_value;
  };

  virtual const std::type_info& type() const override { return typeid(T); };

  virtual void transfer_from(const Property_array_base<Index>& other_base,
                             Index other_index, Index this_index) override {

    CGAL_precondition(other_base.type() == type());
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    CGAL_precondition(std::size_t(other_index) < other.capacity() && std::size_t(this_index) < capacity());
    m_data[this_index] = other.m_data[other_index];
  }

public:

  // todo: there's not really a good reason to use these, maybe they should be removed

  [[nodiscard]] std::size_t size() const { return std::count(m_active_indices.begin(), m_active_indices.end(), true); }

  [[nodiscard]] std::size_t capacity() const { return m_data.size(); }

  const_reference operator[](Index i) const {
    CGAL_precondition(std::size_t(i) < m_data.size());
    return m_data[std::size_t(i)];
  }

  reference operator[](Index i) {
    CGAL_precondition(std::size_t(i) < m_data.size());
    return m_data[std::size_t(i)];
  }

  iterator begin() { return m_data.begin(); }

  iterator end() { return m_data.end(); }

  const_iterator begin() const { return m_data.begin(); }

  const_iterator end() const { return m_data.end(); }

public:

  bool operator==(const Property_array<Index, T>& other) const {
    return &other == this;
  }

  bool operator!=(const Property_array<Index, T>& other) const { return !operator==(other); }

};

// todo: property maps/array handles should go in their own file

// todo: add const/read-only handle
template <typename Index, typename T>
class Property_array_handle {

  std::reference_wrapper<Property_array<Index, T>> m_array;

public:

  // Necessary for use as a boost::property_type
  using key_type = Index;
  using value_type = T;
  using reference = typename std::vector<T>::reference;
  using const_reference = typename std::vector<T>::const_reference;
  using category = boost::lvalue_property_map_tag;

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Property_array_handle(Property_array<Index, T>& array) : m_array(array) {}

  [[nodiscard]] std::size_t size() const { return m_array.get().size(); }

  [[nodiscard]] std::size_t capacity() const { return m_array.get().capacity(); }

  Property_array<Index, T>& array() const { return m_array.get(); }

  // todo: This might not be needed, if the other operator[] is made const
  const_reference operator[](Index i) const { return m_array.get()[i]; }

  reference operator[](Index i) { return m_array.get()[i]; }

  // todo: maybe these can be const, in an lvalue property map?
  iterator begin() { return m_array.get().begin(); }

  iterator end() { return m_array.get().end(); }

  const_iterator begin() const { return m_array.get().begin(); }

  const_iterator end() const { return m_array.get().end(); }

  bool operator==(const Property_array_handle<Index, T>& other) const { return other.m_array.get() == m_array.get(); }

  bool operator!=(const Property_array_handle<Index, T>& other) const { return !operator==(other); }

  inline friend reference get(Property_array_handle<Index, T> p, const Index& i) { return p[i]; }

  inline friend void put(Property_array_handle<Index, T> p, const Index& i, const T& v) { p[i] = v; }

};

template <typename Index = std::size_t>
class Property_container {

  std::multimap<std::string, std::shared_ptr<Property_array_base<Index>>> m_properties;
  std::vector<bool> m_active_indices{};

public:

  template<typename T>
  using Array = Property_array<Index, T>;

  Property_container() = default;

  Property_container(const Property_container<Index>& other) {
    m_active_indices = other.m_active_indices;

    for (auto [name, array] : other.m_properties) {
      // todo: this could probably be made faster using emplace_hint
      m_properties.emplace(
        name,
        array->clone(m_active_indices)
      );
    }
  }

  Property_container(Property_container<Index>&& other) { *this = std::move(other); }

  // This is not exactly an assignment as existing unique properties are kept.
  Property_container<Index>& operator=(const Property_container<Index>& other) {
    m_active_indices = other.m_active_indices;

    for (auto [name, array] : other.m_properties) {
      // search if property already exists
      auto range = m_properties.equal_range(name);
      auto it = range.first;
      for (; it != range.second; it++) {
        if (typeid(*array) == typeid((*it->second)))
          break;
      }

      if (it != range.second)
        it->second->copy(*array);
      else
        m_properties.emplace(name, array->clone(m_active_indices));
    }

    return *this;
  }

  // This is not exactly an assignment as existing unique properties are kept.
  Property_container<Index>& operator=(Property_container<Index>&& other) {
    m_active_indices = std::move(other.m_active_indices);

    for (auto [name, array] : other.m_properties) {
      // search if property already exists
      auto range = m_properties.equal_range(name);
      auto it = range.first;
      for (; it != range.second; it++) {
        if (typeid(*array) == typeid((*it->second)))
          break;
      }

      if (it != range.second)
        it->second->copy(std::move(*array));
      else
        m_properties.emplace(name, array->clone(m_active_indices));
    }

    // The moved-from property map should retain all of its properties, but contain 0 elements
    other.reserve(0);
    return *this;
  }

  template <typename T>
  std::pair<std::reference_wrapper<Property_array<Index, T>>, bool>
  get_or_add_property(const std::string& name, const T default_value = T()) {
    auto range = m_properties.equal_range(name);
    for (auto it = range.first; it != range.second; it++) {
      Property_array<Index, T>* typed_array_ptr = dynamic_cast<Property_array<Index, T>*>(it->second.get());
      if (typed_array_ptr != nullptr)
        return { {*typed_array_ptr}, false };
    }

    auto it = m_properties.emplace(
      name,
      std::make_shared<Property_array<Index, T>>(
        m_active_indices,
        default_value
      )
    );

    return {{*dynamic_cast<Property_array<Index, T>*>(it->second.get())}, true};
  }

  template <typename T>
  Property_array<Index, T>& add_property(const std::string& name, const T default_value = T()) {
    // todo: I'm not settled on the naming, but it's really convenient to have a function like this
    auto [array, created] = get_or_add_property(name, default_value);
    CGAL_precondition(created);
    return array.get();
  }

/*
  // todo: misleading name, maybe it could be add_same_properties?
  void copy_properties(const Property_container<Index>& other) {
    for (auto [name, other_array]: other.m_properties) {
      // If this container doesn't have any property by this name, add it (with the same type as in other)
      if (!property_exists(name))
        m_property_arrays.emplace(name, other_array->empty_clone(m_active_indices));
    }
  }*/

  template <typename T>
  const Property_array<Index, T>& get_property(const std::string& name) const {
    return *(get_property_if_exists<T>(name));
  }

  template <typename T>
  Property_array<Index, T>& get_property(const std::string& name) {
    return *(get_property_if_exists<T>(name));
  }

  template <typename T>
  std::optional<std::reference_wrapper<Property_array<Index, T>>> get_property_if_exists(const std::string& name) const {
    auto range = m_properties.equal_range(name);
    for (auto it = range.first; it != range.second; it++) {
      Property_array<Index, T>* typed_array_ptr = dynamic_cast<Property_array<Index, T>*>(it->second.get());
      if (typed_array_ptr != nullptr)
        return *typed_array_ptr;
    }

    return {};
  }

  template <typename T>
  bool property_exists(const std::string& name) const {
    auto range = m_properties.equal_range(name);

    for (auto it = range.first; it != range.second; it++) {
      Property_array<Index, T>* typed_array_ptr = dynamic_cast<Property_array<Index, T>*>(it->second.get());
      if (typed_array_ptr != nullptr)
        return true;
    }

    return false;
  }

  /*!
   * Removes all properties with the name from the container.
   *
   * @param name
   * @return number of removed properties.
   */
  std::size_t remove_properties(const std::string& name) { return m_properties.erase(name); }

  template <typename T>
  bool remove_property(const Property_array<Index, T>& arrayToRemove) {
    const Property_array_base<Index>* ref = dynamic_cast<const Property_array_base<Index>*>(&arrayToRemove);
    for (auto it = m_properties.begin(); it != m_properties.end(); it++) {
      auto const& [name, array] = *it;
      if (array.get() == ref) {
        m_properties.erase(it);
        return true;
      }
    }
    return false;
  }

  void remove_all_properties_except(const std::vector<std::string>& preserved_names) {
    // todo: if this is used often, it should take a parameter pack instead of a vector
    // A fold expression could then be used in place of std::find for better performance
    for (auto it = m_properties.begin(); it != m_properties.end();) {
      auto const& [name, array] = *it;
      if (std::find(preserved_names.begin(), preserved_names.end(), name) == preserved_names.end())
        it = m_properties.erase(it);
      else
        it++;
    }
  }

  std::vector<std::string> properties() const {
    std::vector<std::string> property_names{};
    for (auto const& [name, _]: m_properties)
      property_names.emplace_back(name);
    return property_names;
  }

  std::size_t num_properties() const { return m_properties.size(); }

/* Deactivated as there may be several Property_maps with different types but the same name.
  const std::type_info& property_type(const std::string& name) const {
    if (auto it = m_property_arrays.find(name); it != m_property_arrays.end())
      return it->second->type();
    else
      return typeid(void);
  }*/

public:

  void reserve(std::size_t n) {
    m_active_indices.resize(n);
    for (auto [name, array]: m_properties)
      array->reserve(n);
  }

  void resize(std::size_t n) {
    reserve(n);
    std::fill(m_active_indices.begin(), m_active_indices.end(), true);
  }

  [[nodiscard]] std::size_t size() const { return std::count(m_active_indices.begin(), m_active_indices.end(), true); }

  [[nodiscard]] std::size_t capacity() const { return m_active_indices.size(); }

  Index emplace_back() {

    // Expand the storage and return the last element
    reserve(capacity() + 1);
    m_active_indices.back() = true;
    auto first_new_index = Index(capacity() - 1);
    reset(first_new_index);
    return first_new_index;
  }

  Index emplace() {

    // If there are empty slots, return the index of one of them and mark it as full
    auto first_unused = std::find_if(m_active_indices.begin(), m_active_indices.end(), [](bool used) { return !used; });
    if (first_unused != m_active_indices.end()) {
      *first_unused = true;
      auto index = Index(std::distance(m_active_indices.begin(), first_unused));
      reset(index);
      return index;
    }

    return emplace_back();
  }

  Index emplace_group_back(std::size_t n) {

    // Expand the storage and return the start of the new region
    reserve(capacity() + n);
    for (auto it = m_active_indices.end() - n; it < m_active_indices.end(); ++it)
      *it = true;
    return Index(capacity() - n);
  }

  Index emplace_group(std::size_t n) {

    auto search_start = m_active_indices.begin();
    while (search_start != m_active_indices.end()) {

      // Find the first unused cell
      auto unused_begin = std::find_if(
        search_start, m_active_indices.end(),
        [](bool used) { return !used; }
      );

      auto unused_end = unused_begin;

      // Determine if the group fits
      if (std::distance(unused_begin, m_active_indices.end()) >= static_cast<typename std::iterator_traits<std::vector<bool>::iterator>::difference_type>(n))
        unused_end = std::find_if(
          unused_begin, (std::min)(unused_begin + n, m_active_indices.end()),
          [](bool used) { return used; }
      );

      // If the discovered range was large enough
      if (std::distance(unused_begin, unused_end) >= static_cast<typename std::iterator_traits<std::vector<bool>::iterator>::difference_type>(n)) {

        // Mark the indices as used, and reset the properties of each of them
        // todo: it would be better to provide a function to set a range
        for (auto it = unused_begin; it < unused_end; ++it) {
          *it = true;
          reset(Index(std::distance(m_active_indices.begin(), it)));
        }

        // Return the first index of the range
        return Index(std::distance(m_active_indices.begin(), unused_begin));
      }

      // If we didn't find a large enough region, continue our search after the end
      search_start = unused_end;
    }

    // If no empty regions were found, expand the storage
    return emplace_group_back(n);
  }

  void swap(Index a, Index b) {
    for (auto [name, array]: m_properties)
      array->swap(a, b);
  }

  void reset(Index i) {
    for (auto [name, array]: m_properties)
      array->reset(i);
  }

  void erase(Index i) {
    m_active_indices[i] = false;
    for (auto [name, array]: m_properties)
      array->reset(i);
  }

  bool is_erased(Index i) const {
    return !m_active_indices[i];
  }

  // todo: I'd prefer to eliminate this, if possible
  void mark_active(Index i) {
    return m_active_indices[i] = true;
  }

  void mark_inactive(Index i) {
    return m_active_indices[i] = false;
  }

  std::vector<Index> active_list() const {
    std::vector<Index> indices;
    for (std::size_t i = 0; i < m_active_indices.size(); ++i)
      if (m_active_indices[i]) indices.emplace_back(i);
    return indices;
  }

  std::vector<Index> inactive_list() const {
    std::vector<Index> indices;
    for (std::size_t i = 0; i < m_active_indices.size(); ++i)
      if (!m_active_indices[i]) indices.emplace_back(i);
    return indices;
  }

  void shrink_to_fit() {
    for (auto [name, array]: m_properties)
      array->shrink_to_fit();
  }

  /*!
   * Adds the elements of the other container to this container for each property which is present in this container.
   *
   * Gaps in both containers are preserved, and all elements of the other container are guaranteed
   * to appear after the elements of this container.
   * Properties in this container which don't appear in the other container are extended with default values.
   * Properties in the other container which don't appear in this one are not included.
   * todo: merge() would be useful as well, but could break contiguous regions in the other container
   *
   * @param other
   */
  void append(const Property_container<Index>& other) {

    m_active_indices.insert(m_active_indices.end(), other.m_active_indices.begin(), other.m_active_indices.end());
    for (auto [name, array]: m_properties) {

      auto range = other.m_properties.equal_range(name);
      auto it = range.first;
      for (; it != range.second; it++) {
        if (typeid(array.get()) == typeid((it->second.get())))
          break;
      }

      if (it != range.second)
        array->append(*it->second.get());
      else
        array->reserve(m_active_indices.size());
    }
  }

/*
  // todo: maybe should be renamed to transfer_from, but I'd rather remove this functionality entirely
  void transfer(const Property_container<Index>& other, Index other_index, Index this_index) {
    CGAL_precondition(other.m_property_arrays.size() == m_property_arrays.size());
    for (auto [name, array]: m_property_arrays) {
      auto other_array = other.m_property_arrays.at(name);
      array->transfer_from(*other_array, other_index, this_index);
    }
  }*/

  // todo: maybe a compress() method?
};

}

#endif // DOXYGEN_RUNNING

#endif //CGAL_PROPERTY_CONTAINTER_H
