#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <CGAL/assertions.h>

#include <map>

#include <boost/optional.hpp>

namespace CGAL::Properties {

template <typename Index>
class Property_array_base {
public:

  virtual ~Property_array_base() = default;

  // todo: Declare virtual functions here, for things which need to be done within the Property container

  virtual void append(const Property_array_base<Index>& other) = 0;

  virtual void reserve(std::size_t n) = 0;

  virtual void swap(Index a, Index b) = 0;

  virtual void reset(Index i) = 0;

};

/*!
 * \brief Indexed storage for arbitrary types
 *
 * @tparam T
 */
template <typename Index, typename T>
class Property_array : public Property_array_base<Index> {

  std::vector<T> m_data;
  const std::vector<bool>& m_active_indices;
  T m_default_value;

public:

  Property_array(const std::vector<bool>& active_indices, const T& default_value) :
    m_data(), m_active_indices(active_indices), m_default_value(default_value) {

    m_data.reserve(active_indices.capacity());
    m_data.resize(active_indices.size(), m_default_value);
  }

  virtual void append(const Property_array_base<Index>& other_base) {
    auto& other = dynamic_cast<const Property_array<Index, T>&>(other_base);
    CGAL_precondition(m_data.size() + other.m_data.size() == m_active_indices.size());
    m_data.insert(m_data.end(), other.m_data.begin(), other.m_data.end());
  }

  virtual void reserve(std::size_t n) override {
    CGAL_precondition(m_active_indices.size() == n);
    m_data.resize(n, m_default_value);
  };

  virtual void swap(Index a, Index b) override {
    CGAL_precondition(a < m_data.size() && b < m_data.size());
    std::iter_swap(m_data.begin() + a, m_data.begin() + b);
  };

  virtual void reset(Index i) override {
    CGAL_precondition(i < m_data.size());
    m_data[i] = m_default_value;
  };

  std::size_t capacity() const { return m_data.size(); }

public:

  const T& operator[](std::size_t i) const {
    CGAL_precondition(i < m_data.size());
    return m_data[i];
  }

  T& operator[](std::size_t i) {
    CGAL_precondition(i < m_data.size());
    return m_data[i];
  }

public:

  bool operator==(const Property_array<Index, T>& other) const {
    return &other == this;
  }

  bool operator!=(const Property_array<Index, T>& other) const { return !operator==(other); }

};

template <typename Index = std::size_t>
class Property_container {

  std::map<std::string, std::shared_ptr<Property_array_base<Index>>> m_property_arrays;
  std::vector<bool> m_active_indices;

public:

  template <typename T>
  std::pair<std::reference_wrapper<Property_array<Index, T>>, bool>
  add(const std::string& name, const T default_value = T()) {
    auto [it, created] = m_property_arrays.emplace(
      name,
      std::make_shared<Property_array<Index, T>>(
        m_active_indices,
        default_value
      )
    );
    auto [key, array] = *it;
    auto& typed_array = dynamic_cast<Property_array<Index, T>&>(*array);
    return {{typed_array}, !created};
  }

  template <typename T>
  const Property_array<Index, T>& get(const std::string& name) const {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<const Property_array<Index, T>&>(*m_property_arrays.at(name));
  }

  template <typename T>
  Property_array<Index, T>& get(const std::string& name) {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<Property_array<Index, T>&>(*m_property_arrays.at(name));
  }

  /*!
   * Removes a property array from the container
   *
   * @param name
   * @return True if a container with this name existed, false otherwise
   */
  bool remove(const std::string& name) { return m_property_arrays.erase(name) == 1; }

  std::size_t n_properties() { return m_property_arrays.size(); }

public:

  void reserve(std::size_t n) {
    m_active_indices.resize(n);
    for (auto [name, array]: m_property_arrays)
      array->reserve(n);
  }

  std::size_t size() const { return std::count(m_active_indices.begin(), m_active_indices.end(), true); }

  std::size_t capacity() const { return m_active_indices.size(); }

  Index emplace() {

    // If there are empty slots, return the index of one of them and mark it as full
    auto first_unused = std::find_if(m_active_indices.begin(), m_active_indices.end(), [](bool used) { return !used; });
    if (first_unused != m_active_indices.end()) {
      *first_unused = true;
      return std::distance(m_active_indices.begin(), first_unused);
    }

    // Otherwise, expand the storage and return the last element
    reserve(capacity() + 1);
    m_active_indices.back() = true;
    // todo: should emplacing an element also reset it to default values?
    reset(capacity() - 1);
    return capacity() - 1;

  }

  Index emplace_group(std::size_t n) {

    auto search_start = m_active_indices.begin();
    while (search_start != m_active_indices.end()) {

      // Find the first unused cell
      auto unused_begin = std::find_if(
        search_start, m_active_indices.end(),
        [](bool used) { return !used; }
      );

      // Determine if the group fits
      auto unused_end = std::find_if(
        unused_begin, std::min(unused_begin + n, m_active_indices.end()),
        [](bool used) { return used; }
      );

      // If the discovered range was large enough
      if (std::distance(unused_begin, unused_end) >= n) {

        // Mark the indices as used, and reset the properties of each of them
        // todo: it would be better to provide a function to set a range
        for (auto it = unused_begin; it < unused_end; ++it) {
          *it = true;
          reset(std::distance(m_active_indices.begin(), it));
        }

        // Return the first index of the range
        return std::distance(m_active_indices.begin(), unused_begin);
      }

      // If we didn't find a large enough region, continue our search after the end
      search_start = unused_end;
    }

    // If no empty regions were found, expand the storage
    reserve(capacity() + n);
    for (auto it = m_active_indices.end() - n; it < m_active_indices.end(); ++it)
      *it = true;
    return capacity() - n;

  }

  void swap(Index a, Index b) {
    for (auto [name, array]: m_property_arrays)
      array->swap(a, b);
  }

  void reset(Index i) {
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  void erase(Index i) {
    m_active_indices[i] = false;
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  /*!
   * Adds the elements of the other container to this container for each property which is present in this container.
   *
   * Gaps are preserved, and all elements of the other container are guaranteed
   * to appear after the elements of this container.
   * todo: merge() would be useful as well, but could break contiguous regions in the other container
   *
   * @param other
   */
  void append(const Property_container<Index>& other) {
    // todo

    m_active_indices.insert(m_active_indices.end(), other.m_active_indices.begin(), other.m_active_indices.end());
    for (auto [name, array]: m_property_arrays) {
      auto it = other.m_property_arrays.find(name);
      if (it != other.m_property_arrays.end())
        array->append(*it->second);
      else
        array->reserve(m_active_indices.size());
    }
  }
};

}

#endif //ORTHTREE_TESTS_PROPERTIES_H
