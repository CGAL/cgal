#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <CGAL/assertions.h>

#include <map>

#include <boost/optional.hpp>

namespace CGAL::Properties {

class Property_array_base {
public:

  virtual ~Property_array_base() = default;

  // todo: Declare virtual functions here, for things which need to be done within the Property container

  virtual void reserve(std::size_t n) = 0;

  virtual void swap(std::size_t a, std::size_t b) = 0;

  virtual void reset(std::size_t i) = 0;

};

/*!
 * \brief Indexed storage for arbitrary types
 *
 * @tparam T
 */
template <typename T>
class Property_array : public Property_array_base {

  std::vector<T> m_data;
  const std::vector<bool>& m_active_indices;
  T m_default_value;

public:

  Property_array(const std::vector<bool>& active_indices, const T& default_value) :
    m_data(), m_active_indices(active_indices), m_default_value(default_value) {

    m_data.reserve(active_indices.capacity());
    m_data.resize(active_indices.size(), m_default_value);
  }

  virtual void reserve(std::size_t n) override {
    CGAL_precondition(m_active_indices.size() == n);
    m_data.resize(n, m_default_value);
  };

  virtual void swap(std::size_t a, std::size_t b) override {
    CGAL_precondition(a < m_data.size() && b < m_data.size());
    std::iter_swap(m_data.begin() + a, m_data.begin() + b);
  };

  virtual void reset(std::size_t i) override {
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

  bool operator==(const Property_array<T>& other) const {
    return &other == this;
  }

  bool operator!=(const Property_array<T>& other) const { return !operator==(other); }

};

class Property_container {

  std::map<std::string, std::shared_ptr<Property_array_base>> m_property_arrays;
  std::vector<bool> m_active_indices;

public:

  template <typename T>
  std::pair<std::reference_wrapper<Property_array<T>>, bool>
  add(const std::string& name, const T default_value = T()) {
    auto [it, created] = m_property_arrays.emplace(
      name,
      std::make_shared<Property_array<T>>(
        m_active_indices,
        default_value
      )
    );
    auto [key, array] = *it;
    auto& typed_array = dynamic_cast<Property_array<T>&>(*array);
    return {{typed_array}, !created};
  }

  template <typename T>
  const Property_array<T>& get(const std::string& name) const {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<const Property_array<T>&>(*m_property_arrays.at(name));
  }

  template <typename T>
  Property_array<T>& get(const std::string& name) {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<Property_array<T>&>(*m_property_arrays.at(name));
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

  std::size_t emplace() {
    // todo: should emplacing an element also reset it to default values?

    // If there are empty slots, return the index of one of them and mark it as full
    auto first_unused = std::find_if(m_active_indices.begin(), m_active_indices.end(), [](bool used) { return !used; });
    if (first_unused != m_active_indices.end()) {
      *first_unused = true;
      return std::distance(m_active_indices.begin(), first_unused);
    }

    // Otherwise, expand the storage and return the last element
    reserve(size() + 1);
    m_active_indices.back() = true;
    return size() - 1;

  }

  void swap(std::size_t a, std::size_t b) {
    for (auto [name, array]: m_property_arrays)
      array->swap(a, b);
  }

  void reset(std::size_t i) {
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  void erase(std::size_t i) {
    m_active_indices[i] = false;
    for (auto [name, array]: m_property_arrays)
      array->reset(i);
  }

  std::size_t size() const { return std::count(m_active_indices.begin(), m_active_indices.end(), true); }

  std::size_t capacity() const { return m_active_indices.size(); }

};

}

#endif //ORTHTREE_TESTS_PROPERTIES_H
