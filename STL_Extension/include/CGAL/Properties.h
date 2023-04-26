#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <CGAL/assertions.h>

#include <map>

#include <boost/optional.hpp>

namespace CGAL::Properties {

template <char... chars>
using String_literal_type = std::integer_sequence<char, chars...>;

template <typename T, T... chars>
constexpr String_literal_type<chars...> operator ""_param() { return {}; }

class Property_array_base {
public:

  virtual ~Property_array_base() = default;

  // todo: Declare virtual functions here for things which need to be done within the Property container

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

public:

  Property_array(const std::vector<bool>& active_indices, const T& default_value) :
    m_active_indices(active_indices), m_data() {}

  bool operator==(const Property_array<T>& other) const {
    return &other == this;
  }

  bool operator!=(const Property_array<T>& other) const { return !operator==(other); }

};

class Property_container {


  std::map<std::string, std::shared_ptr<Property_array_base>> m_property_arrays;
  std::vector<bool> m_active_indices;

public:

  /*
  template <typename T>
  Property_array<T>& add(const std::string& name, const T default_value = T()) {
    return dynamic_cast<Property_array<T>&>(*(*m_property_arrays.emplace(name, std::make_shared<Property_array<T>>(
      m_active_indices,
      default_value)
    ).first).second);
  }

  template <typename T>
  Property_array<T>& get(const std::string& name) {
    CGAL_precondition(m_property_arrays.count(name) != 0);
    return dynamic_cast<Property_array<T>&>(*m_property_arrays[name]);
  }
  */

  template <typename T, typename Name>
  Property_array<T>& get() {

    static Property_array<T> array{m_active_indices, T{}};
    return array;

  }

  template <typename T, typename Name>
  std::pair<std::reference_wrapper<Property_array<T>>, std::unique_lock<std::mutex>> threadsafe_get() {
    static std::mutex m{};
    return {std::reference_wrapper<Property_array<T>>{get<T, Name>()}, std::unique_lock<std::mutex>{m}};
  }


private:

  static constexpr std::size_t next_type_index() {
    static std::atomic<std::size_t> value{0};
    return value++;
  }

  template <typename T, typename Name>
  static constexpr std::size_t type_index() {
    static std::size_t value{next_type_index()};
    return value;
  }

};

}

#endif //ORTHTREE_TESTS_PROPERTIES_H
