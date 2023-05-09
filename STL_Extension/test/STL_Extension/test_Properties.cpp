
#include <CGAL/Properties.h>

using namespace CGAL::Properties;

void test_property_creation() {

  Property_container properties;

  // Should return an integer array which didn't previously exist
  auto [integers, preexisting] = properties.add("integer", 5);
  static_assert(std::is_same_v<decltype(integers), std::reference_wrapper<Property_array<std::size_t, int>>>);
  assert(!preexisting);
  assert(properties.n_properties() == 1);

  auto [floats, _] = properties.add<float>("float");
  static_assert(std::is_same_v<decltype(floats), std::reference_wrapper<Property_array<std::size_t, float>>>);
  assert(properties.n_properties() == 2);

  // get() should retreive the same arrays
  assert(integers.get() == properties.get<int>("integer"));
  assert(floats.get() == properties.get<float>("float"));

  // remove() should delete a property array & return if it existed
  assert(!properties.remove("not-a-real-property"));
  auto removed = properties.remove("integer");
  assert(removed);
  assert(properties.n_properties() == 1);

}

void test_property_array() {

  Property_container properties;

  auto [integers, integers_existed] = properties.add("integers", 5);

  // Reserve space for 100 elements
  properties.reserve(100);
  assert(properties.capacity() == 100);
  assert(properties.size() == 0);

  // Newly emplaced elements should go at the front
  assert(properties.emplace() == 0);
  assert(properties.emplace() == 1);
  assert(properties.emplace() == 2);
  assert(properties.size() == 3);

  // Make sure that the new elements are equal to the default value
  assert(integers.get()[0] == 5);
  assert(integers.get()[1] == 5);
  assert(integers.get()[2] == 5);

  // Add a new property
  auto [floats, floats_existed] = properties.add("floats", 6.0f);

  // The new property array should already be of the right size
  assert(floats.get().capacity() == 100);
  assert(properties.size() == 3);

  // Pre-existing elements should contain the default value
  assert(floats.get()[0] == 6.0f);
  assert(floats.get()[1] == 6.0f);
  assert(floats.get()[2] == 6.0f);

  // Update values for a few elements
  floats.get()[0] = 1.0f;
  floats.get()[1] = 2.0f;
  floats.get()[2] = 3.0f;
  integers.get()[2] = -2;
  assert(floats.get()[0] == 1.0f);
  assert(floats.get()[1] == 2.0f);
  assert(floats.get()[2] == 3.0f);
  assert(integers.get()[2] == -2);

  // Reset an element, and all of its properties should revert to the defaults
  properties.reset(2);
  assert(floats.get()[2] == 6.0f);
  assert(integers.get()[2] == 5);

  // Erase an element, and the size should be reduced
  properties.erase(1);
  assert(properties.size() == 2);
  assert(properties.capacity() == 100);

  // A newly emplaced element should take the empty slot
  assert(properties.emplace() == 1);
  assert(properties.size() == 3);
  // todo: should the new element have default properties?
  assert(properties.emplace() == 3);
  assert(properties.size() == 4);

  // Swapping a pair of elements swaps all of their properties
  properties.swap(0, 3);
  assert(integers.get()[0] == 5);
  assert(floats.get()[0] == 6.0f);
  assert(integers.get()[3] == 5);
  assert(floats.get()[3] == 1.0f);
}

int main() {

  test_property_creation();
  test_property_array();

  return 0;
}