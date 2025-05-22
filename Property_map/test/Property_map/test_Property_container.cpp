
#include <CGAL/Property_container.h>
#include <CGAL/use.h>

using namespace CGAL::Properties::Experimental;

void test_property_creation() {

  Property_container properties;

  // Should return an integer array which didn't previously exist
  auto [integers, created] = properties.get_or_add_property("integer", 5);
  static_assert(std::is_same_v<decltype(integers), std::reference_wrapper<Property_array<std::size_t, int>>>);
  assert(created);
  assert(properties.num_properties() == 1);

  auto [floats, _] = properties.get_or_add_property<float>("float");
  static_assert(std::is_same_v<decltype(floats), std::reference_wrapper<Property_array<std::size_t, float>>>);
  assert(properties.num_properties() == 2);

  // get() should retrieve the same arrays
  assert(integers.get() == properties.get_property<int>("integer"));
  assert(floats.get() == properties.get_property<float>("float"));

  // remove() should delete a property array & return if it existed
  assert(!properties.remove_properties("not-a-real-property"));
  auto removed = properties.remove_property<int>(integers);
  assert(removed);
  assert(properties.num_properties() == 1);

  // Add a new property
  auto [bools, bools_created] = properties.get_or_add_property<bool>("bools", false);
  static_assert(std::is_same_v<decltype(bools), std::reference_wrapper<Property_array<std::size_t, bool>>>);
  Property_array<std::size_t, bool>& b = bools.get();
  CGAL_USE(b);
}

void test_element_access() {

  Property_container properties;

  auto& integers = properties.add_property("integers", 5);

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
  assert(integers[0] == 5);
  assert(integers[1] == 5);
  assert(integers[2] == 5);

  // Add a new property
  auto& floats = properties.add_property("floats", 6.0f);

  // The new property array should already be of the right size
  assert(floats.capacity() == 100);
  assert(properties.size() == 3);

  // Pre-existing elements should contain the default value
  assert(floats[0] == 6.0f);
  assert(floats[1] == 6.0f);
  assert(floats[2] == 6.0f);

  // Update values for a few elements
  floats[0] = 1.0f;
  floats[1] = 2.0f;
  floats[2] = 3.0f;
  integers[2] = -2;
  assert(floats[0] == 1.0f);
  assert(floats[1] == 2.0f);
  assert(floats[2] == 3.0f);
  assert(integers[2] == -2);

  // Reset an element, and all of its properties should revert to the defaults
  properties.reset(2);
  assert(floats[2] == 6.0f);
  assert(integers[2] == 5);

  // Erase an element, and the size should be reduced
  properties.erase(1);
  assert(properties.size() == 2);
  assert(properties.capacity() == 100);
  assert(properties.active_list().size() == 2);
  assert(properties.inactive_list().size() == 98);

  // A newly emplaced element should take the empty slot
  assert(properties.emplace() == 1);
  assert(properties.size() == 3);
  // todo: should the new element have default properties?
  assert(properties.emplace() == 3);
  assert(properties.size() == 4);

  // Swapping a pair of elements swaps all of their properties
  properties.swap(0, 3);
  assert(integers[0] == 5);
  assert(floats[0] == 6.0f);
  assert(integers[3] == 5);
  assert(floats[3] == 1.0f);

}

void test_emplace_group() {

  Property_container properties;

  auto& a = properties.add_property("a", 5);
  CGAL_USE(a);
  // Insert a group of 100 elements
  properties.emplace_group(100);
  assert(properties.size() == 100);

  // Eliminate a few regions
  properties.erase(3);
  assert(properties.is_erased(3));
  assert(properties.size() == 99);
  for (int i = 20; i < 25; ++i)
    properties.erase(i);
  assert(properties.is_erased(23));
  assert(properties.size() == 94);
  for (int i = 50; i < 80; ++i)
    properties.erase(i);
  assert(properties.is_erased(53));
  assert(properties.size() == 64);

  // A group of size 4 should only fit in the empty region fo size 5
  assert(properties.emplace_group(4) == 20);
  assert(properties.size() == 68);
  assert(properties.capacity() == 100);

  // A group of size 16 should only fit in the empty region fo size 30
  assert(properties.emplace_group(16) == 50);
  assert(properties.size() == 84);
  assert(properties.capacity() == 100);

  // Another group of size 16 should require the storage to expand, because the largest empty region is mostly full now
  assert(properties.emplace_group(16) == 100);
  assert(properties.size() == 100);
  assert(properties.capacity() == 116);

}

void test_append() {

  // Create a pair of property containers with similar contents
  Property_container properties_a, properties_b;
  properties_a.add_property("ints", 1);
  properties_b.add_property("ints", 2);
  properties_a.add_property("floats", 3.0f);
  properties_b.add_property("floats", 4.0f);

  // One container will also contain an extra property
  properties_a.add_property("bools", true);

  // Add some values to both property sets
  properties_a.emplace_group(10);
  properties_b.emplace_group(5);
  assert(properties_a.size() == 10);
  assert(properties_b.size() == 5);

  // Add the second group to the end of the first
  properties_a.append(properties_b);
  assert(properties_a.size() == 15);
  assert(properties_b.size() == 5);

  // Initialized values from the second group should appear after those of the first
  assert(properties_a.get_property<int>("ints")[5] == 1);
  assert(properties_a.get_property<int>("ints")[12] == 2);
  assert(properties_a.get_property<float>("floats")[5] == 3.0f);
  assert(properties_a.get_property<float>("floats")[12] == 4.0f);

  // Additional properties in the first group should have expanded too, and been filled with defaults
  // note: the property array must be const, because non const operator[] doesn't work for vector<bool>!
  assert(std::as_const(properties_a).get_property<bool>("bools")[12] == true);
}

void test_constructors() {

  // Default constructor should have no properties
  Property_container<std::size_t> a{};
  assert(a.num_properties() == 0);

  // Copy constructor should duplicate all properties
  auto& a_ints = a.add_property("ints", 0);
  auto& a_floats = a.add_property("floats", 0.0f);
  a.emplace_group(10);
  a.get_property<int>("ints")[3] = 1;
  a.get_property<float>("floats")[3] = 1.0f;
  Property_container<std::size_t> b{a};
  assert(b.num_properties() == a.num_properties() && b.num_properties() == 2);
  assert(b.get_property<int>("ints")[3] == a.get_property<int>("ints")[3] && b.get_property<int>("ints")[3] == 1);
  assert(b.get_property<float>("floats")[3] == a.get_property<float>("floats")[3] && b.get_property<float>("floats")[3] == 1.0f);

  // Copy-assignment operator should do effectively the same thing as the copy constructor
  Property_container<std::size_t> c;
  c = a;
  assert(c.num_properties() == a.num_properties() && c.num_properties() == 2);
  assert(c.get_property<int>("ints")[3] == a.get_property<int>("ints")[3] && c.get_property<int>("ints")[3] == 1);
  assert(c.get_property<float>("floats")[3] == a.get_property<float>("floats")[3] && c.get_property<float>("floats")[3] == 1.0f);

  // Copied property containers should not be synced with the original
  a.add_property("more_ints", 2);
  assert(a.num_properties() == 3);
  assert(b.num_properties() == 2);
  assert(c.num_properties() == 2);
  a.get_property<int>("ints")[4] = 2;
  assert(a.get_property<int>("ints")[4] == 2);
  assert(b.get_property<int>("ints")[4] == 0);
  assert(c.get_property<int>("ints")[4] == 0);

  // Copy assignment should not invalidate previously obtained array references,
  // but it should update their values
  auto &b_ints = b.get_property<int>("ints");
  auto &b_floats = b.get_property<float>("floats");
  assert(b_ints[4] == 0);
  b = a;
  assert(b.num_properties() == 3);
  assert(b_ints[4] == 2);

  // Move assignment shouldn't invalidate references either
  Property_container<std::size_t> d{c};
  auto &d_ints = d.get_property<int>("ints");
  assert(d_ints[4] == 0);
  d = std::move(a);
  assert(d.num_properties() == 3);
  assert(d_ints[4] == 2);

  // Moved-from should be empty
  // All properties are preserved, though
  assert(a.num_properties() == 3);
  assert(a.size() == 0);
  assert(a_ints.capacity() == 0);
  assert(a_floats.capacity() == 0);

  // Move constructor should behave like move assignment
  Property_container<std::size_t> e{std::move(b)};
  assert(e.num_properties() == 3);
  assert(b.num_properties() == 3);
  assert(b.size() == 0);
  assert(b_ints.capacity() == 0);
  assert(b_floats.capacity() == 0);
}


int main() {

  test_property_creation();
  test_element_access();
  test_emplace_group();
  test_append();
  test_constructors();

  return 0;
}
