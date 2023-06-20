// STL includes.
#include <unordered_map>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include "include/utils.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>

// Custom Neighbor_query and Region_type classes for region growing.
namespace Custom {

  // An object that stores indices of all its neighbors.
  struct Object {
    std::vector<std::vector<Object>::const_iterator> neighbors;
    bool operator==(const Object& obj) const {
      return neighbors == obj.neighbors;
    }
  };

  // A range of objects.
  using Objects = std::vector<Object>; // std::list<Object> works as well

  // The Neighbor_query functor that accesses neighbors stored in
  // the object struct above.
  class Neighbor_query {
    const Objects& m_objects;

    using Item = typename Objects::const_iterator;
    using Region = std::vector<Item>;

  public:
    Neighbor_query(const Objects& objects) :
    m_objects(objects)
    { }

    void operator()(
      const Item &query,
      std::vector<Item>& neighbors) const {

      for (auto it = m_objects.begin(); it != m_objects.end(); it++) {
        if (it == query) {
          neighbors = query->neighbors;
          return;
        }
      }
    }
  };

  // The Region_type class, where the function is_part_of_region() verifies
  // a very specific condition that the first and second objects in the
  // range are in fact neighbors; is_valid_region() function always
  // returns true after the first call to the update() function.
  // These are the only functions that have to be defined.
  class Region_type {
    bool m_is_valid = false;
    const std::vector<Object>& m_input;

  public:
    Region_type(const std::vector<Object> &input) : m_input(input) { }

    using Primitive = std::size_t;

    using Item = std::vector<Object>::const_iterator;
    using Region = std::vector<Item>;

    using Region_unordered_map = std::unordered_map<Item, std::size_t, CGAL::Shape_detection::internal::hash_item<Item> >;
    using Region_index_map = boost::associative_property_map<Region_unordered_map>;

    Region_index_map region_index_map() {
      Region_index_map index_map(m_region_map);
      return index_map;
    }

    bool is_part_of_region(
      const Item query,
      const Region& region) const {

      if (region.size() == 0) return false;

      auto it = m_input.begin();

      if (query == it || query == (it + 1))
        return true;

      return false;
    }

    inline bool is_valid_region(const Region&) const {
      return m_is_valid;
    }

    Primitive primitive() {
      return Primitive();
    }

    bool update(const Region&) {
      m_is_valid = true;
      return m_is_valid;
    }

  private:
    Region_unordered_map m_region_map;
  };

} // namespace Custom

// Typedefs.
using Object         = Custom::Object;
using Objects        = Custom::Objects;
using Neighbor_query = Custom::Neighbor_query;
using Region_type    = Custom::Region_type;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

int main() {

  // Define a range of objects, where the first two objects form
  // the first region, while the third object forms the second region.
  // Note that Objects is a random access container here, however the
  // same algorithm/example can work with other containers, e.g. std::list.
  Objects objects(3);
  auto it = objects.begin();

  // Region 1.
  objects[0].neighbors.push_back(it+1);
  objects[1].neighbors.push_back(it);

  // The single third object constitutes the second region.

  std::cout << "* number of input objects: " << objects.size() << std::endl;
  assert(objects.size() == 3);

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query = Neighbor_query(objects);
  Region_type    region_type    = Region_type(objects);

  // Create an instance of the region growing class.
  Region_growing region_growing(
    objects, neighbor_query, region_type);

  // Run the algorithm.
  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));
  std::cout << "* number of found regions: " << regions.size() << std::endl;
  assert(regions.size() == 2);
  return EXIT_SUCCESS;
}
