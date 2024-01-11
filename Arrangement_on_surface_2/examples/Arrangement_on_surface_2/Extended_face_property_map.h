#ifndef EXTENDED_FACE_PROPERTY_MAP_H
#define EXTENDED_FACE_PROPERTY_MAP_H

// A property map that reads/writes the information to/from the extended face.
template <typename Arrangement, class Type> class Extended_face_property_map {
public:
  using Face_handle = typename Arrangement::Face_handle;

  // Boost property type definitions.
  using category = boost::read_write_property_map_tag;
  using value_type = Type;
  using reference = value_type&;
  using key_type = Face_handle;

  // The get function is required by the property map concept.
  friend reference get(const Extended_face_property_map& /* map */, key_type key)
  { return key->data(); }

  // The put function is required by the property map concept.
  friend void put(const Extended_face_property_map& /* map */,
                  key_type key, value_type val)
  { key->set_data(val); }
};

#endif
