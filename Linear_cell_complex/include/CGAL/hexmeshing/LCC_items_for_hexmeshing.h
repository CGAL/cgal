
#ifndef LCC_ITEMS_FOR_HEXMESHING_H
#define LCC_ITEMS_FOR_HEXMESHING_H

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/hexmeshing/Hexmeshing_generic_point.h>
#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>

namespace CGAL::internal::Hexmeshing {
/**
 * @brief Enumeration representing the type of volume in the hexahedral mesh
 * 
 * This enum class defines the different states a volume can be in during
 * the hexahedral mesh generation and refinement process.
 */
enum class VolumeType {
  NONE,         // Newly created vol_attribute
  REFINEMENT,   // Previously NONE volumes that was selected for refinement
  ID_EXPANSION, // Expansion of the Identified layer of volumes
  IDENTIFIED    // Volumes identified for refinement
};

class LCCItemsForHexmeshing
{
public:
  /**
   * @brief A wrapper class for Linear Cell Complex dart attributes
   * @tparam Storage The storage type for the Linear Cell Complex
   * 
   * This class defines the attributes associated with vertices, faces, and volumes
   * in the hexahedral mesh.
   */
  template < class Storage >
  struct Dart_wrapper
  {
    /**
     * @brief Attribute class for vertices with point information
     * 
     * Extends CGAL::Cell_attribute_with_point to store vertex coordinates and IDs.
     * Each vertex can be identified by a unique ID, which can be a simple number
     * or a composite of up to three numbers.
     */
    struct VertexAttr : public Cell_attribute_with_point<Storage> {
      using Base = Cell_attribute_with_point<Storage>;
      VertexAttr():
        Base(),
        id(max_id)
      {}

      VertexAttr(const Point& p):
        Base(p),
        id(max_id)
      {}
      // The id can be a simple number, or a pair up to three numbers
      // Pair identifiers are placeholders until they can be identified with a single size_t number

      static constexpr size_t max_id = std::numeric_limits<size_t>::max();
      size_t id = max_id;

      Vector normal;
    };

    /**
     * @brief Attribute values for volumes in the hexahedral mesh
     * 
     * Stores information about the state and properties of a volume cell,
     * including its iteration state, type, area ID, ownership, and connected component ID.
     */
    struct VolumeAttrValue {
      static constexpr uint max_cc_id = std::numeric_limits<uint>::max();

      char iteration = -1;
      VolumeType type = VolumeType::NONE;

      AreaId area_id = {0,0,0}; // Used for ghost cells;
      bool owned = true; // Disable this when creating others constraint area
      size_t cc_id = max_cc_id;

      Point centroid;
      double fraction;
      Vector gradient;
    };

    /**
     * @brief Attribute values for faces in the hexahedral mesh
     * 
     * Stores information about face properties including template ID,
     * plane orientation, plane ID, and connected component ID.
     */
    struct FaceAttrValue {
      static constexpr uint max_plane_id = std::numeric_limits<uint>::max();
      static constexpr uint max_cc_id = std::numeric_limits<uint>::max();

      char template_id = 0;
      std::bitset<3> plane;
      size_t plane_id = max_plane_id;
      size_t cc_id = max_cc_id;

      Segment dual_edge;
      Point intersection;
      Vector normal;
    };

    typedef Cell_attribute<Storage, FaceAttrValue> FaceAttr;
    typedef Cell_attribute<Storage, VolumeAttrValue> VolumeAttr;
    typedef std::tuple<VertexAttr, void, FaceAttr, VolumeAttr> Attributes;
  };
};

using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItemsForHexmeshing>;
using Dart_handle = typename LCC::Dart_handle;
using Vertex_handle = typename LCC::Vertex_attribute_handle;
using size_type = typename LCC::size_type;

using DartInfo = LCCItemsForHexmeshing::Dart_wrapper<LCC::Storage>;

// Do not use this, only to satisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::VolumeAttrValue& first, const DartInfo::VolumeAttrValue& second){
  return true;
}

// Do no use this, only to statisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::FaceAttrValue& first, const DartInfo::FaceAttrValue& second){
  return true;
}

}




#endif