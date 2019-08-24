#ifndef CGAL_FACEWIDTH_H
#define CGAL_FACEWIDTH_H

#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_topology/internal/Generic_map_selector.h>
#include <CGAL/Surface_mesh_topology/internal/Shortest_noncontractible_cycle.h>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

template <class Mesh_>
class Facewidth {

public:

  using Mesh_original = Mesh_;

  struct Default_weight_functor {
    using Weight_t = unsigned int;
    template <class T>
    Weight_t operator() (T) const { return 1; }
  };

  using Gmap                       = typename Generic_map_selector<Mesh_original>::Generic_map;
  using Gmap_wrapper               = Generic_map_selector<Mesh_original>;
  using Dart_handle_original       = typename Gmap_wrapper::Dart_handle_original;
  using Dart_handle                = typename Gmap::Dart_handle;
  using Dart_const_handle          = typename Gmap::Dart_const_handle;
  using size_type                  = typename Gmap::size_type;
  using Path                       = CGAL::Surface_mesh_topology::Path_on_surface<Mesh_original>;
  using SNC                        = Shortest_noncontractible_cycle<Gmap>;
  
  Facewidth(Mesh_original& gmap)
  {
    typename Gmap_wrapper::Origin_to_copy_map origin_to_radial;
    Gmap_wrapper::copy(m_radial_map, gmap, origin_to_radial, m_copy_to_origin);
    m_copy_to_origin.clear();
    Gmap_wrapper::copy(m_gmap, gmap, m_origin_to_copy, m_copy_to_origin);
    // Initialize 0-attributes for m_gmap
    for (auto it = m_gmap.darts().begin(), itend = m_gmap.darts().end(); it != itend; ++it) {
      if (m_gmap.template attribute<0>(it)==NULL)
      { m_gmap.template set_attribute<0>(it, m_gmap.template create_attribute<0>()); }
    }
    // Assign values
    int counter = 0;
    for (auto it = m_gmap.template one_dart_per_cell<0>().begin(),
              itend = m_gmap.template one_dart_per_cell<0>().end(); it != itend; ++it)
    {
      m_vertex_list.push_back(it);
      m_gmap.template info<0>(it) = counter;
      ++counter;
    }

    // m_face_list contains dart handles of m_gmap
    for (auto it = m_gmap.template one_dart_per_cell<2>().begin(),
              itend = m_gmap.template one_dart_per_cell<2>().end(); it != itend; ++it)
    { m_face_list.push_back(it); }

    // Create edge_list
    std::vector<Dart_handle> edge_list;
    for (auto it = m_radial_map.template one_dart_per_cell<1>().begin(),
              itend = m_radial_map.template one_dart_per_cell<1>().end(); it != itend; ++it)
    {
      edge_list.push_back(it);
    }

    // m_radial_map hasn't been changed so far

    // face_list contains dart handles of m_radial_map
    std::vector<Dart_handle> face_list;
    for (auto it = m_radial_map.template one_dart_per_cell<2>().begin(),
              itend = m_radial_map.template one_dart_per_cell<2>().end(); it != itend; ++it)
    { face_list.push_back(it); }
    
    // Adding "centroids"
    std::vector<Dart_handle> centroids;
    for (auto it : face_list) // face_list contains dart handles of m_radial_map
    {
      auto new_vertex = m_radial_map.insert_cell_0_in_cell_2(it);
      centroids.push_back(new_vertex);
    }

    // Initialize 1-attributes of m_radial_map
    for (auto it = m_radial_map.darts().begin(), itend = m_radial_map.darts().end(); it != itend; ++it) {
      if (m_radial_map.template attribute<1>(it)==NULL)
      { m_radial_map.template set_attribute<1>(it, m_radial_map.template create_attribute<1>()); }
    }

    // Assign values
    for (int i = 0; i < centroids.size(); ++i) {
      auto u = centroids[i];
      bool first_run = true;
      for (Dart_handle it = u; first_run || it != u; it = m_radial_map.next(m_radial_map.opposite2(it))) {
        first_run = false;
        m_radial_map.template info<1>(it) = i;
      }
    }

    // Initialize 0-attributes of m_radial_map
    for (auto it = m_radial_map.darts().begin(), itend = m_radial_map.darts().end(); it != itend; ++it) {
      if (m_radial_map.template attribute<0>(it)==NULL)
      { m_radial_map.template set_attribute<0>(it, m_radial_map.template create_attribute<0>()); }
    }
    for (auto it = m_radial_map.darts().begin(), itend = m_radial_map.darts().end(); it != itend; ++it) {
      m_radial_map.template info<0>(it) = -1;
    }
    typedef typename Gmap::template Attribute_handle<0>::type Attribute_handle_0;
    for (Attribute_handle_0 att_it = m_gmap.template attributes<0>().begin(), 
                            att_itend = m_gmap.template attributes<0>().end(); att_it != att_itend; ++att_it)
    {
      auto it_radial = origin_to_radial[m_copy_to_origin[att_it->dart()]];
      m_radial_map.template info<0>(it_radial) = att_it->info();
    }

    // Remove the marked edges of m_radial_map
    for (auto dh : edge_list) 
    {
      CGAL_assertion(m_radial_map.template is_removable<1>(dh));
      m_radial_map.template remove_cell<1>(dh);
    }
  }

  ~Facewidth()
  {
  }
  
  std::vector<Dart_handle_original> compute_facewidth()
  {
    m_cycle.clear();
    // Find edgewidth of the radial map
    SNC snc_to_find_facewidth(m_radial_map);
    Path_on_surface<Gmap> edgewidth_of_radial_map = snc_to_find_facewidth.compute_edgewidth();

    int last_vertex_index = -1, first_vertex_index = -1;
    int last_face_index = -1;
    for (int i = 0, n = edgewidth_of_radial_map.length(); i <= n; i++)
    {
      Dart_const_handle dh = edgewidth_of_radial_map[i % n];
      int face_index = m_radial_map.template info<1>(dh);
      if (m_radial_map.template info<0>(dh) == -1) dh = m_radial_map.next(dh);
      CGAL_assertion(m_radial_map.template info<0>(dh) != -1);
      int vertex_index = m_radial_map.template info<0>(dh);
      
      if (last_face_index == face_index)
      {
        m_cycle.push_back(m_copy_to_origin[m_vertex_list[last_vertex_index]]);  
        m_cycle.push_back(m_copy_to_origin[m_face_list[face_index]]); 
      }

      // if (first_vertex_index == -1) first_vertex_index = vertex_index;
      last_vertex_index = vertex_index;
      last_face_index = face_index;
    }
    return m_cycle;
  }

private:
  Gmap m_gmap, m_radial_map;
  std::vector<Dart_handle> m_vertex_list, m_face_list;
  std::vector<Dart_handle_original> m_cycle;
  typename Gmap_wrapper::Origin_to_copy_map m_origin_to_copy;
  typename Gmap_wrapper::Copy_to_origin_map m_copy_to_origin;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif
