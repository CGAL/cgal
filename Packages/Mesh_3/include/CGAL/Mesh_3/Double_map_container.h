#ifndef CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H
#define CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H

#include <set>
#include <CGAL/Double_map.h>

namespace CGAL {

  namespace Mesh_3 {

    /** Example of a \c NON-REMOVABLE container:
        one cannot remove elements from it, but the front. */
    template <typename Element, class Quality>
    class Double_map_container 
    {
      // --- private datas ---
      Double_map<Element, Quality> m;

    public:
      bool empty() const
      {
        return m.empty();
      }

      Element get_next_element()
      {
        CGAL_assertion(!m.empty());
        
        return m.front()->first;

      }

      void add_element(const Element& e, const Quality& q)
      {
        m.insert(e, q);
      }

      void remove_next_element()
      {
        m.pop_front();
      }

      void remove_element(const Element& e)
      {
        m.erase(e);
      }
    }; // end Double_map_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H
