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
#ifdef DEBUG
#ifndef NO_DOUBLE_MAP_DEBUG
	std::cerr << "get_next_element(" << &*(m.front()->second) << ")\n";
#endif
#endif
        return m.front()->second;

      }

      void add_element(const Element& e, const Quality& q)
      {
#ifdef DEBUG
#ifndef NO_DOUBLE_MAP_DEBUG
	std::cerr << "add_element(" << &*e << ")\n";
#endif
#endif
        m.insert(e, q);
      }

      void remove_next_element()
      {
#ifdef DEBUG
#ifndef NO_DOUBLE_MAP_DEBUG
	std::cerr << "pop_front(" << &*(m.front()->second) << ")\n";
#endif
#endif
        m.pop_front();
      }

      void remove_element(Element& e)
      {
#ifdef DEBUG
#ifndef NO_DOUBLE_MAP_DEBUG
	std::cerr << "remove_element(" << &*e << ")\n";
#endif
#endif
        m.erase(e);
      }
    }; // end Double_map_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H
