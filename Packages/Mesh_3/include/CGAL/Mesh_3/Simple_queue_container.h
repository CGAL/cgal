#ifndef CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
#define CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H

#include <queue>

namespace CGAL {

  namespace Mesh_3 {

    /** Example of a \c NON-REMOVABLE container:
        one cannot remove elements from it, but the front. */
    template <typename Element>
    class Simple_queue_container 
    {
      // --- private datas ---
      std::queue<Element> q;

    public:
      bool empty() const
      {
        return q.empty();
      }

      Element& get_next_element()
      {
        return q.front();
      }

      void add_element(const Element& e)
      {
        q.push(e);
      }

      void remove_next_element()
      {
        q.pop();
      }
    }; // end Simple_queue_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
