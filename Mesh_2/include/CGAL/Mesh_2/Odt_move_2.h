#ifndef CGAL_ODT_MOVE_2_H
#define CGAL_ODT_MOVE_2_H

#include <CGAL/Mesh_2/Uniform_sizing_field_2.h>

namespace CGAL
{
namespace Mesh_2
{
  template<typename CDT,
           typename SizingField = Uniform_sizing_field<CDT> >
  class Odt_move_2
  {
    typedef typename CDT::Vertex_handle          Vertex_handle;
    typedef typename CDT::Geom_traits::Point_2   Point_2;

  public:
    typedef SizingField Sizing_field;

  public:
    Point_2 operator()(const CDT& cdt, Vertex_handle v) const
    {
      //todo
      return CGAL::ORIGIN;
    }
  };

} //end namespace Mesh_2
} //end namespace CGAL

#endif //CGAL_ODT_MOVE_2_H
