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
    typedef typename CDT::Geom_traits::Vector_2  Vector_2;
    typedef typename CDT::Geom_traits::Segment_2 Segment;

  public:
    typedef SizingField Sizing_field;

  public:
    Vector_2 operator()(Vertex_handle v,
          const CDT& cdt,
          const Sizing_field& m_sizing_field = Sizing_field()) const
    {
      // Calculate the average of circumcenters incident to current vertex.
      CGAL_assertion(!cdt.is_infinite(v));
      //v should not be constrained, o.w. it is not allowed to move
      CGAL_assertion(!cdt.are_there_incident_constraints(v));

      typename CDT::Cvd_cell cell = cdt.dual(v);
      CGAL_assertion(!cell.is_infinite());
      CGAL_assertion(cell.is_simply_ccw_oriented());

      if(cell.is_empty())
        return CGAL::NULL_VECTOR;

      double area = 0.;
      Vector_2 move = CGAL::NULL_VECTOR;
      for(std::size_t i = 0; i < cell.number_of_vertices(); ++i)
      {
        move = move + Vector_2(CGAL::ORIGIN, cell.point(i));
      }

      return move;
    }

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("ODT"); }
#endif
  };

} //end namespace Mesh_2
} //end namespace CGAL

#endif //CGAL_ODT_MOVE_2_H
