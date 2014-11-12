#ifndef CGAL_LLOYD_MOVE_2_H
#define CGAL_LLOYD_MOVE_2_H

#include <CGAL/Mesh_2/Uniform_sizing_field_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_voronoi_diagram_2.h>

namespace CGAL
{
namespace Mesh_2
{
  template<typename CDT,
           typename SizingField = Uniform_sizing_field<typename CDT::Geom_traits> >
  class Lloyd_move_2
  {
    typedef typename CDT::Vertex_handle          Vertex_handle;
    typedef typename CDT::Geom_traits::Point_2   Point_2;
    typedef typename CDT::Cvd_cell               Cvd_cell;

    typedef typename CDT::Geom_traits::Segment_2  Segment;
    typedef typename CDT::Geom_traits::Triangle_2 Triangle;
    typedef typename CDT::Geom_traits::FT         FT;
    typedef typename CDT::Geom_traits::Vector_2   Vector_2;

  public:
    typedef SizingField Sizing_field;

  public:
    Vector_2 operator()(Vertex_handle v,
          const CDT& cdt,
          const Sizing_field& sizing_field = Sizing_field()) const
    {
      if(cdt.are_there_incident_constraints(v))
        return CGAL::NULL_VECTOR;

      Point_2 p = v->point();
      Vector_2 move = CGAL::NULL_VECTOR;
      FT sum_masses(0);

      Cvd_cell cell = cdt.dual(v);
      if(cell.is_infinite() || cell.is_empty())
        return CGAL::NULL_VECTOR; //center of mass is at infinity!

      CGAL_assertion(cell.number_of_vertices() > 2);

      typename Cvd_cell::segment_iterator sit = cell.segments_begin();
      for( ; sit != cell.segments_end(); ++sit)
      {
        Segment s = *sit;
        Triangle tri(p, s.source(), s.target());
        Point_2 tri_centroid = CGAL::centroid(tri);

        // Compute mass
        FT density = density_2d(tri_centroid, sizing_field);
        FT abs_area = CGAL::abs(tri.area());
        FT mass = abs_area * density;

        move = move + mass * Vector_2(p, tri_centroid);
        sum_masses += mass;
      }

      CGAL_assertion(sum_masses != 0.0);
      return move / sum_masses;
    }

  public:
    void before_move(CDT& cdt)
    {
      update_blind_faces(cdt);
    }

    void after_move(CDT& cdt)
    {
      //update inside/outside tags
      typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt);
      mesher.mark_facets();
    }

    void after_all_moves(CDT& cdt)
    {
      update_blind_faces(cdt);
    }

  private:
    FT density_2d(const Point_2& p,
                  const Sizing_field& sizing_field) const
    {
      FT s = sizing_field(p);
      CGAL_assertion( s > 0. );

      // 1 / s^(d+2)
      return ( 1/(s*s*s*s) );
    }

    void update_blind_faces(CDT& cdt)
    {
      //update blindness
      CGAL::Constrained_voronoi_diagram_2<CDT> cvd(&cdt);
      cvd.tag_faces_blind();
    }

#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  public:
    static std::string name() { return std::string("Lloyd"); }
#endif

  };

} //end namespace Mesh_2
} //end namespace CGAL

#endif //CGAL_LLOYD_MOVE_2_H
