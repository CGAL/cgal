#ifndef CGAL_LLOYD_OPTIMIZE_MESH_2_H
#define CGAL_LLOYD_OPTIMIZE_MESH_2_H

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/Constrained_voronoi_diagram_2.h>

namespace CGAL
{

  template<typename CDT>
  void lloyd_optimize_mesh_2(CDT& cdt,
                             const unsigned int nb_iterations)
  {
    typedef typename CDT::Geom_traits Gt;
    typedef CGAL::Lipschitz_sizing_field_2<Gt> Lip_sizing;

    std::set<typename Gt::Point_2> points;
    for(typename CDT::Finite_vertices_iterator
      vit = cdt.finite_vertices_begin();
      vit != cdt.finite_vertices_end();
      ++vit)
        points.insert(vit->point());
    Lip_sizing size(points.begin(), points.end());
    
    typedef CGAL::Mesh_2::Lloyd_move_2<CDT, Lip_sizing> Mf;
    CGAL::Mesh_2::Mesh_global_optimizer_2<CDT, Mf> lloyd(cdt);
    lloyd.set_sizing_field(size);

    CGAL::Constrained_voronoi_diagram_2<CDT> cvd(&cdt);
    for(unsigned int i = 0; i < nb_iterations; ++i)
    {
      cvd.tag_faces_blind();
      lloyd(1/*nb_iterations*/);
    }
    cvd.tag_faces_blind();//update blindness

    //update inside/outside info
    typedef CGAL::Lipschitz_sizing_field_criteria_2<CDT, Lip_sizing>
      Lip_criteria;
    CGAL::Delaunay_mesher_2<CDT, Lip_criteria> mesher(cdt);
    mesher.mark_facets();
  }

} //end namespace CGAL

#endif
