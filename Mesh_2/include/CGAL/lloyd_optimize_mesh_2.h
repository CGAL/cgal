#ifndef CGAL_LLOYD_OPTIMIZE_MESH_2_H
#define CGAL_LLOYD_OPTIMIZE_MESH_2_H

#include <CGAL/Mesh_2/Mesh_global_optimizer_2.h>
#include <CGAL/Mesh_2/Lloyd_move_2.h>
#include <CGAL/Mesh_2/Lipschitz_sizing_field_2.h>

namespace CGAL
{

  template<typename CDT>
  void lloyd_optimize_mesh_2(CDT& cdt,
                             const unsigned int nb_iterations = 1,
                             const double convergence_ratio = 0.001)
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
    CGAL::Mesh_2::Mesh_global_optimizer_2<CDT, Mf> lloyd(cdt,
                                                         convergence_ratio);
    lloyd.set_sizing_field(size);
    lloyd(nb_iterations);
  }

} //end namespace CGAL

#endif
