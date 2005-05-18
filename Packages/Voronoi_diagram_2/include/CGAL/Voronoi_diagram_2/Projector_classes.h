#ifndef CGAL_VORONOI_DIAGRAM_2_PROJECTOR_CLASSES_H
#define CGAL_VORONOI_DIAGRAM_2_PROJECTOR_CLASSES_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//------------------------------------------------------------------------

template<class DG>
struct Tds_project {
  typedef DG   Dual_graph;

  typedef typename DG::Triangulation_data_structure
  Dual_graph_data_structure;

  const Dual_graph_data_structure& operator()(const Dual_graph& dg) const
  {
    return dg.tds();
  }

  const Dual_graph_data_structure& operator()(const Dual_graph* dg) const
  {
    return dg->tds();
  }
};

//------------------------------------------------------------------------

template<class DG>
struct Ds_project {
  typedef DG   Dual_graph;

  typedef typename DG::Data_structure
  Dual_graph_data_structure;

  const Dual_graph_data_structure& operator()(const Dual_graph& dg) const
  {
    return dg.data_structure();
  }

  const Dual_graph_data_structure& operator()(const Dual_graph* dg) const
  {
    return dg->data_structure();
  }
};

//------------------------------------------------------------------------

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_PROJECTOR_CLASSES_H
