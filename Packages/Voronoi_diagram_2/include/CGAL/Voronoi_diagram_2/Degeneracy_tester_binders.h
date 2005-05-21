#ifndef CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H
#define CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=======================================================================

template<class VDA>
class Edge_degeneracy_tester_binder
{
 private:
  typedef typename VDA::Edge_degeneracy_tester  EDT;

 public:
  typedef typename EDT::result_type               result_type;
  typedef Arity_tag<1>                            Arity;

  Edge_degeneracy_tester_binder(const VDA* vda = NULL) : vda_(vda) {}

  template<class A>
  bool operator()(const A& a) const {
    CGAL_precondition( vda_ != NULL );
    return vda_->edge_tester()(vda_->dual(), a);
  }

 private:
  const VDA* vda_;
};

//=======================================================================

template<class VDA>
class Face_degeneracy_tester_binder
{
 private:
  typedef typename VDA::Face_degeneracy_tester  FDT;

 public:
  typedef typename FDT::result_type               result_type;
  typedef Arity_tag<1>                            Arity;

  Face_degeneracy_tester_binder(const VDA* vda = NULL) : vda_(vda) {}

  template<class A>
  bool operator()(const A& a) const {
    CGAL_precondition( vda_ != NULL );
    return vda_->face_tester()(vda_->dual(), vda_->edge_tester(), a);
  }

 private:
  const VDA* vda_;
};

//=======================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_DEGENERACY_TESTER_BINDERS_H
