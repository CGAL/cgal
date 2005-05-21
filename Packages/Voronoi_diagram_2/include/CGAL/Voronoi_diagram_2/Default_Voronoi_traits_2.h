#ifndef CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
#define CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Cached_degeneracy_testers.h>
#include <CGAL/Handle_for_virtual.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class DG, class ET, class FT>
class Default_Voronoi_traits_2
{
 private:
  typedef Default_Voronoi_traits_2<DG,ET,FT>   Self;

 public:
  typedef DG  Dual_graph;
  typedef ET  Edge_degeneracy_tester;
  typedef FT  Face_degeneracy_tester;

  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================

template<class DG, class ET, class FT>
class Default_cached_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_cached_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

  typedef Cached_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester_base,
					Edge_degeneracy_tester>
  Face_degeneracy_tester;


 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};


//=========================================================================
//=========================================================================


template<class DG, class ET, class FT>
class Default_ref_counted_Voronoi_traits_2
{
 private:
  typedef ET  Edge_degeneracy_tester_base;
  typedef FT  Face_degeneracy_tester_base;

  typedef Default_ref_counted_Voronoi_traits_2<DG,ET,FT>  Self;

 public:
  typedef DG           Dual_graph;

  typedef Ref_counted_edge_degeneracy_tester<Edge_degeneracy_tester_base>
  Edge_degeneracy_tester;

  typedef Ref_counted_face_degeneracy_tester<Face_degeneracy_tester_base,
					     Edge_degeneracy_tester>
  Face_degeneracy_tester;

 public:
  const Edge_degeneracy_tester& edge_degeneracy_tester_object() const {
    return e_tester_;
  }

  const Face_degeneracy_tester& face_degeneracy_tester_object() const {
    return f_tester_;
  }

 protected:
  Edge_degeneracy_tester e_tester_;
  Face_degeneracy_tester f_tester_;
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_VORONOI_TRAITS_2_H
