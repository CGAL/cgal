#include <CGAL/convex_hull_3.h> 
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/Nary_union.h>
#include <CGAL/Nef_3/Relabel_volume.h>
#include <CGAL/Nef_3/Mark_bounded_volumes.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>
#include <CGAL/Nef_3/Bounding_box_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

namespace CGAL {

template <typename pIt, class InpIt, class ForIt, class OutIt, class AdBiFct>
  OutIt fold_indices_polyhedron( const pIt points, InpIt first1, InpIt beyond1,
				 ForIt first2, ForIt beyond2,
				 OutIt result,
				 AdBiFct fct) {
  for ( ; first1 != beyond1; ++first1) {
    for ( ForIt i = first2; i != beyond2; ++i) {
      *result++ = fct( *(points+(*first1)), i->point());
    }
  }
  return result;
}

struct Add_points {
    template <class Point>
    Point operator()( const Point& p, const Point& q) const {
        using CGAL::ORIGIN;
        return ORIGIN + (p-ORIGIN) + (q-ORIGIN);
    }
};

template<typename Nef_polyhedron>
class Trunk_offset {
  
  typedef typename Nef_polyhedron::Kernel Kernel;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;
  typedef Polyhedron_3<Kernel> Polyhedron;
  typedef typename Polyhedron::Vertex_const_iterator PVertex_const_iterator;
  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Polyhedron::Halfedge_around_facet_const_circulator
    Halfedge_around_facet_const_circulator;
  typedef typename Nef_polyhedron::Object_handle Object_handle;
  typedef typename Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
  typedef typename Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Nef_polyhedron::Halffacet_const_iterator Halffacet_const_iterator;
  typedef typename Nef_polyhedron::Halffacet_const_handle Halffacet_const_handle;
  typedef typename Nef_polyhedron::SHalfedge_around_facet_const_circulator 
    SHalfedge_around_facet_const_circulator;
  typedef typename Box_intersection_d::Box_d<double,3> Box;
  typedef Bounding_box_3<typename Kernel::Kernel_tag, Kernel> EBox;
  typedef Nary_union<Nef_polyhedron> NUBQ;
  typedef Relabel_volume<Nef_polyhedron> Relabel_volume;
  typedef Mark_bounded_volumes<Nef_polyhedron> Mark_bounded_volumes;

  
  bool has_two_unmarked_volumes(Nef_polyhedron& N) {
    std::cerr << "number of volumes " << N.number_of_volumes() << std::endl;
    int argc = 0;
    char* argv[1];
    QApplication a(argc, argv);
    CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
      new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
    a.setMainWidget(w);
    w->show();
    a.exec();

    Volume_const_iterator ci=++N.volumes_begin();
    while(ci != N.volumes_end() && ci->mark()) ++ci;
    return ci != N.volumes_end();
  }

 public:
  Trunk_offset() {}

  template<typename p_it, typename f_it>
    Nef_polyhedron operator()(p_it pbegin, p_it pend,
			      f_it fbegin, f_it fend,
			      const Polyhedron& P) {
    
    Polyhedron Ptemp;
    Add_points add;
    NUBQ nary_union;

    for(f_it f = fbegin; f != fend; ++f) {
      std::vector<Point> points;
      fold_indices_polyhedron( pbegin,
			       f->first, f->second,
			       P.vertices_begin(), P.vertices_end(),
			       back_inserter( points),
			       add);
      convex_hull_3( points.begin(), points.end(), Ptemp); 
      nary_union.add_polyhedron(Nef_polyhedron(Ptemp));
    }
    
    Nef_polyhedron result = nary_union.get_union();
    CGAL_assertion_msg(result.number_of_volumes() == 3,
		       "not implemented, yet");
    Relabel_volume rv;
    result.delegate(rv);
    CGAL_assertion_msg(result.number_of_volumes() == 2,
		       "not implemented, yet");
    result = (!result).closure();
    return result;
  }
};

} //namespace CGAL
