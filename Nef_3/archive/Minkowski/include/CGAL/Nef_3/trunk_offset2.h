#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/Nary_union_by_queue.h>
#include <CGAL/Nef_3/Nary_union_by_pq.h>
//#include <CGAL/Nef_3/Nary_union_by_pq_me.h>
#include <CGAL/Nef_3/Relabel_volume.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>
#include <vector>

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

template<typename Trunk_box>
struct Callback {
  bool intersection;
  Callback() : intersection(false) {}
  void operator()( Trunk_box& box0, Trunk_box& box1 ) {intersection=true;}
  bool intersect() const {return intersection;}
};


template<typename Kernel>
class Trunk_offset {

  typedef typename Kernel::Point_3 Point;
  typedef Polyhedron_3<Kernel> Polyhedron;
  typedef Nef_polyhedron_3<Kernel> Nef_polyhedron;
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
  typedef Nary_union_by_queue<Nef_polyhedron> NUBQ;
  typedef Relabel_volume<Nef_polyhedron> Relabel_volume;

  class Trunk_box : public Box_intersection_d::Box_d< double, 3 > {

    typedef std::pair<double, double> double_pair;
    typedef Box_intersection_d::box_limits<double> box_limits;

    void extend( const Point& p) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( p.x() );
      q[1] = CGAL::to_interval( p.y() );
      q[2] = CGAL::to_interval( p.z() );
      Box_intersection_d::Box_d< double, 3 >::extend(q);
    }

  public:
   Trunk_box(Facet_const_handle f) {
      Halfedge_around_facet_const_circulator e(f->facet_begin()), end(e);
      CGAL_For_all(e,end)
        extend(e->vertex()->point());
    }

   Trunk_box(Halffacet_const_iterator f) {
      SHalfedge_around_facet_const_circulator e(f->facet_cycles_begin()), end(e);
      CGAL_For_all(e,end)
        extend(e->source()->source()->point());
    }
  };

  int mod, suf;

 public:
  Trunk_offset(int modulo=1, int suffix=0) : mod(modulo), suf(suffix) {}

  template<typename p_it, typename f_it>
    Nef_polyhedron operator()(p_it pbegin, p_it pend,
                              f_it fbegin, f_it fend,
                              const Polyhedron& P) {

    Box pbox;
    for(PVertex_const_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( v->point().x() );
      q[1] = CGAL::to_interval( v->point().y() );
      q[2] = CGAL::to_interval( v->point().z() );
      pbox.extend(q);
    }

    Box nbox;
    for(p_it pit=pbegin; pit!=pend(); ++pit) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( pit->x() );
      q[1] = CGAL::to_interval( pit->y() );
      q[2] = CGAL::to_interval( pit->z() );
      nbox.extend(q);
    }

    int dx = pbox.max_coord(0) - pbox.min_coord(0);
    int dy = pbox.max_coord(1) - pbox.min_coord(1);
    int dz = pbox.max_coord(2) - pbox.min_coord(2);

    int nx = (nbox.max_coord(0) - nbox.min_coord(0)) / (dx * 2);
    int ny = (nbox.max_coord(1) - nbox.min_coord(1)) / (dy * 2);
    int nz = (nbox.max_coord(2) - nbox.min_coord(2)) / (dz * 2);

    std::cerr << dx << " " << dy << " " << dz << std::endl;

    return Nef_polyhedron();

    Polyhedron Ptemp;
    Add_points add;
    NUBQ nubq;

    while(suf > 0) {
      --fend;
      --suf;
      std::vector<Point> points;
      fold_indices_polyhedron( pbegin,
                               fend->first, fend->second,
                               P.vertices_begin(), P.vertices_end(),
                               back_inserter( points),
                               add);
      convex_hull_3( points.begin(), points.end(), Ptemp);
      nubq.add_polyhedron(Nef_polyhedron(Ptemp));
    }

    {
      int i=0;
      for(f_it f = fbegin; f != fend; ++f, ++i) {
        if(i%mod!=0) continue;
        //        std::cerr << i << std::endl;
        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);
        convex_hull_3( points.begin(), points.end(), Ptemp);
        nubq.add_polyhedron(Nef_polyhedron(Ptemp));
      }
    }

    Nef_polyhedron result = nubq.get_union();

    while(result.number_of_volumes() != 3 && mod != 1) {
      //    while(mod != 1) {
      mod /= 2;
      int j = 0;
      for(f_it f = fbegin; f != fend; ++f, ++j) {
        if(j%mod!=0 || (j%(2*mod))==0) continue;
        //        std::cerr << j << std::endl;
        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);
        convex_hull_3( points.begin(), points.end(), Ptemp);
        nubq.add_polyhedron(Nef_polyhedron(Ptemp));
      }
      result = nubq.get_union();
    }

    Relabel_volume rv(true);
    result.delegate(rv);

    if(mod == 1)
      return result;

    NUBQ nubq2;
    nubq2.add_polyhedron(result);

    Box box;
    for(Vertex_const_iterator nv = result.vertices_begin(); nv != result.vertices_end(); ++nv) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( nv->point().x() );
      q[1] = CGAL::to_interval( nv->point().y() );
      q[2] = CGAL::to_interval( nv->point().z() );
      box.extend(q);
    }

    int k=0;
    for(f_it f = fbegin; f != fend; ++f, ++k) {
      if(k%(mod)==0) continue;
      //      std::cerr << k << std::endl;
      std::vector<Point> points;
      fold_indices_polyhedron( pbegin,
                               f->first, f->second,
                               P.vertices_begin(), P.vertices_end(),
                               back_inserter( points),
                               add);
      bool added = false;
      typedef typename std::vector<Point>::const_iterator p_it2;
      for(p_it2 pit = points.begin(); pit != points.end(); ++pit) {
        std::pair<double, double> q = CGAL::to_interval(pit->x());
        if(q.second < box.min_coord(0) || q.first > box.max_coord(0)) continue;
        q = CGAL::to_interval(pit->y());
        if(q.second < box.min_coord(1) || q.first > box.max_coord(1)) continue;
        q = CGAL::to_interval(pit->z());
        if(q.second < box.min_coord(2) || q.first > box.max_coord(2)) continue;

        Object_handle o = result.locate(*pit);
        Volume_const_iterator c;
        if(assign(c,o) && c->mark() == false) {
          convex_hull_3( points.begin(), points.end(), Ptemp);
          nubq2.add_polyhedron(Nef_polyhedron(Ptemp));
          added=true;
          break;
        }
      }
      if(!added) {
        convex_hull_3( points.begin(), points.end(), Ptemp);
        std::vector<Trunk_box> a, b;
        Halffacet_const_iterator fN;
        Facet_const_iterator fP;
        CGAL_forall_facets( fN, result)  a.push_back( Trunk_box( fN ) );
        for(fP = Ptemp.facets_begin(); fP!=Ptemp.facets_end(); ++fP)
          b.push_back( Trunk_box( fP ) );
        Callback<Trunk_box> callback;
        box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),callback);
        if(callback.intersect())
          nubq2.add_polyhedron(Nef_polyhedron(Ptemp));
      }
    }

    return nubq2.get_union();
  }
};

} //namespace CGAL
