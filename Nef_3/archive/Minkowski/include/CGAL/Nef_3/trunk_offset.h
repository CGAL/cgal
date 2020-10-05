#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/Nary_union.h>
#include <CGAL/Nef_3/Relabel_volume.h>
#include <CGAL/Nef_3/Mark_bounded_volumes.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>
#include <CGAL/Nef_3/Bounding_box_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>

namespace CGAL {

template <typename Kernel, typename K2,
          typename Vertex,
          typename Coordinate>
class Smaller_than
{
public:
  Smaller_than(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  bool operator()( const Vertex& v1, const Vertex& v2) {
    switch(coord) {
    case 0: return CGAL::compare_x(v1, v2) == SMALLER;
    case 1: return CGAL::compare_y(v1, v2) == SMALLER;
    case 2: return CGAL::compare_z(v1, v2) == SMALLER;
    default: CGAL_error();
    }
    return false;
  }

private:
  Coordinate coord;
};

template <typename K2, typename Vertex, typename Coordinate>
  class Smaller_than<CGAL::Lazy_kernel<typename K2::EK>, K2, Vertex, Coordinate>
{
public:
  Smaller_than(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  bool operator()( const Vertex& v1, const Vertex& v2) {
    switch(coord) {
    case 0: return CGAL::to_interval(v1.x()).second <
                           CGAL::to_interval(v2.x()).first;
    case 1: return CGAL::to_interval(v1.y()).second <
                           CGAL::to_interval(v2.y()).first;
    case 2: return CGAL::to_interval(v1.z()).second <
                           CGAL::to_interval(v2.z()).first;
    default: CGAL_error();
    }
    return false;
  }

private:
  Coordinate coord;
};

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
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;
  typedef Polyhedron_3<Kernel> Polyhedron;
  typedef typename Polyhedron::Vertex_const_iterator PVertex_const_iterator;
  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle Edge_const_handle;
  typedef typename Polyhedron::Halfedge_const_iterator Edge_const_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_const_circulator
    Halfedge_around_facet_const_circulator;
  typedef typename Nef_polyhedron::Object_handle Object_handle;
  typedef typename Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
  typedef typename Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Nef_polyhedron::Halffacet_const_iterator Halffacet_const_iterator;
  typedef typename Nef_polyhedron::Halffacet_const_handle Halffacet_const_handle;
  typedef typename Nef_polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Nef_polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
    SHalfedge_around_facet_const_circulator;
  typedef Bounding_box_3<Tag_true, Kernel> EBox;
  typedef Nary_union<Nef_polyhedron> NUBQ;
  typedef Relabel_volume<Nef_polyhedron> Relabel_volume;
  typedef Mark_bounded_volumes<Nef_polyhedron> Mark_bounded_volumes;

  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef SNC_intersection<SNC_structure> SNC_intersection;

  typedef Smaller_than<
    Kernel,
    Kernel,
    Point,
    int> Smaller_;

  class Trunk_box : public Box_intersection_d::Box_d< double, 3 > {

    Nef_polyhedron N;

    void extend( const Point& p) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( p.x() );
      q[1] = CGAL::to_interval( p.y() );
      q[2] = CGAL::to_interval( p.z() );
      Box_intersection_d::Box_d< double, 3 >::extend(q);
    }

  public:
    Halffacet_const_handle nf;
    Segment_3 seg;

    Trunk_box(Halffacet_const_handle f, Nef_polyhedron& N_) {
      N = N_;
      nf = f;
      SHalfedge_around_facet_const_circulator e(f->facet_cycles_begin()), end(e);
      CGAL_For_all(e,end)
        extend(e->source()->source()->point());
    }

    Trunk_box(Edge_const_handle e, Nef_polyhedron& N_) {
      N = N_;
      seg = Segment_3(e->vertex()->point(),
                      e->opposite()->vertex()->point());
      extend(e->vertex()->point());
      extend(e->opposite()->vertex()->point());
    }

    Trunk_box(Halffacet_const_handle f) {
      nf = f;
      SHalfedge_around_facet_const_circulator e(f->facet_cycles_begin()), end(e);
      CGAL_For_all(e,end)
        extend(e->source()->source()->point());
    }

    Trunk_box(Halfedge_const_handle e) {
      seg = Segment_3(e->source()->point(),
                      e->twin()->source()->point());
      extend(e->source()->point());
      extend(e->twin()->source()->point());
    }

    Nef_polyhedron get_polyhedron() { return N;}
  };


  template<typename Trunk_box, typename Nary_union>
  struct CallbackNfSeg {
    Nary_union& nu;
    Unique_hash_map<Halffacet_const_handle, bool> inserted;
    int count;

    CallbackNfSeg(Nary_union& nu_) : nu(nu_), inserted(false), count(0) {}

    void operator()( Trunk_box& box0, Trunk_box& box1 ) {
      if(inserted[box1.get_polyhedron().halffacets_begin()]) {
        return;
      }

      SNC_intersection is;
      Point ip;
      if(!is.does_intersect_internally(box1.seg, box0.nf, ip)) {
        return;
      }

      inserted[box1.get_polyhedron().halffacets_begin()] = true;
      Nef_polyhedron N(box1.get_polyhedron());
      std::cerr << "count " << ++count << std::endl;
      nu.add_polyhedron(N);
    }
  };

  template<typename Trunk_box, typename Nary_union>
  struct CallbackSegPf {
    Nary_union& nu;
    Unique_hash_map<Halffacet_const_handle, bool> inserted;
    int count;

    CallbackSegPf(Nary_union& nu_) : nu(nu_), inserted(false), count(0) {}

    void operator()( Trunk_box& box0, Trunk_box& box1 ) {
      if(inserted[box1.get_polyhedron().halffacets_begin()]) {
        return;
      }

      SNC_intersection is;
      Point ip;
      if(!is.does_intersect_internally(box0.seg, box1.nf, ip))
        return;

      inserted[box1.get_polyhedron().halffacets_begin()] = true;
      Nef_polyhedron N(box1.get_polyhedron());
      std::cerr << "count " << ++count << std::endl;
      nu.add_polyhedron(N);
    }
  };

  /*
  template<typename Trunk_box, typename Nary_union>
  struct CallbackSegSeg {
    Nary_union& nu;
    Unique_hash_map<Halffacet_const_handle, bool> inserted;
    int count;

    Callback(Nary_union& nu_) : nu(nu_), inserted(false), count(0) {}

    void operator()( Trunk_box& box0, Trunk_box& box1 ) {
      if(inserted[box1.get_polyhedron().halffacets_begin()]) {
        std::cerr << "already inserted" << std::endl;
        return;
      }

      SNC_intersection is;
      Point ip;
      if(!is.does_intersect_internally(box1.seg, box0.nf, ip))
        return;

      inserted[box1.get_polyhedron().halffacets_begin()] = true;
      Nef_polyhedron N(box1.get_polyhedron());
      std::cerr << "count " << ++count << std::endl;
      nu.add_polyhedron(N);
    }
  };
  */

  bool has_two_unmarked_volumes(Nef_polyhedron& N) {
    std::cerr << "number of volumes " << N.number_of_volumes() << std::endl;
    /*
    int argc = 0;
    char* argv[1];
    QApplication a(argc, argv);
    CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
      new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
    a.setMainWidget(w);
    w->show();
    a.exec();
    */
    Volume_const_iterator ci=++N.volumes_begin();
    while(ci != N.volumes_end() && ci->mark()) ++ci;
    return ci != N.volumes_end();
  }

  int mod, step, max_points;

  template <typename p_it, typename f_it, typename Hash>
  Nef_polyhedron first_approximation_with_mod(p_it pbegin, p_it pend,
                                              f_it fbegin, f_it fend,
                                              const Polyhedron& P, Hash& Done) {
    Polyhedron Ptemp;
    Add_points add;
    NUBQ nary_union;

    {
      std::cerr << "modulo " << mod << std::endl;
      int i=0;
      for(f_it f = fbegin; f != fend; ++f, ++i) {
        if(i%mod!=0) continue;
        std::cerr << i << std::endl;
        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);
        convex_hull_3( points.begin(), points.end(), Ptemp);
        nary_union.add_polyhedron(Nef_polyhedron(Ptemp));
        Done[i] = true;
      }
    }

    Nef_polyhedron result = nary_union.get_union();

    while(!has_two_unmarked_volumes(result) && mod != 1) {
      mod /= 2;
      std::cerr << "modulo " << mod << std::endl;
      int j = 0;
      for(f_it f = fbegin; f != fend; ++f, ++j) {
        if(j%mod!=0) continue;
        if(j%(2*mod)==0) continue;
        std::cerr << j << std::endl;
        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);
        convex_hull_3( points.begin(), points.end(), Ptemp);
        nary_union.add_polyhedron(Nef_polyhedron(Ptemp));
        Done[j] = true;
      }
      result = nary_union.get_union();
    }

    std::cerr << "found first approximation " << std::endl;

    if(mod > 1) {
      Relabel_volume rv(true);
      result.delegate(rv);

      CGAL_assertion(result.number_of_volumes() == 2);
      CGAL_assertion(result.volumes_begin()->mark());
      CGAL_assertion(!(++result.volumes_begin())->mark());
    }

    return result;
  }

  template <typename p_it>
  Nef_polyhedron first_approximation_from_bbox(p_it pbegin, p_it pend,
                                               const Polyhedron& P) {

    //    typename Kernel::FT q[3];

    p_it curr = pbegin;


    EBox bbt(*curr, *curr);
    for(++curr; curr != pend; ++curr) {
      bbt = bbt + EBox(*curr, *curr);
    }
    /*
    std::cerr << "bbt " << bbt.min_coord(0)
              << ", " << bbt.min_coord(1)
              << ", " << bbt.min_coord(2)
              << " - " << bbt.max_coord(0)
              << ", " << bbt.max_coord(1)
              << ", " << bbt.max_coord(2)
              << std::endl;
    */

    std::cerr << "bbt " << bbt.get_min()
              << " - " << bbt.get_max() << std::endl;

    EBox bbp;
    for(PVertex_const_iterator pvi = P.vertices_begin();
        pvi != P.vertices_end(); ++pvi) {
      bbp = bbp + EBox(pvi->point(), pvi->point());
    }

    /*
    std::cerr << "bbp " << bbp.min_coord(0)
              << ", " << bbp.min_coord(1)
              << ", " << bbp.min_coord(2)
              << " - " << bbp.max_coord(0)
              << ", " << bbp.max_coord(1)
              << ", " << bbp.max_coord(2)
              << std::endl;
    */

    std::cerr << "bbp " << bbp.get_min()
              << " - " << bbp.get_max() << std::endl;

    Nef_polyhedron cube;
    std::ifstream in("unit_cube_at_origin.nef3");
    in >> cube;

    /*
    Vector vec_scale = Vector((bbt.max_coord(0) - bbt.min_coord(0)) -
      (bbp.max_coord(0) - bbp.min_coord(0)),
      (bbt.max_coord(1) - bbt.min_coord(1)) -
      (bbp.max_coord(1) - bbp.min_coord(1)),
      (bbt.max_coord(2) - bbt.min_coord(2)) -
      (bbp.max_coord(2) - bbp.min_coord(2)));
    */

    Vector vec_scale((bbt.get_max() - bbt.get_min()) -
                     (bbp.get_max() - bbp.get_min()));

    std::cerr << "scale " << vec_scale << std::endl;

    Aff_transformation_3 scale(vec_scale.hx(), 0, 0,
                               0, vec_scale.hy(), 0,
                               0, 0, vec_scale.hz(),
                               vec_scale.hw());
    cube.transform(scale);
    /*
    Vector vec = Vector(bbt.min_coord(0),
                        bbt.min_coord(1),
                        bbt.min_coord(2));
    */
    Vector vec((bbt.get_min() - CGAL::ORIGIN) +
               (bbp.get_max() - CGAL::ORIGIN));
    std::cerr << "translate " << vec << std::endl;

    Aff_transformation_3 trans(CGAL::TRANSLATION, vec);
    cube.transform(trans);
    cube = (!cube).closure();
    return cube;
  }

  template <typename p_it>
  Nef_polyhedron compute_convex_hull
    (p_it pbegin, p_it pend) {

    Polyhedron CV;
    convex_hull_3( pbegin, pend, CV);
    return Nef_polyhedron(CV);
  }

  template <typename p_it>
  Nef_polyhedron first_approximation_from_convex_hull
    (p_it pbegin, p_it pend, const Polyhedron& P) {

    Nef_polyhedron NCV = compute_convex_hull(pbegin, pend);
    NCV = !NCV;

    EBox bbp;
    for(PVertex_const_iterator pvi = P.vertices_begin();
        pvi != P.vertices_end(); ++pvi) {
      bbp = bbp + EBox(pvi->point(), pvi->point());
    }

    std::cerr << "bbp " << bbp.get_min()
              << " - " << bbp.get_max() << std::endl;

    Nef_polyhedron result(NCV);

    Vector vec(bbp.get_max() - CGAL::ORIGIN);
    std::cerr << "translate " << vec << std::endl;
    Aff_transformation_3 trans(CGAL::TRANSLATION, vec);
    result.transform(trans);
    result = result.join(NCV);
    return result;
  }

  template <typename p_it, typename f_it>
    Nef_polyhedron first_approximation_from_kdtree
    (p_it pbegin, p_it pend,
     f_it fbegin, f_it fend, const Polyhedron& P) {
    NUBQ nary_union;
    std::vector<Point> point_copies;
    for(p_it curr = pbegin; curr != pend; ++curr)
      point_copies.push_back(*curr);

    recursive_approximation_from_kdtree
      (pbegin, point_copies.begin(), point_copies.end(),
       fbegin, fend, nary_union, 0);

    Nef_polyhedron NCV(nary_union.get_union());
    Mark_bounded_volumes mbv(true);
    NCV.delegate(mbv);

    CGAL_assertion(NCV.number_of_volumes() == 2);
    CGAL_assertion(!NCV.volumes_begin()->mark());
    CGAL_assertion((++NCV.volumes_begin())->mark());

    PVertex_const_iterator pvi = P.vertices_begin();
    FT q[3];
    q[0] = pvi->point().x();
    q[1] = pvi->point().y();
    q[2] = pvi->point().z();
    EBox bbp(q);
    for(++pvi; pvi != P.vertices_end(); ++pvi) {
      bbp.extend(pvi->point());
    }

    //    std::cerr << "bbp " << bbp.min_coord()
    //                  << " - " << bbp.get_max() << std::endl;

    Nef_polyhedron result(NCV);
    Vector vec(bbp.max_coord(0),
               bbp.max_coord(1),
               bbp.max_coord(2));
    std::cerr << "translate " << vec << std::endl;
    Aff_transformation_3 trans(CGAL::TRANSLATION, vec);
    result.transform(trans);
    result = result.intersection(NCV);

    //    std::cout << NCV << result;

    return result;
  }

  template <typename c_it, typename p_it, typename f_it, class NaryUnion>
  void recursive_approximation_from_kdtree
    (c_it cbegin, p_it pbegin, p_it pend, f_it fbegin, f_it fend,
     NaryUnion& nunion, int depth) {
    typedef typename std::iterator_traits<f_it>::value_type value_type;
    typedef typename value_type::first_type intp;
    std::cerr << "recursive " << std::distance(pbegin, pend) << std::endl;
    if(std::distance(pbegin, pend) < max_points) {
      std::list<Point> tmp_points;
      for(f_it curr = fbegin; curr != fend; ++curr)
        for(intp ip=curr->first; ip != curr->second; ++ip)
          tmp_points.push_back(*(cbegin+*ip));
      nunion.add_polyhedron(compute_convex_hull
                            (tmp_points.begin(),
                             tmp_points.end()));
    } else {
      Smaller_ smaller_(depth%3);
      std::nth_element(pbegin,
                       pbegin+std::distance(pbegin, pend)/2,
                       pend,
                       smaller_);
      p_it median = pbegin+ std::distance(pbegin, pend)/2;

      std::list<value_type> f1, f2;
      for(f_it curr = fbegin; curr != fend; ++curr) {
        bool b1(false), b2(false);
        for(intp ip=curr->first; ip != curr->second; ++ip)
          if(smaller_(*median, *(cbegin+*ip)))
            b2|=true;
          else
            b1|=true;
        if(b1)
          f1.push_back(*curr);
        if(b2)
          f2.push_back(*curr);
      }

      recursive_approximation_from_kdtree
        (cbegin, pbegin, median, f1.begin(), f1.end(), nunion, depth+1);
      recursive_approximation_from_kdtree
        (cbegin, median, pend, f2.begin(), f2.end(), nunion, depth+1);

    }
  }

  /*
  template <typename Nary_union>
  class Facet_check_by_range_query {

    const Nef_polyhedron& N;
    Nary_union& nu;

  public:
    Facet_check_by_range_query(const Nef_polyhedron& N_,
                               Nary_union& nu_)
      : N(N_), nu(nu_) {}

      void add_box(const std::vector<Point>& points,
                   const Box& box)
      {
        if(N.in_range(box)) {
          Polyhedron Ptemp;
          convex_hull_3( points.begin(), points.end(), Ptemp);
          nu.add_polyhedron(Nef_polyhedron(Ptemp));
        }
      }

      void add_polyhedron() {}
  }
  */

  template <typename Nary_union>
  class Facet_check_by_pl_and_rs {

    std::vector<Trunk_box> a, p_edges, p_facets;
    const Nef_polyhedron& N;
    Nary_union& nu;

  public:
    Facet_check_by_pl_and_rs(const Nef_polyhedron& N_,
                             Nary_union& nu_)
      : N(N_), nu(nu_) {}

    void add_box(const std::vector<Point>& points)
    {
      typedef typename std::vector<Point>::const_iterator p_it2;
      for(p_it2 pit = points.begin(); pit != points.end(); ++pit) {
        /*
        std::pair<double, double> q = CGAL::to_interval(pit->x());
        if(q.second < box.min_coord(0) ||
           q.first > box.max_coord(0)) continue;
        q = CGAL::to_interval(pit->y());
        if(q.second < box.min_coord(1) ||
           q.first > box.max_coord(1)) continue;
        q = CGAL::to_interval(pit->z());
        if(q.second < box.min_coord(2) ||
           q.first > box.max_coord(2)) continue;
        */

        Object_handle o = N.locate(*pit);
        Volume_const_iterator c;
        if(assign(c,o) && c->mark() == false) {
          Polyhedron Ptemp;
          convex_hull_3( points.begin(), points.end(), Ptemp);
          nu.add_polyhedron(Nef_polyhedron(Ptemp));
          return;
        }
      }

      Polyhedron Ptmp;
      convex_hull_3( points.begin(), points.end(), Ptmp);
      Nef_polyhedron Ntmp(Ptmp);
      Edge_const_iterator eP;
      for(eP = Ptmp.edges_begin(); eP!=Ptmp.edges_end(); ++eP)
        p_edges.push_back( Trunk_box( eP, Ntmp) );
      Halffacet_const_iterator fP;
      for(fP = Ntmp.halffacets_begin(); fP!=Ntmp.halffacets_end(); ++fP)
        if(fP->is_twin())
           p_facets.push_back( Trunk_box( fP, Ntmp));
    }

    void add_polyhedra() {
      std::cerr << " box_intersection" << std::endl;

      a.clear();
      Halffacet_const_iterator fN;
      CGAL_forall_facets( fN, N)
        a.push_back( Trunk_box( fN ) );

      CallbackNfSeg<Trunk_box, Nary_union> callback0(nu);
      box_intersection_d( a.begin(), a.end(),
                          p_edges.begin(), p_edges.end(),
                          callback0);

      a.clear();
      Halfedge_const_iterator eN;
      CGAL_forall_edges( eN, N)
        a.push_back( Trunk_box( eN) );

      CallbackSegPf<Trunk_box, Nary_union> callback1(nu);
      box_intersection_d( a.begin(), a.end(),
                          p_facets.begin(), p_facets.end(),
                          callback1);
      a.clear();

    }
  };


  template <typename Nary_union>
  class Facet_check_by_pl_and_rs2 {

    std::vector<Trunk_box> a, b, p_edges, p_facets;
    const Nef_polyhedron& N;
    Nary_union& nu;

  public:
    Facet_check_by_pl_and_rs2(const Nef_polyhedron& N_,
                             Nary_union& nu_)
      : N(N_), nu(nu_)
      {
        Halffacet_const_iterator fN;
        CGAL_forall_facets( fN, N)
          a.push_back( Trunk_box( fN ) );

        Halfedge_const_iterator eN;
        CGAL_forall_edges( eN, N)
          b.push_back( Trunk_box( eN) );
      }

    void add_box(const std::vector<Point>& points)
    {
      typedef typename std::vector<Point>::const_iterator p_it2;
      for(p_it2 pit = points.begin(); pit != points.end(); ++pit) {
        /*
        std::pair<double, double> q = CGAL::to_interval(pit->x());
        if(q.second < box.min_coord(0) ||
           q.first > box.max_coord(0)) continue;
        q = CGAL::to_interval(pit->y());
        if(q.second < box.min_coord(1) ||
           q.first > box.max_coord(1)) continue;
        q = CGAL::to_interval(pit->z());
        if(q.second < box.min_coord(2) ||
           q.first > box.max_coord(2)) continue;
        */

        Object_handle o = N.locate(*pit);
        Volume_const_iterator c;
        if(assign(c,o) && c->mark() == false) {
          Polyhedron Ptemp;
          convex_hull_3( points.begin(), points.end(), Ptemp);
          nu.add_polyhedron(Nef_polyhedron(Ptemp));
          return;
        }
      }

      /*
      Polyhedron Ptmp;
      convex_hull_3( points.begin(), points.end(), Ptmp);
      Nef_polyhedron Ntmp(Ptmp);

      Edge_const_iterator eP;
      for(eP = Ptmp.edges_begin(); eP!=Ptmp.edges_end(); ++eP)
        p_edges.push_back( Trunk_box( eP, Ntmp) );
      CallbackNfSeg<Trunk_box, Nary_union> callback0(nu);
      box_intersection_d( a.begin(), a.end(),
                          p_edges.begin(), p_edges.end(),
                          callback0);
      p_edges.clear();

      Halffacet_const_iterator fP;
      for(fP = Ntmp.halffacets_begin(); fP!=Ntmp.halffacets_end(); ++fP)
        if(fP->is_twin())
           p_facets.push_back( Trunk_box( fP, Ntmp));
      CallbackSegPf<Trunk_box, Nary_union> callback1(nu);
      box_intersection_d( b.begin(), b.end(),
                          p_facets.begin(), p_facets.end(),
                          callback1);
      p_facets.clear();
      */
    }

    void add_polyhedra() {
      std::cerr << " box_intersection" << std::endl;

    }
  };

  template<typename Point_iterator>
  bool intersects(Point_iterator begin, Point_iterator end, const EBox& box)
  {
    Point_iterator pi;
    bool result = false;
    for(pi = begin; pi != end; ++pi) {
      //      std::cerr << "points " << *pi << std::endl;
      if(pi->x() > box.min_coord(0)) {
        result = true;
        break;
      }
    }
    if(!result) return false;

    result = false;
    for(pi = begin; pi != end; ++pi)
      if(pi->y() > box.min_coord(1)) {
        result = true;
        break;
      }
    if(!result) return false;

    result = false;
    for(pi = begin; pi != end; ++pi)
      if(pi->z() > box.min_coord(2)) {
        result = true;
        break;
      }
    if(!result) return false;

    result = false;
    for(pi = begin; pi != end; ++pi)
      if(pi->x() < box.max_coord(0)) {
        result = true;
        break;
      }
    if(!result) return false;

    result = false;
    for(pi = begin; pi != end; ++pi)
      if(pi->y() < box.max_coord(1)) {
        result = true;
        break;
      }
    if(!result) return false;

    result = false;
    for(pi = begin; pi != end; ++pi)
      if(pi->z() < box.max_coord(2)) {
        result = true;
        break;
      }
    return result;
  }

  template <typename Point_iterator>
  bool intersects2(const Polyhedron& P, Point_iterator pbegin, Point_iterator pend)
  {
    Point_iterator pi;
    Facet_const_iterator fi;
    for(fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
      bool result = false;
      for(pi = pbegin; pi != pend; ++pi)
        if(fi->plane().oriented_side(*pi) == CGAL::ON_NEGATIVE_SIDE) {
          result = true;
          break;
        }
      if(!result) return false;
    }

    return true;
  }

 public:
  Trunk_offset(int modulo=1, int step_=0, int mp = 500)
    : mod(modulo), step(step_), max_points(mp) {}

  template<typename p_it, typename f_it>
    Nef_polyhedron operator()(p_it pbegin, p_it pend,
                              f_it fbegin, f_it fend,
                              const Polyhedron& P) {

    std::map<int, bool> Done;
    int j = 0;
    for(f_it f = fbegin; f != fend; ++f, ++j) Done[j] = false;


    //    Nef_polyherdon result =
    //      first_approximation_from_convex_hull(pbegin, pend, P);
    //    Nef_polyhedron result =
    //      first_approximation_from_kdtree(pbegin, pend,
    //                                      fbegin, fend, P);


    Nef_polyhedron result =
      first_approximation_with_mod(pbegin, pend, fbegin, fend, P, Done);

    if(mod == 1) {
      CGAL_assertion(result.number_of_volumes() == 2);
      CGAL_assertion(!result.volumes_begin()->mark());
      CGAL_assertion((++result.volumes_begin())->mark());

      return result;
    }

    CGAL_assertion(result.number_of_volumes() == 2);
    CGAL_assertion(result.volumes_begin()->mark());
    CGAL_assertion(!(++result.volumes_begin())->mark());

    typedef NUBQ Nary_union;
    Nary_union nary_union2;
    nary_union2.add_polyhedron(result);

    int tmp_mod = mod;

    Add_points add;
    bool first = false;

    /*
    do {
      int k=0;
      int count = 0;
      if(!first) mod/=step;
      std::cerr << "modulo " << mod << std::endl;

      for(f_it f = fbegin; f != fend; ++f, ++k) {
        if(k%(mod)!=0 || (!first && k%(mod*step)==0)) continue;
        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);

        typedef typename std::vector<Point>::const_iterator p_it2;
        for(p_it2 pit = points.begin(); pit != points.end(); ++pit) {

          Object_handle o = result.locate(*pit);
          Volume_const_iterator c;
          if(assign(c,o) && c->mark() == false) {
            Polyhedron Ptemp;
            convex_hull_3( points.begin(), points.end(), Ptemp);
            nary_union2.add_polyhedron(Nef_polyhedron(Ptemp));
            Done[k] = true;
            std::cerr << "added facet " << k << ", " << ++count << std::endl;
            break;
          }
        }
        //        facet_check.add_box(points);
      }
      //      facet_check.add_polyhedra();

      first = false;
      result = nary_union2.get_union();

      CGAL_assertion(result.volumes_begin()->mark());

    } while(mod != 1);


    mod = tmp_mod;
    */

    int count = 0;
    do {

      EBox box;
      for(Vertex_const_iterator nv = result.vertices_begin(); nv != result.vertices_end(); ++nv) {
        box.extend(nv->point());
      }

      std::list<Point> boxPoints;
      boxPoints.push_back(Point(box.min_coord(0), box.min_coord(1), box.min_coord(2)));
      boxPoints.push_back(Point(box.min_coord(0), box.min_coord(1), box.max_coord(2)));
      boxPoints.push_back(Point(box.min_coord(0), box.max_coord(1), box.min_coord(2)));
      boxPoints.push_back(Point(box.min_coord(0), box.max_coord(1), box.max_coord(2)));
      boxPoints.push_back(Point(box.max_coord(0), box.min_coord(1), box.min_coord(2)));
      boxPoints.push_back(Point(box.max_coord(0), box.min_coord(1), box.max_coord(2)));
      boxPoints.push_back(Point(box.max_coord(0), box.max_coord(1), box.min_coord(2)));
      boxPoints.push_back(Point(box.max_coord(0), box.max_coord(1), box.max_coord(2)));

      std::cerr << "box " << box.min_coord(0)
                << ", " << box.min_coord(1)
                << ", " << box.min_coord(2)
                << " - " << box.max_coord(0)
                << ", " << box.max_coord(1)
                << ", " << box.max_coord(2) << std::endl;

      int ii = 0;
      mod/=step;
      std::cerr << "modulo " << mod << std::endl;

      for(f_it f = fbegin; f != fend; ++f, ++ii) {
        if((ii%mod) != 0) continue;
        if(Done[ii]) continue;

        std::vector<Point> points;
        fold_indices_polyhedron( pbegin,
                                 f->first, f->second,
                                 P.vertices_begin(), P.vertices_end(),
                                 back_inserter( points),
                                 add);

        Polyhedron Ptemp;
        convex_hull_3( points.begin(), points.end(), Ptemp);
        if(intersects(points.begin(), points.end(), box) ||
           intersects2(P, boxPoints.begin(), boxPoints.end())) {
          nary_union2.add_polyhedron(Nef_polyhedron(Ptemp));
          std::cerr << "added facet " << ii << ", " << ++count << std::endl;
        } else {
          /*
            Nef_polyhedron N(Ptemp);
            Nef_polyhedron empty = N - result;
            if(!empty.is_empty()) {
            std::cerr << empty();
            CGAL_assertion(false);
            }
          */
        }
      }

      result = nary_union2.get_union();
      CGAL_assertion(result.volumes_begin()->mark());

    } while(mod != 1);

    result = !result;

    CGAL_assertion(result.number_of_volumes() == 2);
    CGAL_assertion(!result.volumes_begin()->mark());
    CGAL_assertion((++result.volumes_begin())->mark());

    return result;
  }
};

} //namespace CGAL
