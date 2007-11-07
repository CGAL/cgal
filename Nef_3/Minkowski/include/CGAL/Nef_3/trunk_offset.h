#include <CGAL/convex_hull_3.h> 
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/Nary_union_by_queue.h>
#include <CGAL/Nef_3/Nary_union_by_pq.h>
//#include <CGAL/Nef_3/Nary_union_by_pq_me.h>
#include <CGAL/Nef_3/Nary_union_using_distinct_uniter.h>
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

CGAL_BEGIN_NAMESPACE

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

template<typename Trunk_box>
struct Callback {
  bool intersection;
  Callback() : intersection(false) {}
  void operator()( Trunk_box& box0, Trunk_box& box1 ) {intersection=true;}
  bool intersect() const {return intersection;}
};


template<typename Nef_polyhedron>
class Trunk_offset {
  
  typedef typename Nef_polyhedron::Kernel Kernel;
  typedef typename Kernel::FT FT;
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
  //  typedef typename Box_intersection_d::Box_d<FT,3> EBox;
  typedef Bounding_box_3<typename Kernel::Kernel_tag, Kernel> EBox;
  typedef Nary_union_by_queue<Nef_polyhedron> NUBQ;
  typedef Nary_union_using_distinct_uniter<Nef_polyhedron> NUUDU;
  typedef Relabel_volume<Nef_polyhedron> Relabel_volume;
  typedef Mark_bounded_volumes<Nef_polyhedron> Mark_bounded_volumes;

  typedef Smaller_than<
    Kernel,
    Kernel,
    Point, 
    int> Smaller_;
  
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

  int mod, step, max_points;

  template <typename p_it, typename f_it>
  Nef_polyhedron first_approximation_with_mod(p_it pbegin, p_it pend,
					      f_it fbegin, f_it fend,
					      const Polyhedron& P) {
    Polyhedron Ptemp;
    Add_points add;
#ifdef CGAL_NEF3_NARY_UNION_USING_DU
    NUUDU nary_union;
#else
    NUBQ nary_union;
#endif

    {
      int i=0;
      for(f_it f = fbegin; f != fend; ++f, ++i) {
	if(i%mod!=0) continue;
	//	std::cerr << i << std::endl;
	std::vector<Point> points;
	fold_indices_polyhedron( pbegin,
				 f->first, f->second,
				 P.vertices_begin(), P.vertices_end(),
				 back_inserter( points),
				 add);
	convex_hull_3( points.begin(), points.end(), Ptemp); 
	nary_union.add_polyhedron(Nef_polyhedron(Ptemp));
      }
    }
    
    Nef_polyhedron result = nary_union.get_union();
    
    while(!has_two_unmarked_volumes(result) && mod != 1) {
      //      std::cerr << result;
      mod /= 2;
      int j = 0;
      for(f_it f = fbegin; f != fend; ++f, ++j) {
	if(j%mod!=0 || (j%(2*mod))==0 || (mod!=1)) continue;
	std::cerr << j << std::endl;
	std::vector<Point> points;
	fold_indices_polyhedron( pbegin,
				 f->first, f->second,
				 P.vertices_begin(), P.vertices_end(),
				 back_inserter( points),
				 add);
	convex_hull_3( points.begin(), points.end(), Ptemp); 
	nary_union.add_polyhedron(Nef_polyhedron(Ptemp));
      }
      result = nary_union.get_union();
    }

    Relabel_volume rv(true);
    result.delegate(rv);    
    
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
#ifdef CGAL_NEF3_NARY_UNION_USING_DU
    NUUDU nary_union;
#else
    NUBQ nary_union;
#endif
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
    result = result.intersection(NCV);

    //    std::cout << NCV << result;

    /*
    int argc = 0;
    char* argv[1];
    QApplication a(argc, argv);
    CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
      new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
    a.setMainWidget(w);
    w->show();
    a.exec();
    */

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
  
 public:
  Trunk_offset(int modulo=1, int step_=0, int mp = 500) 
    : mod(modulo), step(step_), max_points(mp) {}

  template<typename p_it, typename f_it>
    Nef_polyhedron operator()(p_it pbegin, p_it pend,
			      f_it fbegin, f_it fend,
			      const Polyhedron& P) {
    
    //    Nef_polyhedron result = 
    //      first_approximation_from_convex_hull(pbegin, pend, P);
    Nef_polyhedron result = 
      first_approximation_from_kdtree(pbegin, pend, 
				      fbegin, fend, P);

    CGAL_assertion(result.number_of_volumes() == 2);
    CGAL_assertion(!result.volumes_begin()->mark());
    CGAL_assertion((++result.volumes_begin())->mark());

    result = !result;
    result = result.interior();

    CGAL_assertion(result.number_of_volumes() == 2);
    CGAL_assertion(result.volumes_begin()->mark());
    CGAL_assertion(!(++result.volumes_begin())->mark());

    if(mod == 1)
      return result;

#ifdef CGAL_NEF3_NARY_UNION_USING_DU
    NUUDU nary_union2;
#else
    NUBQ nary_union2;
#endif
    nary_union2.add_polyhedron(result);

    Box box;
    for(Vertex_const_iterator nv = result.vertices_begin(); nv != result.vertices_end(); ++nv) {
      std::pair<double, double> q[3];
      q[0] = CGAL::to_interval( nv->point().x() );
      q[1] = CGAL::to_interval( nv->point().y() );
      q[2] = CGAL::to_interval( nv->point().z() );
      box.extend(q);
    }

    Add_points add;
    bool first = true;
    do {
      int k=1;
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
	    Polyhedron Ptemp;
	    convex_hull_3( points.begin(), points.end(), Ptemp);
	    Nef_polyhedron Ntemp(Ptemp);
	    nary_union2.add_polyhedron(Ntemp);
	    added=true;
	    break;
	  }
	}
	
	/*
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
	*/
      }
      first = false;      
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

CGAL_END_NAMESPACE
