#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H 1

#include <CGAL/basic.h>
#include <set>
#include <CGAL/Segment_Delaunay_graph_storage_site_2.h>
#include <CGAL/Segment_Delaunay_graph_simple_storage_site_2.h>


CGAL_BEGIN_NAMESPACE

template<class STraits>
class Construct_storage_site_2
{
public:
  typedef STraits                                    Storage_traits;
  typedef typename Storage_traits::Storage_site_2    Storage_site_2;
  typedef typename Storage_traits::Point_handle      Point_handle;
  typedef typename Storage_traits::Geom_traits       Geom_traits;

  typedef Storage_site_2                             result_type;
  //  struct Arity {};

protected:
  typedef typename Geom_traits::Intersections_tag    ITag;

  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4);
  }

  inline
  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4,
			const Point_handle& h5,
			const Point_handle& h6, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4, h5, h6);
  }

  inline
  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4,
			bool is_first_exact, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4,
						    is_first_exact);
  }

  
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&, const Tag_false&) const {
    CGAL_assertion( false );
    return Storage_site_2();
  }

  inline
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&, const Tag_false&) const {
    CGAL_assertion( false );
    return Storage_site_2();
  }

  inline
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			bool is_first_exact, const Tag_false&) const {
    CGAL_assertion( false );
    return Storage_site_2();
  }

public:
  inline
  result_type operator()(const Point_handle& h) const {
    return Storage_site_2::construct_storage_site_2(h);
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2) const {
    return Storage_site_2::construct_storage_site_2(h1, h2);
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4) const {
    return construct(h1, h2, h3, h4, ITag());
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4,
			 const Point_handle& h5,
			 const Point_handle& h6) const {
    return construct(h1, h2, h3, h4, h5, h6, ITag());
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4,
			 bool is_first_exact) const {
    return construct(h1, h2, h3, h4, is_first_exact, ITag());
  }

  // constructs the point of intersection
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1) const {
    CGAL_precondition( ss0.is_segment() && ss1.is_segment() );
    return Storage_site_2::construct_storage_site_2
      ( ss0.source_of_supporting_site(),
	ss0.target_of_supporting_site(),
	ss1.source_of_supporting_site(),
	ss1.target_of_supporting_site() );
  }

  Storage_site_2
  split_on_point_first_subsegment(const Storage_site_2& s,
				  const Storage_site_2& p) const
  {
    // Splits the first storage site which is a segment using the
    // second storage site which is an exact point
    // Two new storage sites are created, corresponding to the two
    // subsegments
    CGAL_precondition( s.is_segment() && p.is_point() );

    // computing the first sub-segment
    if ( s.is_input(0) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	return operator()(s.source_of_supporting_site(), p.point());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   true );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	return operator()( s.source_of_supporting_site(),
			   p.point(),
			   s.source_of_crossing_site(0),
			   s.target_of_crossing_site(0),
			   false );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   s.source_of_crossing_site(0),
			   s.target_of_crossing_site(0),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site() );
      }
    }
  }

  Storage_site_2
  split_on_point_second_subsegment(const Storage_site_2& s,
				   const Storage_site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    // computing the second sub-segment
    if ( s.is_input(1) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	return operator()(p.point(), s.target_of_supporting_site());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   false );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	return operator()( p.point(),
			   s.target_of_supporting_site(),
			   s.source_of_crossing_site(1),
			   s.target_of_crossing_site(1),
			   true );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   s.source_of_crossing_site(1),
			   s.target_of_crossing_site(1) );
      }
    }
  }

  Storage_site_2
  construct_first_subsegment(const Storage_site_2& ss0,
			     const Storage_site_2& ss1,
			     const Tag_true&) const
  {
    if ( ss0.is_input(0) ) {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(), true );
    } else {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss0.source_of_crossing_site(0),
	  ss0.target_of_crossing_site(0),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site() );
    }
  }

  Storage_site_2
  construct_second_subsegment(const Storage_site_2& ss0,
			      const Storage_site_2& ss1,
			      const Tag_true&) const
  {
    if ( ss0.is_input(1) ) {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(), false );
    } else {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(),
	  ss0.source_of_crossing_site(1),
	  ss0.target_of_crossing_site(1) );
    }
  }

  Storage_site_2
  construct_first_subsegment(const Storage_site_2&,
			     const Storage_site_2&,
			     const Tag_false&) const
  {
    CGAL_assertion( false );
    return Storage_site_2();
  }

  Storage_site_2
  construct_second_subsegment(const Storage_site_2&,
			      const Storage_site_2&,
			      const Tag_false&) const
  {
    CGAL_assertion( false );
    return Storage_site_2();
  }

  // constructs the subsegment with supporting segment ss0 and
  // endpoints the point of intersection of ss1 and ss0; the boolean
  // determines if the first or segment subsegment is constructed
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1,
			 bool first) const {
    //    CGAL_precondition( ss0.is_segment() && ss1.is_segment() );
    CGAL_precondition( ss0.is_segment() );
    if ( ss1.is_point() ) {
      if ( first ) {
	return split_on_point_first_subsegment(ss0, ss1);
      } else {
	return split_on_point_second_subsegment(ss0, ss1);
      }
    }

    if ( first ) {
      return construct_first_subsegment(ss0, ss1, ITag());
    } else {
      return construct_second_subsegment(ss0, ss1, ITag());
    }
  }

};

//----------------------------------------------------------------------
//----------------------------------------------------------------------


template<class STraits>
class Storage_site_split_2
{
  typedef STraits                                      Storage_traits_2;
  typedef typename Storage_traits_2::Storage_site_2    Storage_site_2;
  typedef typename Storage_traits_2::Point_handle      Point_handle;
  typedef typename Storage_traits_2::Geom_traits       Geom_traits;

  typedef typename Geom_traits::Intersection_tag       Intersection_tag;

  typedef Storage_site_2                               result_type;
  struct Arity {};

  Storage_site_split_2(const Storage_traits_2& straits)
    : straits_(straits) {}

private:
  std::pair<Storage_site_2,Storage_site_2>
  split_on_point(const Storage_site_2& s,
		 const Storage_site_2& p,
		 const Tag_false&) const
  {
    // Splits the first storage site which is a segment using the
    // second storage site which is an exact point
    // Two new storage sites are created, corresponding to the two
    // subsegments
    CGAL_precondition( s.is_segment() && p.is_point() );
    CGAL_precondition( p.is_input() );

    typename Storage_traits_2::Construct_storage_site_2 ctr =
      straits_.construct_storage_site_2_objec();

    Storage_site_2 s1 = ctr(s.source_of_supporting_site(), p.point());
    Storage_site_2 s2 = ctr(p.point(), s.target_of_supporting_site());

    return std::make_pair(s1, s2);
  }

  std::pair<Storage_site_2,Storage_site_2>
  split_on_point(const Storage_site_2& s,
		 const Storage_site_2& p,
		 const Tag_true&) const
  {
    // Splits the first storage site which is a segment using the
    // second storage site which is an exact point
    // Two new storage sites are created, corresponding to the two
    // subsegments
    CGAL_precondition( s.is_segment() && p.is_point() );

    typename Storage_traits_2::Construct_storage_site_2 ctr =
      straits_.construct_storage_site_2_objec();

    // computing the first sub-segment
    Storage_site_2 s1;
    if ( s.is_input(0) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	s1 = ctr(s.source_of_supporting_site(), p.point());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  straits_.geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	s1 = ctr( supp0.source_of_supporting_site(),
		  supp0.target_of_supporting_site(),
		  supp1.source_of_supporting_site(),
		  supp1.target_of_supporting_site(),
		  true );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	s1 = ctr( s.source_of_supporting_site(),
		  p.point(),
		  s.source_of_crossing_site(0),
		  s.target_of_crossing_site(0),
		  false );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  straits_.geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	s1 = ctr( supp0.source_of_supporting_site(),
		  supp0.target_of_supporting_site(),
		  s.source_of_crossing_site(0),
		  s.target_of_crossing_site(0),
		  supp1.source_of_supporting_site(),
		  supp1.target_of_supporting_site() );
      }
    }

    Storage_site_2 s2;
    // computing the second sub-segment
    if ( s.is_input(1) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	s2 = ctr(p.point(), s.source_of_supporting_site());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  straits_.geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	s2 = ctr( supp0.source_of_supporting_site(),
		  supp0.target_of_supporting_site(),
		  supp1.source_of_supporting_site(),
		  supp1.target_of_supporting_site(),
		  false );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	s2 = ctr( p.point(),
		  s.target_of_supporting_site(),
		  s.source_of_crossing_site(1),
		  s.target_of_crossing_site(1),
		  true );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  straits_.geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	s2 = ctr( supp0.source_of_supporting_site(),
		  supp0.target_of_supporting_site(),
		  supp1.source_of_supporting_site(),
		  supp1.target_of_supporting_site(),
		  s.source_of_crossing_site(1),
		  s.target_of_crossing_site(1) );
      }
    }
    return std::make_pair(s1, s2);
  }

public:
  std::pair<Storage_site_2,Storage_site_2>
  operator()(const Storage_site_2& t1, const Storage_site_2& t2) const
  {
    CGAL_precondition( t1.is_segment() );
    CGAL_precondition( t2.is_point() );
    return split_on_point(t1, t2, Intersection_tag());
  }

private:
  const Storage_traits_2& straits_;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace CGALi {

  template<class Gt, class USE_SIMPLE_STORAGE_SITE_Tag>
  struct SDGST2_Which_storage_site;

  // use the simple storage site
  template<class Gt>
  struct SDGST2_Which_storage_site<Gt,Tag_false>
  {
    typedef Gt         Geom_traits;
    typedef Tag_false  Storage_site_tag;

    typedef
    Segment_Delaunay_graph_simple_storage_site_2<Geom_traits>
    Storage_site_2;
  };

  // use the full storage site
  template<class Gt>
  struct SDGST2_Which_storage_site<Gt,Tag_true>
  {
    typedef Gt         Geom_traits;
    typedef Tag_true   Storage_site_tag;

    typedef
    Segment_Delaunay_graph_storage_site_2<Geom_traits>
    Storage_site_2;
  };


} // namespace CGALi


//----------------------------------------------------------------------

template<class Gt>
class Storage_traits_2
{
public:
  typedef Gt                                       Geom_traits;
  typedef typename Geom_traits::Point_2            Point_2;
  typedef typename Geom_traits::Site_2             Site_2;
  typedef std::set<Point_2>                        Point_container;
  typedef typename Point_container::iterator       Point_handle;

private:
  typedef Storage_traits_2<Geom_traits>            Self;
  typedef typename Geom_traits::Intersections_tag  ITag;

public:
  typedef typename
  CGALi::SDGST2_Which_storage_site<Self,ITag>::Storage_site_2
  Storage_site_2;

  typedef Construct_storage_site_2<Self>      Construct_storage_site_2;

  // MK::FIGURE OUT HOW TO PASS A REFERENCE TO GEOM_TRAITS AND HAVE
  // DEFAULT CONSTRUCTOR AS WELL IF POSSIBLE
  Storage_traits_2(const Geom_traits& gt = Geom_traits()) : gt_(gt) {}

  inline const Geom_traits& geom_traits() const { return gt_; }

  inline Construct_storage_site_2
  construct_storage_site_2_object() const {
    return Construct_storage_site_2();
  }

private:
  Geom_traits gt_;
};



CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H
