// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_simple_polygons_2.h
// package       : bops (2.2)
// source        : include/CGAL/bops_simple_polygons_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_SIMPLE_POLYGONS_2_H
#define CGAL_BOPS_SIMPLE_POLYGONS_2_H

//#define CGAL__BOPS_DEBUG_ON
//#define CGAL__DCEL_DEBUG_ON
//#define _DCEL__V2E_DEBUG_ON
//#define CGAL__INTERSECTING_POLYGONS_DEBUG_ON

#include <list>
#include <vector>
#include <CGAL/bops_dcel.h>
#include <CGAL/nsquare_intersecting.h>
#include <CGAL/intersecting_polygons.h>

CGAL_BEGIN_NAMESPACE

template<class I>
class Bops_Polygons_2 : public I {
public:
	typedef typename I::R                         R;
	//typedef 		 I::NT FT;
	typedef typename I::Point                     _Point_2;
	typedef typename I::Segment                   Segment_2;
	typedef typename I::Object                    Object;
	typedef typename I::Output_polygon_container  Polygon_Container;
	typedef typename I::Output_polygon            _Polygon_2;
	typedef typename I::Output_object_container   Output_container;

	typedef typename Output_container::const_iterator    result_iterator;
	typedef typename Output_container::size_type         size_type;

	result_iterator begin()           const { return _result.begin(); }
	result_iterator end()             const { return _result.end(); }
	const Output_container& result()  const { return _result; }
	size_type size()                  const { return _result.size(); }
	bool empty()                      const { return _result.empty(); }
	
	void add_to_result( const _Polygon_2& pgon) {
		_result.push_back(I::Make_object(pgon) );
	}
	
	void add_to_result( const std::list<_Point_2>& l) {
		_result.push_back(I::Make_object(_Polygon_2(l.begin(),l.end())) );
	}
	
	virtual bool operation() = 0;
	virtual ~Bops_Polygons_2() {}
protected:
	/*
	 *	Intersection-Type of Polygons A,B
	 */
	enum Intersection_type {
		is_empty         = 0,
		is_identical     = 1,
		A_is_subset_of_B = 2,
		B_is_subset_of_A = 3,
		is_intersection  = 4
	};

	Output_container  _result;  // std::list<Object> _result;
};



template<class I>
class Bops_Simple_Polygons_2 : public Bops_Polygons_2<I> {
public:

  typedef Bops_dcel<I>                			Dcel;
  typedef typename Dcel::const_faces_iterator       face_iterator;
  typedef typename Dcel::const_edges_iterator       edge_iterator;
  typedef typename Dcel::const_vertices_iterator    vertex_iterator;
      
  typedef _intersecting_polygons<R, Polygon_Container> Intersect_Polygons;

  bool operation() {
    if( !_init ) return false;
    create_dcel();
    unmark ();
    perform(); // call of child-object
    return empty();
  }


  Bops_Simple_Polygons_2() {
    marked_vector_init();
    _pgon_intersection_type= is_empty;
    _init= false;
  }


  bool do_intersect() const {
    return _init ? ( _pgon_intersection_type > 0) : false; }

  Bops_Simple_Polygons_2(const _Polygon_2& pgon1, const _Polygon_2& pgon2)
    : _pgon1(pgon1), _pgon2(pgon2)
  {
      _inter_res= Intersect_Polygons(pgon1,pgon2);
      if( _inter_res.size() > 0 )
          _pgon_intersection_type= is_intersection; // Intersection
      else
          _pgon_intersection_type= calc_intersection_type();

#   ifdef CGAL__BOPS_DEBUG_ON
    cout << "_inter_res.size()="  << _inter_res.size() << endl;
    cout << "intersection_type= " << _pgon_intersection_type << endl;
#   endif // CGAL__BOPS_DEBUG_ON

    marked_vector_init();
    _init= true;
  }
 
     
  virtual ~Bops_Simple_Polygons_2() {}

protected:

  virtual void perform(void) {};

  int calc_intersection_type(int) const {
     typename I::Bbox a_box= I::get_Bbox(_pgon1);
     typename I::Bbox b_box= I::get_Bbox(_pgon2);

     if( a_box == b_box ) {
       if( _pgon1 == _pgon2 ) return is_identical;
     }
     if( !I::do_overlap(a_box, b_box) ) return is_empty;

     if( I::box_is_contained_in_box( a_box, b_box) )
       return A_is_subset_of_B;
     if( I::box_is_contained_in_box( b_box, a_box) )
       return B_is_subset_of_A;

     typename Polygon_Container::const_iterator it;
     int sum1= 0, n1= _pgon1.size();
     for( it= _pgon1.vertices_begin(); it != _pgon1.vertices_end(); it++)
       sum1 += I::has_on_bounded_side(_pgon2, *it) ? +1 : -1;

     if( sum1 ==  n1 || sum1 == -n1 ) {
       int sum2= 0, n2= _pgon2.size();
       for( it= _pgon2.vertices_begin(); it != _pgon2.vertices_end(); it++)
         sum2 += I::has_on_bounded_side(_pgon1,*it) ? +1 : -1;
       if( sum2 ==  n2|| sum2 == -n2) {
         if( sum1 == -n1 && sum2 == -n2 ) // polygons are separated
           // iff no intersections occur, otherwise return 4
           return !_inter_res.size() ? is_empty : is_intersection;
         if( sum1 == n1 && sum2 == -n2 ) // A is subset B
           return A_is_subset_of_B; 
         if( sum1 == -n1 && sum2 == n2 ) // B is subset A
           return B_is_subset_of_A; 
       }
     }

     return is_intersection; // intersections occur
  }

  Intersection_type calc_intersection_type(void) const {
    return (Intersection_type)calc_intersection_type(0);
  }

  void create_dcel(void) {
    /* built up the graph (step 2 in README) */
    std::vector<_Point_2> ptlst;
    std::list< std::pair<int,int> > edlst;

    _inter_res.get_graph_information(ptlst,edlst);
#   ifdef CGAL__BOPS_DEBUG_ON
      cout << "after _inter_res.get_graph_information(ptlst,edlst);" << endl;
      print(cout, "ptlst:", ptlst.begin(), ptlst.end());
      print(cout, "edlst:", edlst.begin(), edlst.end());
      cout << flush;
#   endif //CGAL__BOPS_DEBUG_ON

    dcel.insert(edlst, ptlst);
#   ifdef CGAL__BOPS_DEBUG_ON
      cout << endl << "after dcel.insert(edlst, ptlst)" << endl;
      print(cout, "E", dcel.begin(), dcel.end() );
      print(cout, "V", dcel.vertex_begin(), dcel.vertex_end() );
      print(cout, "F", dcel.face_begin(), dcel.face_end() );
      cout << flush;
#   endif
    marked_vector_init();

    /* coloring the dcel */
    std::list<_Point_2> pts_on_A;
    std::list<_Point_2> pts_on_B;


    pts_on_A = _inter_res.get_color_informationA();
    dcel.colorize(pts_on_A, _RED);
    pts_on_B = _inter_res.get_color_informationB();
    dcel.colorize(pts_on_B, _BLACK);


#   ifdef CGAL__BOPS_DEBUG_ON
      cout << endl << "after dcel.colorize()" << endl;
      print(cout, "E", dcel.begin(), dcel.end() );
      print(cout, "V", dcel.vertex_begin(), dcel.vertex_end() );
      print(cout, "F", dcel.face_begin(), dcel.face_end() );
      cout << flush;
#   endif
  }

  int walk_around(face_iterator face, std::list<_Point_2>& result) const
  {
    /*
       traverses a face in the DCEL and puts the vertices as points in
       into a list. These points define a polygon.
       Additionally, the visited edges will be marked in the marked-vector.
       The size of the marked-vector should be at least the number of edges
       in the DCEL, i.e. marked.size() >= dcel.number_of_edges
       (The result is always be ordered clockwise.)
    */
  
    int n= 0;
    vertex_iterator v0, v; 
    edge_iterator e, a= dcel.begin(face);
    mark(a);
    v0= ((*a).F1() == face) ? (*a).V1() : (*a).V2();
    result.push_back(dcel.point(v0)); n++;
    mark(v0);
    for( e= dcel.next(a,face); e != a; e= dcel.next(e, face) ) {
      mark(e);
      v= (*e).opposite_vertex( v0 );
      mark(v);
      result.push_back(dcel.point(v)); n++;
      v0= v;
    }
  
    return n;
  }


  _Polygon_2 walk_around(face_iterator face, bool ccw = true) const {
    /*
       traverses a face in the DCEL and puts the vertices as points in
       into a polygon.
       Additionally, the visited edges will be marked in the marked-vector.
       The size of the marked-vector should be at least the number of edges
       in the DCEL, i.e. marked.size() >= dcel.number_of_edges
       (The result is always be ordered clockwise.)
    */
    std::list<_Point_2> pt_list;
    if( walk_around( face, pt_list) )
      return ccw ? _Polygon_2(pt_list.rbegin(), pt_list.rend()) :
                   _Polygon_2(pt_list.begin(),  pt_list.end()) ;
  
    return _Polygon_2();
  }

  void marked_vector_init() {
    marked= std::vector<bool>(dcel.number_of_edges());
    marked_vertex= std::vector<bool>(dcel.number_of_vertices());
  }

  void unmark () {
    std::fill(marked.begin(), marked.end(), false);
    std::fill(marked_vertex.begin(), marked_vertex.end(), false);
  }

  void mark   (vertex_iterator v) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked_vertex[(*v).index()]= true;
  }
  void unmark (vertex_iterator v) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked_vertex[(*v).index()]= false;
  }

  bool is_unmarked (vertex_iterator v) const {
    return marked_vertex[(*v).index()] == false;
  }
  bool is_marked (vertex_iterator v) const { return !is_unmarked(v); }

  bool is_unmarked (edge_iterator e) const {
    return marked[(*e).index()] == false;
  }
  bool is_marked (edge_iterator e) const { return !is_unmarked(e); }

  void mark   (edge_iterator e) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked[(*e).index()]= true;
  }
  void unmark (edge_iterator e) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked[(*e).index()]= false;
  }

  void mark   (int index) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked[index]= true;
  }
  void unmark (int index) const {
    const_cast<Bops_Simple_Polygons_2<I>*>(this)->marked[index]= false;
  }

  bool _init;
  _Polygon_2 _pgon1, _pgon2;
  Intersection_type _pgon_intersection_type;
  Dcel dcel;
  std::vector<bool> marked;        // marked list for edges
  std::vector<bool> marked_vertex; // marked list for vertices

  Intersect_Polygons _inter_res;
};


template<class I>
struct Bops_Simple_Polygons_2_Intersection
       : public Bops_Simple_Polygons_2<I>
{
  typedef typename Bops_Simple_Polygons_2<I>::_Polygon_2 _Polygon_2;
  Bops_Simple_Polygons_2_Intersection() {}

  Bops_Simple_Polygons_2_Intersection(
       const _Polygon_2& pgon1, const _Polygon_2& pgon2)
       : Bops_Simple_Polygons_2<I>( pgon1, pgon2) {
  }

  void perform(void);
};


template<class I>
struct Bops_Simple_Polygons_2_Difference
       : public Bops_Simple_Polygons_2<I>
{
  typedef typename Bops_Simple_Polygons_2<I>::_Polygon_2 _Polygon_2;
  Bops_Simple_Polygons_2_Difference() {}

  Bops_Simple_Polygons_2_Difference(
       const _Polygon_2& pgon1, const _Polygon_2& pgon2)
       : Bops_Simple_Polygons_2<I>( pgon1, pgon2) {
  }

  void perform(void);
};


template<class I>
struct Bops_Simple_Polygons_2_Union
       : public Bops_Simple_Polygons_2<I>
{
  typedef typename Bops_Simple_Polygons_2<I>::_Polygon_2 _Polygon_2;
  Bops_Simple_Polygons_2_Union() {}

  Bops_Simple_Polygons_2_Union(
       const _Polygon_2& pgon1, const _Polygon_2& pgon2)
       : Bops_Simple_Polygons_2<I>( pgon1, pgon2) {
  }

  void perform(void);
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_simple_polygons_2.C>
#endif

#endif /* CGAL_BOPS_SIMPLE_POLYGONS_2_H */
