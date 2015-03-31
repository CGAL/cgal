#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H

#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/spatial_sort.h>

// Needed by remove to fill the hole.
#include <CGAL/Periodic_3_regular_triangulation_remove_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <boost/unordered_set.hpp>


namespace CGAL
{
template < class Gt,
		   class Tds = Triangulation_data_structure_3 < Triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3<> >,
		                                                Regular_triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3<> > >
		 >
class Periodic_3_Regular_triangulation_3 : public Periodic_3_triangulation_3<Gt,Tds>
{
	typedef Periodic_3_Regular_triangulation_3<Gt,Tds>  Self;

public:
	typedef Periodic_3_triangulation_3<Gt,Tds>  Base;

	typedef Gt                                    Geometric_traits;
	typedef Tds                                   Triangulation_data_structure;

	typedef Geometric_traits                      Geom_traits;
	typedef typename Gt::FT                       FT;

	typedef typename Gt::Weighted_point_3         Weighted_point;
	typedef typename Gt::Bare_point               Bare_point;
	typedef typename Gt::Segment_3                Segment;
	typedef typename Gt::Triangle_3               Triangle;
	typedef typename Gt::Tetrahedron_3            Tetrahedron;

	typedef typename Base::Periodic_point         Periodic_point;
	typedef typename Base::Periodic_segment       Periodic_segment;
	typedef typename Base::Periodic_triangle      Periodic_triangle;
	typedef typename Base::Periodic_tetrahedron   Periodic_tetrahedron;

	typedef typename Base::Cell_handle            Cell_handle;
	typedef typename Base::Vertex_handle          Vertex_handle;

	typedef typename Base::Cell                   Cell;
	typedef typename Base::Vertex                 Vertex;
	typedef typename Base::Facet                  Facet;
	typedef typename Base::Edge                   Edge;

	typedef typename Base::Cell_circulator        Cell_circulator;
	typedef typename Base::Facet_circulator       Facet_circulator;
	typedef typename Base::Cell_iterator          Cell_iterator;
	typedef typename Base::Facet_iterator         Facet_iterator;
	typedef typename Base::Edge_iterator          Edge_iterator;
	typedef typename Base::Vertex_iterator        Vertex_iterator;

	typedef typename Base::All_cells_iterator     All_cells_iterator;
	typedef typename Base::All_facets_iterator    All_facets_iterator;
	typedef typename Base::All_edges_iterator     All_edges_iterator;
	typedef typename Base::All_vertices_iterator  All_vertices_iterator;

  typedef typename Base::size_type              size_type;
	typedef typename Base::difference_type        difference_type;

	typedef typename Base::Locate_type            Locate_type;
	typedef typename Base::Iterator_type          Iterator_type;

	typedef typename Base::Offset                 Offset;
	typedef typename Base::Iso_cuboid             Iso_cuboid;
	typedef typename Base::Covering_sheets        Covering_sheets;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
using Base::cw;
using Base::ccw;
using Base::domain;
using Base::geom_traits;
using Base::int_to_off;
using Base::number_of_sheets;
using Base::number_of_vertices;
using Base::number_of_edges;
using Base::number_of_facets;
using Base::number_of_cells;
using Base::cells_begin;
using Base::cells_end;
using Base::vertices_begin;
using Base::vertices_end;
using Base::facets_begin;
using Base::facets_end;
using Base::tds;
using Base::next_around_edge;
using Base::vertex_triple_index;
using Base::mirror_vertex;
using Base::orientation;
using Base::insert_dummy_points;
using Base::swap;
using Base::is_1_cover;
using Base::is_virtual;
using Base::point;
#endif

// For strict-ansi compliance
using Base::adjacent_vertices;
using Base::combine_offsets;
using Base::get_offset;
using Base::get_original_vertex;
using Base::get_location_offset;
using Base::get_neighbor_offset;
using Base::incident_edges;
using Base::incident_facets;
using Base::incident_cells;
using Base::is_valid_conflict;
using Base::locate;
using Base::periodic_point;
using Base::segment;

private:
  struct Cell_handle_hash : public std::unary_function<Cell_handle, std::size_t>
  {
    std::size_t operator()(const Cell_handle& ch) const { return boost::hash<typename Cell_handle::pointer>()(&*ch); }
  };
  boost::unordered_set<Cell_handle, Cell_handle_hash> cells_with_too_big_orthoball;

  class Cover_manager
  {
    Periodic_3_Regular_triangulation_3& tr;

  public:
    Cover_manager (Periodic_3_Regular_triangulation_3& tr)
  : tr(tr)
  {}

    void create_initial_triangulation()
    {
      tr.create_initial_triangulation();
    }

    template <class CellIt>
    void delete_too_long_edges(const CellIt begin, const CellIt end)
    {
      tr.delete_too_long_edges(begin, end);
    }

    template <class CellIt>
    void insert_too_long_edges(Vertex_handle v, const CellIt begin, const CellIt end)
    {
      tr.insert_too_long_edges(v, begin, end);
    }

    bool can_be_converted_to_1_sheet () const
    {
      return tr.can_be_converted_to_1_sheet();
    }

    bool update_cover_data_during_management (Cell_handle new_ch, const std::vector<Cell_handle>& new_cells)
    {
      return tr.update_cover_data_during_management(new_ch, new_cells);
    }
  };

public:
	/** @name Creation */ //@{
	Periodic_3_Regular_triangulation_3 (const Iso_cuboid& domain = Iso_cuboid(0, 0, 0, 1, 1, 1),
			                            const Geometric_traits& gt = Geometric_traits())
	: Base(domain, gt)
	{
	}

  // copy constructor duplicates vertices and cells
  Periodic_3_Regular_triangulation_3 (const Periodic_3_Regular_triangulation_3& tr)
  : Base(tr)
  {
    CGAL_triangulation_expensive_postcondition( is_valid() );
  }

  void create_initial_triangulation()
  {
    CGAL_triangulation_assertion( cells_with_too_big_orthoball.empty() );

    for (Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
      cells_with_too_big_orthoball.insert(iter);
  }

  template <class CellIt>
  void delete_too_long_edges(CellIt begin, const CellIt end)
  {
    for (; begin != end; ++begin)
    {
      typename boost::unordered_set<Cell_handle>::iterator iter = cells_with_too_big_orthoball.find(*begin);
      if (iter != cells_with_too_big_orthoball.end())
      {
        cells_with_too_big_orthoball.erase(iter);
      }
    }
  }

  FT squared_orthoball_radius (Cell_handle cell)
  {
    Periodic_point p0 = periodic_point(cell, 0);
    Periodic_point p1 = periodic_point(cell, 1);
    Periodic_point p2 = periodic_point(cell, 2);
    Periodic_point p3 = periodic_point(cell, 3);

    typename Geometric_traits::Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3
                                          = geom_traits().construct_weighted_circumcenter_3_object();

    Bare_point weighted_circumcenter = construct_weighted_circumcenter_3(
        p0.first,  p1.first,  p2.first,  p3.first,
        p0.second, p1.second, p2.second, p3.second);
    Weighted_point pt = point(p0);
    FT ao_2 = squared_distance(static_cast<const Bare_point&>(pt), weighted_circumcenter);
    FT io_2 = ao_2 - pt.weight();
    return io_2;
  }

  template <class CellIt>
  void insert_too_long_edges(Vertex_handle /*v*/, CellIt begin, const CellIt end)
  {
    FT threshold = FT(1)/FT(64) * (domain().xmax()-domain().xmin());
    for (; begin != end; ++begin)
    {
      if (squared_orthoball_radius(*begin) >= threshold)
      {
        cells_with_too_big_orthoball.insert(*begin);
      }
    }
  }

  bool can_be_converted_to_1_sheet () const
  {
    return cells_with_too_big_orthoball.empty();
  }

  bool update_cover_data_during_management (Cell_handle new_ch, const std::vector<Cell_handle>& new_cells)
  {
    bool result = false;
    FT threshold = FT(1)/FT(64) * (domain().xmax() - domain().xmin());

    if (squared_orthoball_radius(new_ch) >= threshold)
    {
      if (is_1_cover())
      {
        tds().delete_cells(new_cells.begin(), new_cells.end());
        this->convert_to_27_sheeted_covering();
        result = true;
      }
      else
        cells_with_too_big_orthoball.insert(new_ch);
    }

    return result;
  }

  virtual void update_cover_data_after_converting_to_27_sheeted_covering ()
  {
    FT threshold = FT(1)/FT(64) * (domain().xmax()-domain().xmin());
    for (Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
    {
      if (squared_orthoball_radius(iter) >= threshold)
      {
        cells_with_too_big_orthoball.insert(iter);
      }
    }
  }

  virtual void update_cover_data_after_setting_domain () {}

  virtual void reinsert_hidden_points_after_converting_to_1_sheeted (std::vector<Weighted_point>& hidden_points)
  {
    while (hidden_points.size())
    {
      insert(hidden_points.back());
      hidden_points.pop_back();
    }
  }

  /** @name Insertion */ //@{
   Vertex_handle insert(const Weighted_point& p, Cell_handle start = Cell_handle()) {
     Conflict_tester tester(p, this);
     Point_hider hider(this);
     Cover_manager cover_manager(*this);
     assert(p.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ));
     return Base::insert_in_conflict(p, start, tester, hider, cover_manager);
   }

   Vertex_handle insert(const Weighted_point& p, Locate_type lt, Cell_handle c,
        int li, int lj) {
      Conflict_tester tester(p, this);
      Point_hider hider(this);
      Cover_manager cover_manager(*this);
      assert(p.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ));
      return Base::insert_in_conflict(p,lt,c,li,lj, tester,hider,cover_manager);
    }
   //@}

   void remove(Vertex_handle v)
   {
     typedef CGAL::Periodic_3_regular_triangulation_remove_traits_3< Gt > P3removeT;
     typedef CGAL::Regular_triangulation_3< P3removeT >
       Euclidean_triangulation;
     typedef Vertex_remover< Euclidean_triangulation > Remover;
     P3removeT remove_traits(domain());
     Euclidean_triangulation tmp(remove_traits);
     Remover remover(this, tmp);
     Conflict_tester ct(this);
     Cover_manager cover_manager(*this);

     Base::remove(v, remover, ct, cover_manager);

     // Re-insert the points that v was hiding.
     for (typename Remover::Hidden_points_iterator hi = remover.hidden_points_begin();
          hi != remover.hidden_points_end();
          ++hi)
     {
         insert(*hi);
     }

     CGAL_triangulation_expensive_assertion(is_valid());
   }

protected:
	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
	}

	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r,
			                  const Offset &o1, const Offset &o2, const Offset &o3)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r, o1, o2, o3) == SMALLER;
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s,
			                                    const Offset& o1, const Offset& o2, const Offset& o3, const Offset& o4) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s, o1,o2,o3,o4);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r,
			                                   const Offset& o1, const Offset& o2, const Offset& o3) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}

//protected:
public:

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q) const
  {
      CGAL_triangulation_precondition(this->equal(p, q));
      return geom_traits().power_test_3_object()(p, q);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
       const Weighted_point &r, const Weighted_point &s,
       const Weighted_point &t, const Offset &o_p,
       const Offset &o_q, const Offset &o_r, const Offset &o_s,
       const Offset &o_t) const
  {
      return geom_traits().power_test_3_object()(p, q, r, s, t, o_p, o_q, o_r, o_s, o_t);
  }

  Oriented_side side_of_oriented_power_sphere(const Weighted_point &p, const Weighted_point &q,
      const Weighted_point &r, const Weighted_point &s, const Weighted_point &t, const Offset &o_p,
      const Offset &o_q, const Offset &o_r, const Offset &o_s,
      const Offset &o_t) const
  {
    return power_test(p,q,r,s,t,o_p,o_q,o_r,o_s,o_t);
  }

  Bounded_side _side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
      const Offset & offset = Offset(), bool perturb = false) const;

  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

protected:
  // Protected, because inheritors(e.g. periodic triangulation for meshing)
  // of the class Periodic_3_Delaunay_triangulation_3 use this class
  class Conflict_tester;
private:
  class Point_hider;

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  template <class TriangulationR3> struct Vertex_remover;
#else
  template <class TriangulationR3>
  struct Vertex_remover
  {
    typedef TriangulationR3      Triangulation_R3;

    typedef typename std::vector<Weighted_point>::iterator Hidden_points_iterator;

    typedef Triple < Vertex_handle, Vertex_handle, Vertex_handle > Vertex_triple;

    typedef typename Triangulation_R3::Triangulation_data_structure TDSE;
    typedef typename Triangulation_R3::Cell_handle        CellE_handle;
    typedef typename Triangulation_R3::Vertex_handle      VertexE_handle;
    typedef typename Triangulation_R3::Facet              FacetE;
    typedef typename Triangulation_R3::Finite_cells_iterator Finite_cellsE_iterator;

    typedef Triple< VertexE_handle, VertexE_handle, VertexE_handle >
    VertexE_triple;

    typedef std::map<Vertex_triple,Facet> Vertex_triple_Facet_map;
    typedef std::map<Vertex_triple, FacetE> Vertex_triple_FacetE_map;
    typedef typename Vertex_triple_FacetE_map::iterator
    Vertex_triple_FacetE_map_it;

    Vertex_remover(const Self *t, Triangulation_R3 &tmp_) : _t(t),tmp(tmp_) {}

    const Self *_t;
    Triangulation_R3 &tmp;

    void add_hidden_points(Cell_handle ch) {
      std::copy(ch->hidden_points_begin(), ch->hidden_points_end(),
    std::back_inserter(hidden));
    }

    Hidden_points_iterator hidden_points_begin() {
      return hidden.begin();
    }
    Hidden_points_iterator hidden_points_end() {
      return hidden.end();
    }
    //private:
    // The removal of v may un-hide some points,
    // Space functions output them.
    std::vector<Weighted_point> hidden;
  };
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
};

template < class Gt, class Tds >
Bounded_side Periodic_3_Regular_triangulation_3<Gt,Tds>::
_side_of_power_sphere(const Cell_handle &c, const Weighted_point &q,
    const Offset &offset, bool perturb ) const
{
  Weighted_point p0 = c->vertex(0)->point(),
        p1 = c->vertex(1)->point(),
        p2 = c->vertex(2)->point(),
        p3 = c->vertex(3)->point();
  Offset o0 = this->get_offset(c,0),
        o1 = this->get_offset(c,1),
        o2 = this->get_offset(c,2),
        o3 = this->get_offset(c,3),
        oq = offset;

  CGAL_triangulation_precondition( orientation(p0, p1, p2, p3, o0, o1, o2, o3) == POSITIVE );

  Oriented_side os = ON_NEGATIVE_SIDE;
  os= side_of_oriented_power_sphere(p0, p1, p2, p3, q, o0, o1, o2, o3, oq);

  if (os != ON_ORIENTED_BOUNDARY || !perturb)
    return (Bounded_side) os;

  //We are now in a degenerate case => we do a symbolic perturbation.
  // We sort the points lexicographically.
  Periodic_point pts[5] = {std::make_pair(p0,o0), std::make_pair(p1,o1),
         std::make_pair(p2,o2), std::make_pair(p3,o3),
         std::make_pair(q,oq)};
  const Periodic_point *points[5] ={&pts[0],&pts[1],&pts[2],&pts[3],&pts[4]};

  std::sort(points, points+5,
      typename Base::template Perturbation_order< typename Gt::Compare_xyz_3 >(geom_traits().compare_xyz_3_object()));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  for (int i=4; i>1; --i) {
    if (points[i] == &pts[4]) {
      CGAL_triangulation_assertion(orientation(p0, p1, p2, p3, o0, o1, o2, o3)
          == POSITIVE);
      // since p0 p1 p2 p3 are non coplanar and positively oriented
      return ON_UNBOUNDED_SIDE;
    }
    Orientation o;
    if (points[i] == &pts[3] &&
        (o = orientation(p0, p1, p2, q, o0, o1, o2, oq)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[2] &&
        (o = orientation(p0, p1, q, p3, o0, o1, oq, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[1] &&
        (o = orientation(p0, q, p2, p3, o0, oq, o2, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[0] &&
        (o = orientation(q, p1, p2 ,p3, oq, o1, o2, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
  }

  CGAL_triangulation_assertion(false);
  return ON_UNBOUNDED_SIDE;
}

template < class Gt, class Tds >
bool
Periodic_3_Regular_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if (!Base::is_valid(verbose, level)) {
    if (verbose)
      std::cerr << "Regular: invalid base" << std::endl;
    return false;
  }

  Conflict_tester tester(this);
  if (!is_valid_conflict(tester, verbose, level)) {
    if (verbose)
      std::cerr << "Regular: conflict problems" << std::endl;
    return false;
  }

  if (verbose)
    std::cerr << "Regular valid triangulation" << std::endl;
  return true;
}

template < class GT, class TDS >
bool
Periodic_3_Regular_triangulation_3<GT,TDS>::
is_valid(Cell_handle ch, bool verbose, int level) const {
  bool error = false;
  if (!Base::is_valid(ch, verbose, level)) {
    error = true;
    if (verbose) {
      std::cerr << "geometrically invalid cell" << std::endl;
      for (int i=0; i<4; i++ )
  std::cerr << ch->vertex(i)->point() << ", ";
      std::cerr << std::endl;
    }
  }
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++ vit) {
    for (int i=-1; i<=1; i++)
      for (int j=-1; j<=1; j++)
  for (int k=-1; k<=1; k++) {
    if (periodic_point(ch,0) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
    || periodic_point(ch,1) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
          || periodic_point(ch,2) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
    || periodic_point(ch,3) == std::make_pair(periodic_point(vit).first,
                  periodic_point(vit).second+Offset(i,j,k)) )
      continue;
    if (_side_of_power_sphere(ch, periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k),true)
        != ON_UNBOUNDED_SIDE) {
      error = true;
      if (verbose) {
        std::cerr << "Regular invalid cell" << std::endl;
        for (int i=0; i<4; i++ ) {
    Periodic_point pp = periodic_point(ch,i);
    std::cerr <<"("<<pp.first <<","<<pp.second<< "), ";
        }
        std::cerr << std::endl;
      }
    }
  }
  }
  return !error;
}

template < class GT, class Tds >
class Periodic_3_Regular_triangulation_3<GT,Tds>::Conflict_tester
{
  // stores a pointer to the triangulation,
  // a point, and an offset
  const Self *t;
  Weighted_point p;
  // stores the offset of a point in 27-cover
  mutable Offset o;

public:
  /// Constructor
  Conflict_tester(const Self *_t) : t(_t), p(Weighted_point()) {}
  Conflict_tester(const Weighted_point &pt, const Self *_t) : t(_t), p(pt) { }

  /** The functor
    *
    * gives true if the circumcircle of c contains p
    */
  bool operator()(const Cell_handle c, const Offset &off) const {
    return (t->_side_of_power_sphere(c, p, t->combine_offsets(o, off), true)
        == ON_BOUNDED_SIDE);
  }

  bool operator()(const Cell_handle c, const Weighted_point& pt,
      const Offset &off) const {
    return (t->_side_of_power_sphere(c, pt, o + off, true) == ON_BOUNDED_SIDE);
  }

  int compare_weight(const Weighted_point& p, const Weighted_point& q) const
  {
    return t->power_test(p, q);
  }

  bool test_initial_cell(Cell_handle c, const Offset &off) const
  {
    return (operator()(c, off));
  }

  void set_point(const Weighted_point &_p) {
    p = _p;
  }

  void set_offset(const Offset &off) const {
    o = off;
  }

  const Offset &get_offset() const {
    return o;
  }

  const Weighted_point &point() const {
    return p;
  }

};

template < class GT, class Tds>
class Periodic_3_Regular_triangulation_3<GT,Tds>::Point_hider
{
  Self *t;
  mutable std::vector<Vertex_handle> vertices;
  mutable std::vector<Weighted_point> hidden_points;
  mutable bool is_original_cube;

public:
  Point_hider(Self *tr) : t(tr), is_original_cube(false) {}

  void set_original_cube (bool b) const {
    is_original_cube = b;
  }

  template <class InputIterator>
  inline void set_vertices(InputIterator start, InputIterator end) const
  {
    while (start != end) {
        std::copy((*start)->hidden_points_begin(),
            (*start)->hidden_points_end(),
            std::back_inserter(hidden_points));

      for (int i=0; i<=3; i++) {
        Vertex_handle v = (*start)->vertex(i);
        if (v->cell() != Cell_handle()) {
          vertices.push_back(v);
          v->set_cell(Cell_handle());
        }
      }
      start ++;
    }
  }

  inline void reinsert_vertices(Vertex_handle v)
  {
    Locate_type lt = Locate_type();
    int li=0, lj=0;

    Cell_handle hc = v->cell();
    for (typename std::vector<Vertex_handle>::iterator
        vi = vertices.begin(); vi != vertices.end(); ++vi) {
      if ((*vi)->cell() != Cell_handle()) continue;
      if (is_original_cube)
      {
        hc = t->locate((*vi)->point(), lt, li, lj, hc);
        hc->hide_point((*vi)->point());
      }
      t->delete_vertex(*vi);
    }
    vertices.clear();
      for (typename std::vector<Weighted_point>::iterator
          hp = hidden_points.begin(); hp != hidden_points.end(); ++hp) {
        hc = t->locate(*hp, lt, li, lj, hc);
        hc->hide_point(*hp);
      }
      hidden_points.clear();
  }

  inline Vertex_handle replace_vertex(Cell_handle c, int index, const Weighted_point& p)
  {
    Vertex_handle v = c->vertex(index);
    c->hide_point(v->point());
    v->set_point(p);
    return v;
  }

  inline void hide_point(Cell_handle c, const Weighted_point& p)
  {
    if (is_original_cube)
      c->hide_point(p);
  }

//  inline void hide(Weighted_point&, Cell_handle ) const  // useless?
//  {
//    CGAL_triangulation_assertion(false);
//  }
//
//  inline void do_hide(const Weighted_point&, Cell_handle ) const // useless?
//  {
//    CGAL_triangulation_assertion(false);
//  }

//  template < class Tester >
//  inline bool replace_vertex(const Weighted_point&, Vertex_handle, const Tester&) const // useless?
//  {
//    return true;
//  }
//
//  template <class Conflict_tester>
//  inline void hide_points(Vertex_handle,
//      const Conflict_tester &)
//  {
//    // No points to hide in the Delaunay triangulation.
//  }
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class GT, class Tds>
template <class TriangulationR3>
struct Periodic_3_Regular_triangulation_3<GT,Tds>::Vertex_remover
{
  typedef TriangulationR3      Triangulation_R3;

  typedef typename std::vector<Weighted_point>::iterator Hidden_points_iterator;

  typedef Triple < Vertex_handle, Vertex_handle, Vertex_handle > Vertex_triple;

  typedef typename Triangulation_R3::Triangulation_data_structure TDSE;
  typedef typename Triangulation_R3::Cell_handle        CellE_handle;
  typedef typename Triangulation_R3::Vertex_handle      VertexE_handle;
  typedef typename Triangulation_R3::Facet              FacetE;
  typedef typename Triangulation_R3::Finite_cells_iterator Finite_cellsE_iterator;

  typedef Triple< VertexE_handle, VertexE_handle, VertexE_handle >
  VertexE_triple;

  typedef std::map<Vertex_triple,Facet> Vertex_triple_Facet_map;
  typedef std::map<Vertex_triple, FacetE> Vertex_triple_FacetE_map;
  typedef typename Vertex_triple_FacetE_map::iterator
  Vertex_triple_FacetE_map_it;

  Vertex_remover(const Self *t, Triangulation_R3 &tmp_) : _t(t),tmp(tmp_) {}

  const Self *_t;
  Triangulation_R3 &tmp;

  void add_hidden_points(Cell_handle ch) {
    std::copy(ch->hidden_points_begin(), ch->hidden_points_end(),
  std::back_inserter(hidden));
  }

  Hidden_points_iterator hidden_points_begin() {
    return hidden.begin();
  }
  Hidden_points_iterator hidden_points_end() {
    return hidden.end();
  }
  private:
  // The removal of v may un-hide some points,
  // Space functions output them.
  std::vector<Weighted_point> hidden;
};
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
}// namespace CGAL

#endif
