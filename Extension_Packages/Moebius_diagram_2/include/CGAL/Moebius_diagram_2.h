#ifndef CHR_MOEBIUS_DIAGRAM_2_H
#define CHR_MOEBIUS_DIAGRAM_2_H


#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_3.h>
//#include <CGAL/Conic_2.h>
#include <CGAL/Object.h>
//#include <CGAL/Line_2_Line_2_intersection.h>
#include <CGAL/Moebius_utils.h>


#include <list>

#ifdef CGAL_USE_QT
# include <CGAL/IO/Qt_widget.h>
#endif

#include <CGAL/Moebius_edge_2.h>
#include <CGAL/Moebius_circulator_3.h>

CGAL_BEGIN_NAMESPACE

template <class Gt, class Cb = Triangulation_cell_base_3<Gt> >
class Moebius_cell_base_3 : public Cb
{
  signed char _side;
  char _inter_facet[4];
  char _inter_edge[6];

 public:
  void init () {
    _side = 0;
    for (int i = 0; i < 4; ++i) _inter_facet[i] = (char) 0;
    for (int i = 0; i < 6; ++i) _inter_edge[i] = (char) 0;
  }

  typedef typename Cb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle Vertex_handle;
  typedef typename Tds::Cell_handle Cell_handle;

  template < typename TDS2 >
    struct Rebind_TDS {
      typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
      typedef Moebius_cell_base_3<Gt, Cb2>         Other;
    };



  Moebius_cell_base_3 () 
    : Cb (), _side (), _inter_facet (), _inter_edge() 
  { init (); }

  Moebius_cell_base_3 (Vertex_handle v0, Vertex_handle v1,
		      Vertex_handle v2, Vertex_handle v3)
    : Cb (v0, v1, v2, v3), _side (), _inter_facet (), _inter_edge() 
  { init (); }

  Moebius_cell_base_3(Vertex_handle v0, Vertex_handle v1,
		     Vertex_handle v2, Vertex_handle v3,
		     Cell_handle n0, Cell_handle n1,
		     Cell_handle n2, Cell_handle n3)
    : Cb (v0, v1, v2, v3, n0, n1, n2, n3),
    _side (), _inter_facet (), _inter_edge() { init (); }

  void set_side (const CGAL::Bounded_side &s) {
    _side = (((signed char) s) << 1) | 1;
  }
  
  bool side_is_set () const { 
    return (_side & 1) != 0;
  }

  CGAL::Bounded_side side () const {
    CGAL_precondition (side_is_set ());
    return CGAL::Bounded_side (_side >> 1);
  }

  void set_inter_facet (int f, int n) 
  {
    CGAL_precondition (0 <= f && f < 4);
    _inter_facet[f] = (char) n;
  }

  void set_inter_facet (Vertex_handle f, int n) 
  { set_inter_facet (index (f), n); }

  void set_inter_edge (int e1, int e2, int inter) 
  { _inter_edge[edge_index (e1, e2)] = (char) inter; }

  void set_inter_edge (Vertex_handle e1, Vertex_handle e2, int inter) 
  { set_inter_edge (index (e1), index (e2), inter); }

  int inter_facet (int f)
  {
    CGAL_precondition (0 <= f && f < 4);
    return _inter_facet[f];
  }

  int inter_facet (Vertex_handle v)
  { return inter_facet (index (v)); }

  int inter_edge (int e1, int e2)
  { return _inter_edge[edge_index (e1, e2)]; }

  int inter_edge (Vertex_handle e1, Vertex_handle e2)
  { return inter_edge (index (e1), index (e2)); }
private:
  static int edge_index (int i, int j)
  {
    static const int tab[4][4] = { {-1, 0, 1, 2 },
				   { 0,-1, 3, 4 },
				   { 1, 3,-1, 5 },
				   { 2, 4, 5,-1 } };
    CGAL_precondition (0 <= i && i < 4);
    CGAL_precondition (0 <= j && j < 4);
    CGAL_precondition (i != j);
      
    return tab[i][j];
  }

};

//#define CHR_HIERARCHY

#ifdef CHR_HIERARCHY
template <class Tr>
class Moebius_no_hierarchy : public Tr {};
#endif

#ifdef CHR_HIERARCHY
template <class Gt, template <class T> class H = Moebius_no_hierarchy>
#else
template <class Gt>
#endif
class Moebius_diagram_2
{
 public:
  typedef typename Gt::Regular_traits Regular_traits;

  typedef Triangulation_vertex_base_3<Regular_traits> Vb;
  typedef Moebius_cell_base_3<Regular_traits> Cb;
  typedef Triangulation_data_structure_3<Vb, Cb>  Regular_tds;

#ifdef CHR_HIERARCHY
  typedef Regular_triangulation_3<Regular_traits, Regular_tds> Base_RT_3;
  typedef H<Base_RT_3> RT_3;
#else
  typedef Regular_triangulation_3<Regular_traits, Regular_tds> RT_3;
#endif

  typedef typename Gt::Point_2 Point;
  typedef typename Gt::Bare_point_2 Bare_point;
  typedef typename Gt::Segment_2 Segment;
  typedef typename Gt::Line_2 Line;
  typedef typename Gt::Circle_2 Circle;

  typedef Moebius_circulator_3<Regular_tds> Circulator;
  typedef Moebius_edge_2<Gt> Moebius_edge;
  typedef std::list<Moebius_edge> Diagram;
  typedef typename Diagram::iterator Moebius_edge_iterator;
  typedef typename Diagram::const_iterator Moebius_edge_const_iterator;

  typedef typename RT_3::Vertex_handle Vertex_handle;
  typedef typename RT_3::Vertex_iterator Vertex_iterator;
  typedef typename RT_3::Cell_handle Cell_handle;
  typedef typename RT_3::Cell_circulator Cell_circulator;
  typedef typename RT_3::Edge Edge;
  typedef typename RT_3::Facet Facet;
  typedef typename RT_3::Facet_circulator Facet_circulator;
  typedef typename RT_3::Finite_edges_iterator Finite_edges_iterator;

 private:
  Gt _gt;
  RT_3 _rt;
  bool _initialized, _built, _constructed;
  Diagram _dia;

 public:
  Moebius_diagram_2 ()
    : _gt (), _rt (), _initialized (false), _built (false), 
    _constructed (false), _dia () 
  {}

  Moebius_diagram_2 (const Moebius_diagram_2<Gt> &base)
    : _gt(base.gt()), _rt(), _initialized(base.initialized()), 
    _built (base.built ()), _constructed (base.constructed ()), 
    _dia(base.diagram ())
  {     
    if (_initialized) _rt = base.rt ();
  }

  template <class InputIterator>
  Moebius_diagram_2 (const InputIterator &first, const InputIterator &last)
    : _gt (), _rt (), _initialized (false), _built (false), 
    _constructed (false), _dia ()
  {
    init (first, last);
    build ();
    construct ();
  }

  template <class InputIterator>
    int init (const InputIterator &first, const InputIterator &last)
    {
      CGAL_precondition (! _initialized && ! _built);
      InputIterator i = first;
      int n = 0;
      while (i != last) {
	_rt.insert (*i);

	++i;
	++n;
      }

      _initialized = true;
      return n;
    }

  void clear () {
    _rt.clear ();
    _dia.clear ();
    _initialized = _built = _constructed = false;
  }

  bool initialized () const { return _initialized; }
  bool built () const { return _built; }
  bool constructed () const { return _constructed; }
  const Gt &gt () const  { return _gt; }
  const RT_3 &rt () const { CGAL_precondition (_initialized); return _rt; }

  bool is_finite (const Vertex_handle &v) const 
  { 
    return ! _rt.is_infinite (v); 
}

  bool is_finite (const Edge &e) const { return ! _rt.is_infinite (e); }
  bool is_finite (const Facet &f) const { return ! _rt.is_infinite (f); }
  bool is_finite (const Cell_handle &c) const { return ! _rt.is_infinite (c); }


  bool inter_facet_is_set (const Facet &f) const {
    return f.first->inter_facet_is_set (f.second);
  }

  int inter_facet (const Facet &f) const
    {
      return f.first->inter_facet (f.second);
    }

  int inter_edge (const Edge &e) const
    {
      return e.first->inter_edge (e.second, e.third);
    }

  Finite_edges_iterator finite_edges_begin () const { 
    CGAL_precondition (_built);
    return rt().finite_edges_begin ();
  }

  Finite_edges_iterator finite_edges_end () const { 
    CGAL_precondition (_built);
    return rt().finite_edges_end ();
  }

  Moebius_edge dual (const Edge &edge) {
    CGAL_assertion (built ());
    Moebius_edge me = dual_ (edge);
    me.finalize ();
    return me;
  }

  void construct () {
    Finite_edges_iterator i = finite_edges_begin (), end = finite_edges_end ();
    while (i != end) {
      if (inter_edge (*i) > 0)
	_dia.push_back (dual (*i));
      ++i;
    }
    _constructed = true;
  }

  const Diagram &diagram () const { return _dia; }

  Moebius_edge_iterator moebius_edges_begin () const {
    CGAL_precondition (constructed ());
    return _dia.begin ();
  }

  Moebius_edge_iterator moebius_edges_end () const {
    CGAL_precondition (constructed ());
    return _dia.end ();
  }

  Moebius_edge_const_iterator moebius_edges_const_begin () const {
    CGAL_precondition (constructed ());
    return _dia.begin ();
  }

  Moebius_edge_const_iterator moebius_edges_const_end () const {
    CGAL_precondition (constructed ());
    return _dia.end ();
  }

  void build () {
    CGAL_precondition (_initialized && ! _built);
    switch (dimension ()) {
    case -2: case -1: case 0:
      // nothing to do: no finite edge
      break;
    case 1: {
      TRACE ("build(): dimension 1\n");
      Finite_edges_iterator i = rt().finite_edges_begin ();
      Finite_edges_iterator end = rt().finite_edges_end ();
      for (; i != end; ++i) build_edge_1D (*i);
    } break;
    case 2: {
      TRACE ("build(): dimension 2\n");
      Finite_edges_iterator i = rt().finite_edges_begin ();
      Finite_edges_iterator end = rt().finite_edges_end ();
      for (; i != end; ++i) build_edge_2D (*i);
    } break;
    case 3: {
      TRACE ("build(): dimension 3\n");
      Finite_edges_iterator i = rt().finite_edges_begin ();
      Finite_edges_iterator end = rt().finite_edges_end ();
      for (; i != end; ++i) build_edge_3D (*i);
    } break;
    default:
      // oops
      CGAL_assertion (false);
    }
    _built = true;
  }
  

 private:

  int dimension () const { return rt().dimension (); }


  void build_edge_1D (const Edge &edge) {
    TRACE ("build_edge_1D\n");
    CGAL_assertion (dimension () == 1);

    Vertex_handle p = edge.first->vertex (edge.second);
    Vertex_handle q = edge.first->vertex (edge.third);
    
    if (moebius_orientation (p, q) == CGAL::ZERO) {
      TRACE ("  line\n");
      set_inter_edge (edge, 1);
    } else if (has_circle (p, q)) {
      TRACE ("  circle\n");
      set_inter_edge (edge, 2);
    } else {
      TRACE ("  no intersection\n");
      set_inter_edge (edge, 0);
    }
  }

  void build_edge_2D (const Edge &edge) {
    TRACE ("build_edge_2D...\n");
    CGAL_assertion (dimension () == 2);

    Vertex_handle p = edge.first->vertex (edge.second);
    Vertex_handle q = edge.first->vertex (edge.third);

    CGAL::Orientation o = moebius_orientation (p, q);

    int inter = o == CGAL::ZERO ? 1 : has_circle (p, q) ? 2 : 0;
    if (inter != 0) {
      
      TRACE (" finding the first facet...\n");
      int i;
      for (i = 0; i < 3; ++i) if (i != edge.second && i != edge.third) break;
      CGAL_assertion (i < 3);
      Vertex_handle r = edge.first->vertex (i);
      
      TRACE (" processing the first facet...\n");
      
      Facet f1 (edge.first, 3);
        if (is_finite (r)) {
	int facet = build_facet_2D (o, f1, p, q, r);
	if (facet < 0) {
	  set_inter_edge (edge, 0);
	  TRACE ("build_edge_2D: no intersection.\n");
	  return;
	}
      } 
      TRACE (" findind the other facet...\n");
      
      Cell_handle neighbor = edge.first->neighbor (i);
      int i1 = neighbor->index (edge.first->vertex (edge.second));
      int i2 = neighbor->index (edge.first->vertex (edge.third));
      
      for (i = 0; i < 3; ++i) if (i != i1 && i != i2) break;
      CGAL_assertion (i < 3);
      Vertex_handle s = neighbor->vertex (i);
      
      TRACE (" processing the other facet...\n");
      Facet f2 (neighbor, 3);
      if (is_finite (s)) {
	int facet = build_facet_2D (o, f2, p, q, s);
	if (facet < 0) {
	  set_inter_edge (edge, 0);
	  TRACE ("build_edge_2D: no intersection.\n");
	  return;
	}
      }
    }
    set_inter_edge (edge, inter);
    TRACE ("build_edge_2D: " << inter << "\n");
  }
  int build_facet_2D (const CGAL::Orientation &o,
		      const Facet &f,
		      const Vertex_handle &p,
		      const Vertex_handle &q,
		      const Vertex_handle &r
		      ) {
    int real = real_build_facet (o, f, p, q, r);
    TRACE ("  build_facet_2D: " << real << "\n");
    if (real == 0) {
      if (side_of_power_edge (o, p, q, r) != CGAL::ON_POSITIVE_SIDE)
	real = -1;
    }
    TRACE ("  real != 0\n");
    if (real < 0) {
      if (! effective_inter_facet_is_set (f))
	set_effective_inter_facet (f, 0);
      return -1;
    }
    TRACE ("  real > 0\n");
    if (effective_inter_facet_is_set (f))
      return effective_inter_facet (f);
    TRACE ("  effective not set\n");
    set_effective_inter_facet (f, real);
    return real;
  }
   

  void build_edge_3D (const Edge &edge) {
    TRACE ("build_edge_3D");
    CGAL_assertion (dimension () == 3);

    Cell_handle cell = edge.first;
    int a = edge.second;
    int b = edge.third;

    Vertex_handle p (cell->vertex(a));
    Vertex_handle q (cell->vertex(b));

    CGAL::Orientation o = moebius_orientation (p, q);

    TRACE (" ("<<p->point()<<") ("<<q->point()<<")\n");
    int inter = o == CGAL::ZERO ? 1 : has_circle (p, q) ? 2 : 0;
    if (inter != 0) {
      // loop over all the incident facet
      bool no_effective = true;
      bool some_real = false;

      Circulator start (edge), i (start);
      do {
	if (is_finite (i.vertex ())) {
	  TRACE ("  loop ("<<i.vertex()->point()<<")\n");
	  int real = real_build_facet (o, i.facet (), 
				       i.source(), i.target(), 
				       i.vertex());
	  if (real == 0) {
	    if (side_of_power_edge (o, i) != CGAL::ON_POSITIVE_SIDE) {
	      real = -1;
	    }
	  }
	  int effective = effective_build_facet (o, real, i);
	  TRACE ("    real " << real << ", eff " << effective << "\n");
	  some_real = some_real || (real > 0);
	  no_effective = no_effective && (effective <= 0);
	  if (real < 0) { 
	    // TODO: continue the loop and set real/effective_inter_facet to 0
	    // all around.
	    CGAL_assertion (no_effective);
	    inter = 0;
	    break;
	  }
	}

	++i;
      } while (i != start);
      if (no_effective && some_real) {
	inter = 0;
      }
    }
    set_inter_edge (edge, inter);
  }

  int effective_build_facet (const CGAL::Orientation &o, 
			     const int &real, 
			     const Circulator &i) {
    if (real < 0) {
      if (! effective_inter_facet_is_set (i.facet ()))
	set_effective_inter_facet (i.facet (), 0);
      return -1;
    }
    if (effective_inter_facet_is_set (i.facet ()))
      return effective_inter_facet (i.facet ());
    Vertex_handle s = i.source(), t = i.target();

    int inter = real;
    switch (real) {
    case 0: break;
    case 1:
      if (side_of_vertex (i.cell(), s, t) == side_of_vertex (i.next_cell(), s, t))
	inter = 0;
      break;
    case 2:
      if (side_of_vertex (i.cell(), s, t) != side_of_vertex (i.next_cell(), s, t))
	inter = 1;
      else if (side_of_vertex (i.cell(), s, t) == CGAL::ON_BOUNDED_SIDE)
	inter = 0;
      else if (same_side_of_vertex (o, i))
	inter = 0;
      break;
     default:
      CGAL_assertion (false);
    }
    set_effective_inter_facet (i.facet (), inter);
    return inter;
  }
  bool same_side_of_vertex (const CGAL::Orientation &o,
			    const Circulator &i) {
    CGAL_assertion (is_finite (i.vertex ()));
    
    if (o == CGAL::ZERO) {
      TRACE ("        same_side (swap)\n");
      return same_side_of_vertex (moebius_orientation (i.vertex(), i.source()),
				   i.vertex(), i.source(), i.prev_vertex(), i.target(), i.next_vertex());
    }
    TRACE ("        same_side (not swap)\n");
    return same_side_of_vertex (o, i.source(), i.target(), i.prev_vertex(), i.vertex(), i.next_vertex());
  }

  bool same_side_of_vertex (const CGAL::Orientation &o,
			    const Vertex_handle &p,
			    const Vertex_handle &q,
			    const Vertex_handle &prev,
			    const Vertex_handle &r,
			    const Vertex_handle &next) {
    CGAL_assertion (o != CGAL::ZERO);
    TRACE ("          orientation = " << o << "\n");
    if (is_finite (prev)) {
      if (is_finite (next)) {
	TRACE ("          finite, finite, calling the predicates\n");
	return circle_side_of_vertex (p, q, r, prev) != circle_side_of_vertex (p, q, r, next);
      }
      TRACE ("          finite, infinite, calling the first predicate\n");
      return circle_side_of_vertex (p, q, r, prev) == CGAL::Oriented_side (o);
    }
    if (is_finite (next)) {
      TRACE ("          infinite, finite, calling the second predicate\n");
      return circle_side_of_vertex (p, q, r, next) == CGAL::Oriented_side (o);
    }
    TRACE ("          infinite, infinite, no need of predicates.\n");
    return false;
  }

  CGAL::Oriented_side side_of_power_edge (const CGAL::Orientation &o, 
					  const Circulator &i) {
    return side_of_power_edge (o, i.source(), i.target(), i.vertex());
  }
  CGAL::Oriented_side side_of_power_edge (const CGAL::Orientation &o,
					  const Vertex_handle &p,
					  const Vertex_handle &q,
					  const Vertex_handle &r) {
    if (o == CGAL::ZERO)
      return CGAL::Oriented_side (moebius_orientation (r, q));
    return circle_side_of_center (p, q, r);
  }

  int real_build_facet (const CGAL::Orientation o, const Facet &f,
			const Vertex_handle &p,
			const Vertex_handle &q,
			const Vertex_handle &r) {
    if (real_inter_facet_is_set (f)) return real_inter_facet (f);
    int inter = real_build_facet (o, p, q, r);
    set_real_inter_facet (f, inter);
    return inter;
  }

  int real_build_facet (const Orientation o,
			const Vertex_handle &p,
			const Vertex_handle &q,
			const Vertex_handle &r) {
    if (o == CGAL::ZERO) {
      if (moebius_orientation (q, r) == CGAL::ZERO) return 1;
      if (circle_cross_line (q, r, p)) return 2;
      return 0;
    }
    
    if (circle_cross_line (p, q, r)) return 2;
    
    return 0;
  }


  // ------- dual ---------------

  Moebius_edge dual_ (const Edge &edge) {
    switch (dimension ()){
    case -1: case 0: return Moebius_edge ();
    case 1: return dual_1D (edge);
    case 2: return dual_2D (edge);
    case 3: return dual_3D (edge);
    default:
      CGAL_assertion (false);
      return Moebius_edge ();
    }
    CGAL_assertion (false);
    return Moebius_edge ();
  }
	

  Moebius_edge bisector (const Edge &edge) {
    Vertex_handle p = edge.first->vertex (edge.second);
    Vertex_handle q = edge.first->vertex (edge.third);
    
    if (moebius_orientation (p, q) == CGAL::ZERO) {
      return Moebius_edge (construct_line (p, q));
    } else if (has_circle (p, q)) {
      return Moebius_edge (construct_circle (p, q));
    } else {
      TRACE ("-!- oops\n");
      return Moebius_edge ();
    }
    CGAL_assertion (false);
    return Moebius_edge ();
  }

  Moebius_edge dual_1D (const Edge &edge) {
    CGAL_assertion (dimension () == 1);
    return bisector (edge);
  }

  Moebius_edge dual_2D (const Edge &edge) {
    CGAL_assertion (dimension () == 2);
    //    std::cerr << "dual_2D: not yet implemented\n";
    //CGAL_assertion (false);
    if (inter_edge (edge) == 0) return Moebius_edge ();

    Moebius_edge medge = bisector (edge);
    
    Vertex_handle p = edge.first->vertex (edge.second);
    Vertex_handle q = edge.first->vertex (edge.third);
    Vertex_handle inf = rt().infinite_vertex ();
    int i;
    for (i = 0; i < 3; ++i) if (i != edge.second && i != edge.third) break;
    CGAL_assertion (i < 3);
    Vertex_handle r = edge.first->vertex (i);
    
    Facet f1 (edge.first, 3);
    if (is_finite (r)) {
      //      const Edge other (f1.first,
      //		f1.first->index (p),
      //		f1.first->index (r));
      dual_facet (medge,
		  side_of_vertex (p, q, inf, r),
		  side_of_vertex (p, q, r, inf),
		  p, q, r, f1);//, bisector (other));
    }

    Cell_handle neighbor = edge.first->neighbor (i);
    int i1 = neighbor->index (edge.first->vertex (edge.second));
    int i2 = neighbor->index (edge.first->vertex (edge.third));
    
    for (i = 0; i < 3; ++i) if (i != i1 && i != i2) break;
    CGAL_assertion (i < 3);
    Vertex_handle s = neighbor->vertex (i);
    Facet f2 (neighbor, 3);
    if (is_finite (s)) {
      dual_facet (medge,
		  side_of_vertex (p, q, inf, s),
		  side_of_vertex (p, q, s, inf),
		  p, q, s, f2);
    }
    return medge;
  }

  Moebius_edge dual_3D (const Edge &edge) {
    CGAL_assertion (dimension () == 3);

    if (inter_edge (edge) == 0) return Moebius_edge ();

    TRACE ("dual_edge_3D ("<<edge.first->vertex(edge.second)->point()<<", "
	   <<edge.first->vertex(edge.third)->point()<<")...\n");

    Moebius_edge medge = bisector (edge);


    Circulator start (edge), i (start);
    Vertex_handle s = i.source(), t = i.target();
    while (side_of_vertex (i.cell(), s, t) != CGAL::ON_UNBOUNDED_SIDE) {
      TRACE ("  vertex inside, passing...\n");
      ++i;
      CGAL_assertion (i != start);
    }
    start = i;
    do {
      if (is_finite (i.vertex ())) {
	dual_facet (medge, 
		    side_of_vertex (i.cell (), s, t),
		    side_of_vertex (i.next_cell (), s, t),
		    i.source(), i.target(), i.vertex(), i.facet());
      }
      ++i;
    } while (i != start);
    TRACE ("dual_edge_3D done.\n");
    return medge;
  }

  void dual_facet (Moebius_edge &medge,
		   const CGAL::Bounded_side side1,
		   const CGAL::Bounded_side side2,
		   const Vertex_handle p,
		   const Vertex_handle q,
		   const Vertex_handle r,
		   const Facet &facet) {
    //		   const Moebius_edge &other) {
    TRACE ("  dual_facet ("<<p->point()<<", "<<q->point()<<", "<<r->point()<<")...\n");
    switch (real_inter_facet (facet)) {
    case 0: TRACE ("    real 0\n"); break;
    case 1:
      TRACE ("    real 1, ");
      switch (effective_inter_facet (facet)) {
      case 0: TRACE ("eff 0\n"); break;
      case 1:
	{
	  TRACE ("eff 1\n");
	  Bare_point v = vertex (p, q, r);
	  if (orientation (p, q, r) == CGAL::POSITIVE) {
	    TRACE ("      start ("<<v<<")\n");
	    medge.start (v);
	  } else { 
	    medge.stop (v);
	    TRACE ("      stop  ("<<v<<")\n");
	  }
	}
	break;
      default:
	CGAL_assertion (false);
      }
      break;
    case 2:
      TRACE ("    real 2, ");
      switch (effective_inter_facet (facet)) {
      case 0: TRACE ("eff 0\n"); break;
      case 1:
	{
	  TRACE ("eff 1\n");
	  CGAL_assertion (side1 != side2);
	  CGAL_assertion (side1 != CGAL::ON_BOUNDARY);
	  CGAL_assertion (side2 != CGAL::ON_BOUNDARY);
	  Segment v = vertices (p, q, r);

	  if (side1 == CGAL::ON_UNBOUNDED_SIDE) {
	    TRACE ("      start ("<<v.source()<<")\t[stop  "<<v.target()<<"]\n");
	    medge.start (v.source ());
	  } else {
	    TRACE ("      stop  ("<<v.target()<<")\t[start "<<v.source()<<"]\n");
	    medge.stop (v.target ());
	  }
	}
	break;
      case 2:
	{
	  TRACE ("eff 2\n");
	  Segment v = vertices (p, q, r);
	  TRACE ("      start ("<<v.source()<<")\t stop ("<<v.target()<<")\n");
	  medge.start (v.source ());
	  medge.stop (v.target ());
	}
	break;
      default:
	CGAL_assertion (false);
      }
      break;
    default:
      CGAL_assertion (false);
    }
  }

  // -- predicates --
  CGAL::Orientation orientation (const Bare_point &p, 
				 const Bare_point &q, 
				 const Bare_point &r) const
    { return gt().orientation_2_object () (p, q, r); }

  CGAL::Orientation orientation (const Vertex_handle &p, 
				 const Vertex_handle &q, 
				 const Vertex_handle &r) const
    { 
      CGAL_assertion (is_finite (p)); 
      CGAL_assertion (is_finite (q)); 
      CGAL_assertion (is_finite (r));
      return orientation (p->point(), q->point(), r->point());
    }

  CGAL::Orientation moebius_orientation (const Point &p, const Point &q) const
    { return gt().moebius_orientation_2_object () (p, q); }
  CGAL::Orientation moebius_orientation (const Vertex_handle &p, 
					 const Vertex_handle &q) const
    {
      CGAL_assertion (is_finite (p)); CGAL_assertion (is_finite (q));
      return moebius_orientation (p->point(), q->point());
    }

  bool has_circle (const Point &p, const Point &q) const
    { return gt().has_circle_2_object () (p, q); }
  bool has_circle (const Vertex_handle &p, const Vertex_handle &q) const
    { 
      CGAL_assertion (is_finite (p)); CGAL_assertion (is_finite (q));
      return has_circle (p->point(), q->point());
    }

  bool circle_cross_line (const Point &p, const Point &q, const Point &r) const
    { TRACE ("        circle_cross_line ("<<p<<", "<<q<<", "<<r<<") = ");
      bool b = gt().circle_cross_line_2_object () (p, q, r);
      TRACE (b << "\n"); return b;
    }
  bool circle_cross_line (const Vertex_handle &p, 
			  const Vertex_handle &q, 
			  const Vertex_handle &r) const
    { 
      CGAL_assertion (is_finite (p)); 
      CGAL_assertion (is_finite (q)); 
      CGAL_assertion (is_finite (r));
      return circle_cross_line (p->point(), q->point(), r->point());
    }
  
  CGAL::Oriented_side circle_side_of_center (const Point &p, 
					     const Point &q, 
					     const Point &r) const
    { TRACE ("        circle_side_of_center ("<<p<<", "<<q<<", "<<r<<") = ");
      CGAL::Oriented_side s = gt().circle_side_of_center_2_object () (p, q, r);
      TRACE (s << "\n"); return s;
    }
  CGAL::Oriented_side circle_side_of_center (const Vertex_handle &p, 
					     const Vertex_handle &q, 
					     const Vertex_handle &r) const
    { 
      CGAL_assertion (is_finite (p)); 
      CGAL_assertion (is_finite (q)); 
      CGAL_assertion (is_finite (r));
      return circle_side_of_center (p->point(), q->point(), r->point());
    }
  

  CGAL::Bounded_side side_of_vertex (const Point &p,const Point &q, 
				     const Point &r, const Point &s)
    { return gt().side_of_vertex_2_object () (p, q, r, s); }
  
  CGAL::Bounded_side side_of_vertex (const Cell_handle &c, 
				     const Vertex_handle &p, 
				     const Vertex_handle &q) {
    if (c->side_is_set ()) return c->side ();
    int ip = c->index (p), iq = c->index (q);
    CGAL::Bounded_side side = side_of_vertex (p, q, 
					      c->vertex (Triangulation_utils_3::next_around_edge (ip, iq)),
					      c->vertex (Triangulation_utils_3::next_around_edge (iq, ip)));
    c->set_side (side);
    return side;
  }


  CGAL::Bounded_side side_of_vertex (const Vertex_handle &p,
				     const Vertex_handle &q,
				     const Vertex_handle &r,
				     const Vertex_handle &s) {
    CGAL::Bounded_side result = side_of_vertex_ (p, q, r, s);
    TRACE (") = " << result << "\n");
    return result;
  }
  CGAL::Bounded_side side_of_vertex_ (const Vertex_handle &p,
				      const Vertex_handle &q,
				      const Vertex_handle &r,
				      const Vertex_handle &s) { 
    TRACE ("    side_of_vertex ("<<p->point()<<", "<<q->point()<<", ");

    CGAL_assertion (is_finite (p));
    CGAL_assertion (is_finite (q));
    
    if (is_finite (r)) {
      TRACE (r->point()<<", ");
      if (is_finite (s)) {
	TRACE (s->point());
	return side_of_vertex (p->point (), q->point (), 
			       r->point (), s->point ());
      }
      TRACE ("<inf>");
      // r finite, s infinite
      if (moebius_orientation (p, q) == CGAL::ZERO && moebius_orientation (q, r) == CGAL::ZERO) {
	return CGAL::Bounded_side (orientation (p, q, r));
      }
      return CGAL::ON_UNBOUNDED_SIDE;
    }
    TRACE ("<inf>, ");
    CGAL_assertion (is_finite (s));
    TRACE (s->point());
    // r infinite, s finite
    if (moebius_orientation (p, q) == CGAL::ZERO 
	&& moebius_orientation (q, s) == CGAL::ZERO) {
      return CGAL::Bounded_side (orientation (p, s, q));
    }
    return CGAL::ON_UNBOUNDED_SIDE;
  }

  CGAL::Oriented_side circle_side_of_center (const Point &p, 
					     const Point &q, 
					     const Point &r, const Point &s) const
    { return gt().circle_side_of_center_2_object () (p, q, r ,s); }

  CGAL::Oriented_side circle_side_of_center (const Vertex_handle &p, 
					     const Vertex_handle &q, 
					     const Vertex_handle &r, 
					     const Vertex_handle &s) const
    { 
      CGAL_assertion (is_finite (p)); 
      CGAL_assertion (is_finite (q)); 
      CGAL_assertion (is_finite (r)); 
      CGAL_assertion (is_finite (s));
      return circle_side_of_center (p->point(), q->point(), 
				    r->point(), s->point());
    }

  CGAL::Oriented_side circle_side_of_vertex (const Point &p, 
					     const Point &q, 
					     const Point &r, 
					     const Point &s) const
    { CGAL::Oriented_side o = gt().circle_side_of_vertex_2_object () (p, q, r, s);
      TRACE ("            circle_side_of_vertex ("<<p<<", "<<q<<", "<<r<<", "<<s<<") = "<<o<<"\n");
      return o;
    }
  CGAL::Oriented_side circle_side_of_vertex (const Vertex_handle &p, 
					     const Vertex_handle &q, 
					     const Vertex_handle &r, const Vertex_handle &s) const
    { 
      CGAL_assertion (is_finite (p)); 
      CGAL_assertion (is_finite (q)); 
      CGAL_assertion (is_finite (r)); 
      CGAL_assertion (is_finite (s));
      return circle_side_of_vertex (p->point(), q->point(), 
				    r->point(), s->point());
    }

  // -- constructors --
  Circle construct_circle (const Point &p, const Point &q)
    { return gt().construct_circle_2_object () (p, q); }

  Circle construct_circle (const Vertex_handle &p, const Vertex_handle &q)
    { return construct_circle (p->point (), q->point ()); }

  Circle construct_circle (const Edge &edge) 
    { return construct_circle (edge.first->vertex (edge.second), 
			       edge.first->vertex (edge.third)); }

  Line construct_line (const Point &p, const Point &q)
    { return gt().construct_line_2_object () (p, q); }

  Line construct_line (const Vertex_handle &p, const Vertex_handle &q)
    { return construct_line (p->point (), q->point ()); }

  Line construct_line (const Edge &edge) 
    { return construct_line (edge.first->vertex (edge.second), 
			     edge.first->vertex (edge.third)); }


  Bare_point vertex (const Circulator &i) {
    return vertex (i.source(), i.target(), i.vertex ());
  }
  Bare_point vertex (const Vertex_handle &p, 
		     const Vertex_handle &q, 
		     const Vertex_handle &r) {
    CGAL_assertion (is_finite (p)); 
    CGAL_assertion (is_finite (q)); 
    CGAL_assertion (is_finite (r));
    return vertex (p->point(), q->point(), r->point());
  }

  Bare_point vertex (const Point &p, const Point &q, const Point &r) {
    CGAL_assertion (moebius_orientation (p, q) == CGAL::ZERO);
    CGAL_assertion (moebius_orientation (q, r) == CGAL::ZERO);
    return gt().construct_vertex_2_object () (p, q, r);
  }


  Segment vertices (const Circulator &i) {
    return vertices (i.source(), i.target(), i.vertex ());
  }
  Segment vertices (const Vertex_handle &p, 
		    const Vertex_handle &q, 
		    const Vertex_handle &r) {
    CGAL_assertion (is_finite (p)); 
    CGAL_assertion (is_finite (q)); 
    CGAL_assertion (is_finite (r));
    return vertices (p->point(), q->point(), r->point());
  }
  Segment vertices (const Point &p, const Point &q, const Point &r) {
    if (moebius_orientation (p, q) == CGAL::ZERO) {
      CGAL_assertion (moebius_orientation (q, r) != CGAL::ZERO);
      return gt().construct_line_vertices_2_object () (p, q, r);
    }
    return gt().construct_circle_vertices_2_object () (p, q, r);
  }

  // -- utils --

  bool real_inter_facet_is_set (const Facet &f) const {
    return inter_facet (f) & 0x08;
  }
  bool effective_inter_facet_is_set (const Facet &f) const {
    return inter_facet (f) & 0x80;
  }

  void set_real_inter_facet (const Facet &f, int n) {
    TRACE ("        set_real_inter_facet ("<<n<<")...\n");
    int i = inter_facet (f);
    int j = (i & 0xF0) | ((n & 0x07) | 0x08);
    set_inter_facet (f, j);
  }
  void set_effective_inter_facet (const Facet &f, int n) {
    TRACE ("        set_effective_inter_facet ("<<n<<")...\n");
    int i = inter_facet (f);
    int j = (((n << 4) & 0x70) | 0x80) | (i & 0x0F);
    set_inter_facet (f, j);
  }
  void set_inter_facet (const Facet &f, int n) {
    CGAL_assertion (dimension () > 1);
    f.first->set_inter_facet (f.second, n);
    if (dimension () > 2)
      f.first->neighbor(f.second)->set_inter_facet (f.first->mirror_index(f.second), n);
  }

  int inter_facet (const Facet &f) {
    return f.first->inter_facet (f.second);
  }
  int real_inter_facet (const Facet &f) {
    return inter_facet (f) & 0x07;
  }
  int effective_inter_facet (const Facet &f) {
    return (inter_facet (f) & 0x70) >> 4;
  }


  void set_inter_edge (const Edge &e, int inter) {
    Vertex_handle e1 (e.first->vertex(e.second));
    Vertex_handle e2 (e.first->vertex(e.third));
    TRACE ("  set_inter_edge " << inter << " "
	   << "(" << e1->point() << ") "
	   << "(" << e2->point() << ")\n");
    switch (dimension ()) {
    case 1:
      e.first->set_inter_edge (e.second, e.third, inter);
      break;
    case 2:
      {
	e.first->set_inter_edge (e.second, e.third, inter);
	int i;
	for (i = 0; i < 3; ++i) if (i != e.second && i != e.third) break;
	CGAL_assertion (i < 3);
	Cell_handle neighbor = e.first->neighbor (i);
	neighbor->set_inter_edge (e1, e2, inter);
      } 
      break;
    case 3:
      {
	Cell_circulator start (_rt.incident_cells (e));
	Cell_circulator c = start;
	do {
	  c->set_inter_edge (e1, e2, inter);
	} while (++c != start);
      }
      break;
    default:
      std::cerr << "Moebius_diagram_2::dimension () == " << dimension () << "\n";
      CGAL_assertion (false);
    }
  }

};
#ifdef CGAL_USE_QT

template <class Gt>
Qt_widget& operator<<(Qt_widget& widget, const Moebius_diagram_2<Gt>& d)
{
  typedef Moebius_diagram_2<Gt> Diagram;
  typedef typename Diagram::Moebius_edge_const_iterator Iterator;

  Iterator i = d.moebius_edges_const_begin (), end = d.moebius_edges_const_end ();
  while (i != end) {
    widget << *i;
    ++i;
  }
  return widget;
}
#endif

CGAL_END_NAMESPACE

#endif// CHR_MOEBIUS_DIAGRAM_2_H
