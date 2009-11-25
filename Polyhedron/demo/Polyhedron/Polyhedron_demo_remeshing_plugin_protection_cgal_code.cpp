#include <CGAL/config.h>
#include "Polyhedron_type.h"
#include "C2t3_type.h"
#include <set>
#include <map>

typedef Tr::Geom_traits::FT FT;

struct Less {
  template <typename Handle>
  bool operator()(const Handle& va, const Handle& vb) const {
    return &*va < &*vb;
  }
};

struct Insert_spheres {
  typedef Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef Polyhedron::size_type size_type;

  typedef Tr::Geom_traits::Point_3 Point_3;
  typedef Tr::Vertex_handle Tr_vertex_handle;

  typedef std::set<Vertex_const_handle, Less> Vertices_set;
  typedef std::map<Vertex_const_handle, size_type, Less> Vertices_counter;

  typedef std::set<Halfedge_const_handle, Less> Border_edges_set;

  Border_edges_set edges_to_consider;
  Vertices_counter border_vertices;
  Vertices_set corner_vertices;

  typedef std::list<Tr::Geom_traits::Point_3> Polyline;
  
  typedef std::vector<Polyline> Polylines;

  typedef std::vector<bool> Polyline_is_a_cycle;

  typedef unsigned Polyline_id;

  typedef std::set<Point_3> Hidden_balls;
  typedef std::multimap<Tr_vertex_handle, Polyline_id> Tr_corner_vertices;
  
  Polylines polylines;
  Polyline_is_a_cycle polyline_is_a_cycle;
  Tr_corner_vertices tr_corner_vertices;
  Hidden_balls hidden_balls;

  C2t3& c2t3;
  Polyhedron* pMesh;
  const FT size;

  Halfedge_const_handle canonical(Halfedge_const_handle he) 
  {
    const Halfedge_const_handle& op = he->opposite();
    if(Less()(he, op))
      return he;
    else 
      return op;
  }

  /** Follow a polyline or a polygon, from the halfedge he. */
  void follow_half_edge(const Halfedge_const_handle he, const bool is_cycle)
  {
    Border_edges_set::iterator it = edges_to_consider.find(canonical(he));
    if(it == edges_to_consider.end()) {
      return;
    }

    Polyline polyline;
    polyline.push_back(he->opposite()->vertex()->point());
    // std::cerr << "Start: " << he->opposite()->vertex()->point() << std::endl;

    Halfedge_const_handle current_he = he;
    do {
      CGAL_assertion(current_he->is_border() || 
                     current_he->opposite()->is_border());
      CGAL_assertion_code(const size_type n = )edges_to_consider.erase(canonical(current_he));
      CGAL_assertion(n > 0);
      Vertex_const_handle v = current_he->vertex();
      polyline.push_back(v->point());
      // std::cerr << v->point() << std::endl;
      if(corner_vertices.count(v) > 0) break;
      Polyhedron::Halfedge_around_vertex_const_circulator 
        loop_he = v->vertex_begin(), end(loop_he);
      ++loop_he;
      // CGAL_assertion((&*loop_he) != (&*current_he) );
      while((&*loop_he) == (&*current_he) ||
            (!loop_he->is_border() && !loop_he->opposite()->is_border()) ) {
        ++loop_he;
        // CGAL_assertion((&*loop_he) != (&*current_he) );
      }
      current_he = loop_he->opposite();
    } while(current_he != he );

    // if(current_he == he)
    //   std::cerr << "New polyline, of size " << polyline.size() << std::endl;
    // else
    //   std::cerr << "New polygon (cycle), of size " << polyline.size() << std::endl;

    polyline_is_a_cycle.push_back(is_cycle);
    polylines.push_back(polyline);
  }

  /** Loop around a corner vertex, and try to follow a polyline of border
      edges, from each incident edge. */
  void loop_around_corner(const Vertex_const_handle v)
  {
    Polyhedron::Halfedge_around_vertex_const_circulator 
      he = v->vertex_begin(), end(he);
    do {
      CGAL_assertion(he->vertex() == v);
      follow_half_edge(he->opposite(), false);
      ++he;
    } while(he != end);
  }

  /** For a non-corner vertex v (that is incident to two border edges),
      mesure the angle between the two edges, and mark the vertex as corner
      edge, if the angle is < 120°. **/
  void mesure_angle(const Vertex_const_handle v)
  {
    Halfedge_const_handle e1;
    Halfedge_const_handle e2;
    Polyhedron::Halfedge_around_vertex_const_circulator he = v->vertex_begin(), end(he);
    // std::cerr << "mesure_handle(" << (void*)(&*v)
    //           << " = " << v->point() << ")";
    bool first = true;
    do {
      CGAL_assertion(he->vertex() == v);
      // std::cerr << he->opposite()->vertex()->point() << std::endl;
      if(he->is_border() || he->opposite()->is_border()) {
        if(first) {
          e1 = he;
          first = false;
        }
        else {
          CGAL_assertion(e2 == Halfedge_const_handle());
          e2 = he;
        }
        // std::cerr << "x";
      }
      // else
        // std::cerr << ".";
      ++he;
    } while(he != end);
      // std::cerr << "\n";
    const Point_3 pv = v->point();
    const Point_3 pa = e1->opposite()->vertex()->point();
    const Point_3 pb = e2->opposite()->vertex()->point();
    const Tr::Geom_traits::Vector_3 av = pv - pa;
    const Tr::Geom_traits::Vector_3 bv = pv - pb;
    const FT sc_prod = av * bv;
    if( sc_prod >= 0 ||
        (sc_prod < 0 && 
         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
    {
      // std::cerr << "Corner (" << pa << ", " << pv
      //           << ", " << pb << ")\n";
      corner_vertices.insert(v);
    }
  }

  const Tr_vertex_handle insert(const Point_3& p, const Polyline_id polyline_id)
  {
    const Tr::Vertex_handle v = c2t3.triangulation().insert(p);
    if(v != Tr_vertex_handle()) {
      v->info() = polyline_id;
    }
    else {
      hidden_balls.insert(p);
    }
    return v;
  }

  void protect(Polyline::const_iterator begin,
               Polyline::const_iterator end,
               const Polyline_id polyline_id) {
    // for(Polyline::const_iterator it
    Polyline::const_iterator end2 = end;
    --end2;
    CGAL_assertion(begin != end2);

    FT distance = 0;
    for(Polyline::const_iterator it = begin;
        it != end2;) {
      const Point& a = *it;
      const Point& b = *++it;
      distance += CGAL_NTS sqrt( CGAL::squared_distance(a, b) );
    }
    // std::cerr << "Distance:   " << distance << std::endl;
    // std::cerr << "Size:       " << size << std::endl;
    const size_type n = static_cast<size_type>(std::ceil(distance / size) + 0.5);
    const FT local_size = distance / n;
    // std::cerr << n << std::endl;
    // std::cerr << "Local size: " << local_size << std::endl;
    // CGAL_assertion(local_size < size);
    const FT r2 = CGAL::square(local_size/1.5);
    Point_3 a(begin->point(), r2);
    Point_3 b(end2->point(), r2);
    const Tr_vertex_handle va = insert(a, polyline_id);
    FT small_distance_to_go = local_size;
    Polyline::const_iterator it = begin;
    while(it != end2) {
      const Point& a = *it;
      const Point& b = *++it;
      // std::cerr << "segment( " << a << ", " << b << ")\n";
      // std::cerr << "small_distance_to_go=" << small_distance_to_go << std::endl;
      const FT d = CGAL_NTS sqrt(squared_distance(a, b));
      // std::cerr << "d=" << d << std::endl;
      FT pos = small_distance_to_go;
      if(pos < d) {
        for(; pos < d;
            pos += local_size)
        {
          const Point p = a +
            pos * ( b - a ) / d;
          insert(Point_3(p, r2), polyline_id);
          // std::cerr << ".";
        }
        // pos -= local_size;
        small_distance_to_go = pos - d;
      }
      else  {
        small_distance_to_go -= d;
      }
      // std::cerr << "\n";
    }
    const Tr_vertex_handle vb = insert(b, polyline_id);

    if(!polyline_is_a_cycle[polyline_id])
    {
      // CGAL_assertion(va == Tr_vertex_handle() || 
      //                vb == Tr_vertex_handle() || 
      //                va != vb);
      if(va != Tr_vertex_handle())
        tr_corner_vertices.insert(std::make_pair(va, polyline_id));
      if(vb != Tr_vertex_handle())
        tr_corner_vertices.insert(std::make_pair(va, polyline_id));
    }

    // std::cerr << "One polyline is protected!\n";
  }

  void separate_balls() {
    Tr& tr = c2t3.triangulation();
    bool restart = false;
    do {
      restart = false;
      for(Tr::Finite_edges_iterator eit = tr.finite_edges_begin(), 
            end = tr.finite_edges_end(); eit != end; ++eit) 
      {
        const Tr_vertex_handle& va = eit->first->vertex(eit->second);
        const Tr_vertex_handle& vb = eit->first->vertex(eit->third);
        if(non_adjacent_but_intersect(va, vb))
        {
          std::cerr << "Balls " << va->point() << " (on polyline "
                    << va->info() << ")"
                    << " and " << vb->point() << " (on polyline "
                    << vb->info() << ")"
                    << " intersect.\n";

          restart = true;
          if(va->point().weight() > vb->point().weight())
            refine_ball(va);
          else
            refine_ball(vb);
          break;
        }
      }
    }
    while(restart);
  }

  bool non_adjacent_but_intersect(const Tr_vertex_handle& va,
                                  const Tr_vertex_handle& vb) const
  {
    typedef Tr_corner_vertices::const_iterator const_iterator;
    typedef std::pair<const_iterator, const_iterator> Range;
    typedef std::set<Polyline_id> P_id_set;
    Range range_a = tr_corner_vertices.equal_range(va);
    Range range_b = tr_corner_vertices.equal_range(vb);

    bool non_adjacent;

    if(range_a.first == range_a.second &&
       range_b.first == range_b.second) 
    {
      non_adjacent = ( va->info() != vb->info() );
    }
    else {
      P_id_set s_a;
      P_id_set s_b;
      for(const_iterator it = range_a.first; it != range_a.second; ++it)
      {
        s_a.insert(it->second);
      }
      s_a.insert(va->info());
      for(const_iterator it = range_b.first; it != range_b.second; ++it)
      {
        s_b.insert(it->second);
      }
      s_b.insert(va->info());
      P_id_set intersection;
      std::set_intersection(s_a.begin(), s_a.end(),
                            s_b.begin(), s_b.end(),
                            std::inserter(intersection, intersection.begin()));
      non_adjacent = intersection.empty();
    }
    if(non_adjacent)
    {
      const Point_3& a = va->point();
      const Point_3& b = vb->point();
      if( CGAL_NTS sqrt( CGAL::squared_distance(a, b) ) <
          CGAL_NTS sqrt( a.weight() ) + CGAL_NTS sqrt( b.weight() ) )
      {
        return true;
      }
    }
    return false;
  }

  void refine_ball(Tr_vertex_handle v) {
    typedef Tr_corner_vertices::iterator iterator;
    typedef std::pair<iterator, iterator> Range;
    Range r = tr_corner_vertices.equal_range(v);

    if(r.first != r.second)
    {
      std::cerr << "Refine vertex ball " << v->point() << "\n";
      Point_3 new_point(v->point().point(), v->point().weight() / 4);
      Polyline_id info = v->info();
      c2t3.triangulation().remove(v);
      Tr_vertex_handle new_v = c2t3.triangulation().insert(new_point);
      new_v->info() = info;
      iterator hint = tr_corner_vertices.begin();
      for(iterator it = r.first; it != r.second; ++it) {
        tr_corner_vertices.insert(hint, std::make_pair(new_v, it->second));
      }
      tr_corner_vertices.erase(r.first, r.second);
    }
    else {
      std::cerr << "Refine ball " << v->point() << "\n";
      Point_3 new_point(v->point().point(), v->point().weight() / 16);
      Polyline_id info = v->info();
      c2t3.triangulation().remove(v);
      Tr_vertex_handle new_v = c2t3.triangulation().insert(new_point);
      if(new_v != Tr_vertex_handle()) 
      {
        std::cerr << "Refined ball " << new_point << " is now hidden!\n";
        new_v->info() = info;
      }
    }
  }

  Insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size)
    : c2t3(c2t3), pMesh(pMesh), size(size)
  {
    // That call orders the set of edges of the polyhedron, so that the
    // border edges are at the end of the sequence of edges.
    pMesh->normalize_border();

    // Iterate over border edges, and find out which vertices are corner
    // vertices (more than two incident border edges).
    for(Polyhedron::Edge_const_iterator 
          eit = pMesh->border_edges_begin (),
          end = pMesh->edges_end();
        eit != end; ++eit)
    {
      edges_to_consider.insert(canonical(eit));
      Polyhedron::Vertex_const_handle v = eit->vertex();
      for(Polyline_id i = 0; i < 2; ++i) {
        if(++border_vertices[v] == 3)
          corner_vertices.insert(v);
        v = eit->opposite()->vertex();
      }

    }
    std::cerr << "Corner vertices: " << corner_vertices.size() << std::endl;
    std::cerr << "Border vertices: " << border_vertices.size() << std::endl;
    
    // Iterate over non-corner border vertices, and mesure the angle.
    for(Vertices_counter::iterator it = border_vertices.begin(),
        end = border_vertices.end(); it != end; ++it)
    {
      const Vertex_const_handle v = it->first;
      if(corner_vertices.count(v) == 0) {
        CGAL_assertion(it->second == 2);
        mesure_angle(v);
      }
    }
    std::cerr << "New corner vertices: " << corner_vertices.size() << std::endl;

    // Follow the polylines...
    for(Vertices_set::iterator it = corner_vertices.begin(),
        end = corner_vertices.end(); it != end; ++it)
    {
      loop_around_corner(*it);
    }

    // ... and the cycles.
    while(! edges_to_consider.empty() ) {
      follow_half_edge(*edges_to_consider.begin(), true);
    }

    // Then insert the protecting balls.
    for(Polyline_id i = 0; i < polylines.size(); ++i)
    {
      protect(polylines[i].begin(), polylines[i].end(), i);
    }

    std::cerr << "Number of polylines: " << polylines.size() << "\n";
    std::cerr << "Number of protecting balls: " 
              << c2t3.triangulation().number_of_vertices() << std::endl;
    std::cerr << "Number of hidden balls: " << hidden_balls.size() << "\n";
    for(Hidden_balls::const_iterator it = hidden_balls.begin(),
          end = hidden_balls.end(); it != end; ++it)
    {
      std::cerr << "Hidden ball " << *it << ", hidden by the ball "
                << c2t3.triangulation().nearest_power_vertex(*it)->point()
                << "\n";
    }

    separate_balls();
  }
};

void insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size) {

  Insert_spheres go(c2t3, pMesh, size);
}
