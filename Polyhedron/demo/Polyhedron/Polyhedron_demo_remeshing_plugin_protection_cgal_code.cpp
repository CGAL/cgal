#include <CGAL/config.h>
#include "Polyhedron_type.h"
#include "C2t3_type.h"
#include <set>
#include <map>
#include <CGAL/tuple.h>

typedef Tr::Geom_traits::FT FT;

struct Less {
  template <typename Handle>
  bool operator()(const Handle& va, const Handle& vb) const {
    return &*va < &*vb;
  }
};

typedef std::vector<Tr::Geom_traits::Point_3> Polyline;
typedef unsigned Polyline_id;

struct Ball_context {
  typedef Tr::Vertex_handle Tr_vertex_handle;

  Tr_vertex_handle pred; // predecessor of the current non-corner ball on a
                         // polyline
  Tr_vertex_handle succ; // successor..
  Polyline_id id; // polyline id
  Polyline::const_iterator iterator; // in the polyline, iterator to the
                                     // point that precedes the ball
                                     // (*it, *(it+1)) is the segment.
  FT position; // distance between the center of the current ball, and
               // halfedge->opposite()->vertex()->point()
  bool is_corner;
};

Vertex_info::Vertex_info() : context_(0) {
}

Vertex_info::~Vertex_info() {
  if(context_)
  {
    delete context_;
    context_ = 0;
  }
}

Ball_context* Vertex_info::context() {
  if(!context_) {
    context_ = new Ball_context;
  }
  return context_;
}

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

  typedef std::vector<Polyline> Polylines;

  typedef std::vector<bool> Polyline_is_a_cycle;
  typedef std::set<Point_3> Hidden_balls;
  typedef CGAL::cpp0x::tuple<Polyline_id, 
                             Polyline::const_iterator,
                             Tr_vertex_handle> Corner_context;
  typedef std::multimap<Tr_vertex_handle, 
                        Corner_context> Tr_corner_vertices;
  
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
    Tr::Vertex_handle v = c2t3.triangulation().insert(p);
    if(v == Tr_vertex_handle()) {
      v = c2t3.triangulation().nearest_power_vertex(p);
      hidden_balls.insert(p);
      std::cerr << "Hidden ball " << p << " (id="
                << polyline_id << "), hidden by the ball "
                << v->point()
                << " (" 
                << print_context(v)
                << ")\n";
    }
    v->info().context()->id = polyline_id;
    return v;
  }

  void cover(Polyline::const_iterator begin,
             Polyline::const_iterator end2,
             const Polyline_id polyline_id,
             FT alpha, // alpha is the desired radius of covering ball
                       // the real radius must be lower.
             Tr_vertex_handle va,
             Tr_vertex_handle vb,
             const FT position_begin,
             const FT position_end)
  {
    std::cerr << "cover(" << va->point() << ", " 
              << vb->point() << ", " << polyline_id << ")\n";
    FT radius_begin = CGAL_NTS sqrt(va->point().weight());
    FT radius_end = CGAL_NTS sqrt(vb->point().weight());

    if(begin == end2 && position_begin >= position_end)
      return;

    FT distance = 0;
    for(Polyline::const_iterator it = begin;
        it != end2; /*'it' is incremented in the body*/ ) {
      const Point& a = *it;
      const Point& b = *++it;
      distance += CGAL_NTS sqrt( CGAL::squared_distance(a, b) );
    }
    distance -= radius_begin;
    distance -= radius_end;
    distance -= position_begin;
    distance += position_end;
    if(distance < 0) return;
    // distance -= (alpha * 2) / 3;
    // std::cerr << "Distance:   " << distance << std::endl;
    // std::cerr << "Size:       " << size << std::endl;

    // The desired radius is alpha. The desired inter-distance (between
    // balls centers) is 4*alpha/3
    // const FT desired_inter_distance = 4 * alpha / 3;

    // Let r be the real radius (r <= alpha).
    // Let inter_distance be the real inter-distance:
    //   inter_distance = 4 * r / 3 <= desired_inter_distance
    // Let n be the number of steps: the number of inserted balls is n+1.
    // Then : 
    //   distance = r / 3 + n * (4 * r /3 ) + r / 3
    //
    //   3* distance = (4 * n + 2 ) * r
    //
    //       1     3 * distance      
    //   n = - * ( ------------ - 2 )
    //       4           r
    //
    // n >= 1 then r <= distance / 2
    const int n = (std::max)
      (1, 
       static_cast<int>(std::ceil(((3 * distance / alpha) - 2) / 4) + 0.5));
    
    const FT r = (3 * distance ) / ( 4 * n + 2) ;
    // std::cerr << n << std::endl;
    // std::cerr << "Local size: " << local_size << std::endl;
    // CGAL_assertion(local_size < size);
    const FT inter_distance = 4 * r / 3;
    const FT r2 = CGAL::square(r);

    std::cerr << "distance = " << distance
              << "\nalpha = " << alpha
              << "\nn = " << n
              << "\nr = " << r
              << std::endl;

    Point_3 a(begin->point(), r2);
    Point_3 b(end2->point(), r2);
    FT small_distance_to_go = position_begin + radius_begin + (r / 3);
    Polyline::const_iterator it = begin;
    bool first = true;
    Tr_vertex_handle first_non_corner;
    Tr_vertex_handle last_non_corner;
    Tr_vertex_handle previous = va;
    int counter = 0;
    while(it <= end2 && counter <= n)//small_distance_to_go < (position_end - radius_end))) 
    {
      const Point& a = *it;
      const Point& b = *(it+1);
      std::cerr << "segment( " << a << ", " << b << ")\n";
      std::cerr << "small_distance_to_go=" << small_distance_to_go << std::endl;
      const FT d = CGAL_NTS sqrt(squared_distance(a, b));
      std::cerr << "d=" << d << std::endl;
      FT pos = small_distance_to_go;
      if(pos < d && counter <= n) {
        for(; pos < d && counter <= n;
            pos += inter_distance)
        {
          std::cerr << "pos=" << pos << "" << std::endl;
          const Point p = a +
            pos * ( b - a ) / d;
          ++counter;
          Tr_vertex_handle current = insert(Point_3(p, r2), polyline_id);
          if(current != Tr_vertex_handle()) {
            current->info().context()->pred = previous;
            current->info().context()->iterator = it;
            current->info().context()->position = pos;
            current->info().context()->is_corner = false;
            previous->info().context()->succ = current;

            // at the end of the loop, last_non_corner is the last non-corner
            // ball on the polyline
            last_non_corner = current;
            if(first) {
              // and first_non_corner is the first non-corner ball on the
              // polyline
              first_non_corner = current;
              first = false;
            }
            previous = current;
          }
        }
        // pos -= inter_distance;
        small_distance_to_go = pos - d;
      }
      else  {
        small_distance_to_go -= d;
      }
      ++it;
      // std::cerr << "\n";
    }
    if(!polyline_is_a_cycle[polyline_id])
    {
      // CGAL_assertion(va == Tr_vertex_handle() || 
      //                vb == Tr_vertex_handle() || 
      //                va != vb);
      // if(va != Tr_vertex_handle())
      if(va->info().context()->is_corner) {
        std::cerr << "Corner vertex va " << va->point()
                  << ", neighbor " << *first_non_corner << std::endl;
        tr_corner_vertices.insert(std::make_pair(va, 
                                                 Corner_context(polyline_id, begin, first_non_corner)));
        CGAL_assertion(first_non_corner->info().context()->pred == va);
      }
      if(vb->info().context()->is_corner) {
        std::cerr << "Corner vertex vb " << vb->point()
                  << ", neighbor " << *last_non_corner << std::endl;
        tr_corner_vertices.insert(std::make_pair(vb, 
                                                 Corner_context(polyline_id, end2, last_non_corner)));
      }
    }
    else 
    {
      vb = va;
    }
    if(last_non_corner != Tr_vertex_handle()) {
      last_non_corner->info().context()->succ = vb;
    }
    // if(vb != Tr_vertex_handle())
    vb->info().context()->pred = last_non_corner;

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
          std::cerr << "Balls " << va->point() << " ("
                    << print_context(va)
                    << ") and " << vb->point() << " ("
                    << print_context(vb) 
                    << ") intersect.\n";

          if(// !va->info().context()->is_corner &&
             va->point().weight() > vb->point().weight()) {
            restart = true;
            refine_ball(va);
            break;
          }
          else // if(!vb->info().context()->is_corner)
          {
            restart = true;
            refine_ball(vb);
            break;
          }
        }
      }
    }
    while(restart);
  }

  std::string print_context(Tr_vertex_handle va) {
    std::stringstream s;
    using CGAL::cpp0x::get;
    if(!va->info().context()->is_corner)
      s << "on polyline " << va->info().context()->id;
    else {
      typedef Tr_corner_vertices::const_iterator const_iterator;
      typedef std::pair<const_iterator, const_iterator> Range;
      Range range_a = tr_corner_vertices.equal_range(va);
      s << "on polylines";
      for(const_iterator it = range_a.first; it != range_a.second; ++it)
        s << " " << get<0>(it->second);
    }
    return s.str();
  }

  bool non_adjacent_but_intersect(const Tr_vertex_handle& va,
                                  const Tr_vertex_handle& vb) const
  {
    using CGAL::cpp0x::get;

    typedef Tr_corner_vertices::const_iterator const_iterator;
    typedef std::pair<const_iterator, const_iterator> Range;
    typedef std::set<Polyline_id> P_id_set;
    bool non_adjacent;

    if(!va->info().context()->is_corner && //range_a.first == range_a.second &&
       !vb->info().context()->is_corner)    // range_b.first == range_b.second) 
    {
      non_adjacent = ( va->info().context()->id != vb->info().context()->id );
    }
    else {
      Range range_a = tr_corner_vertices.equal_range(va);
      Range range_b = tr_corner_vertices.equal_range(vb);

      P_id_set s_a;
      P_id_set s_b;
      for(const_iterator it = range_a.first; it != range_a.second; ++it)
      {
        s_a.insert(get<0>(it->second));
      }
      s_a.insert(va->info().context()->id);
      for(const_iterator it = range_b.first; it != range_b.second; ++it)
      {
        s_b.insert(get<0>(it->second));
      }
      s_b.insert(va->info().context()->id);
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
    typedef Tr_corner_vertices::const_iterator const_iterator;
    typedef Tr_corner_vertices::iterator iterator;
    typedef std::pair<iterator, iterator> Range;
    using CGAL::cpp0x::get;
    if(v->info().context()->is_corner) //r.first != r.second)
    {
      Range r = tr_corner_vertices.equal_range(v);
      std::cerr << "Refine vertex ball " << v->point() << "\n"
                << "Number of adjacent balls: " << distance(r.first, r.second)
                << std::endl;
      Point_3 new_point(v->point().point(), v->point().weight() / 4);
      Polyline_id id = v->info().context()->id;
      c2t3.triangulation().remove(v);
      Tr_vertex_handle new_v = insert(new_point, id);
      new_v->info().context()->is_corner = true;
      for(iterator it = r.first; it != r.second; ++it) {
        Tr_vertex_handle neighbor = get<2>(it->second);
        Ball_context* context = neighbor->info().context();
        std::cerr << "  adjacent ball: " << neighbor->point() 
                  << ", id=" << context->id << std::endl;
        if(v == context->pred)
        {
          context->pred = new_v;
        }
        else {
          CGAL_assertion(v == context->succ);
          context->succ = new_v;
        }
      }

      iterator hint = tr_corner_vertices.begin();
      for(iterator it = r.first; it != r.second; ++it) {
        hint = tr_corner_vertices.insert(hint, 
                                         std::make_pair(new_v, it->second));
      }
      for(iterator it = r.first; it != r.second; ++it) {
        refine_ball(get<2>(it->second));
      }
      tr_corner_vertices.erase(r.first, r.second);
    }
    else {
      std::cerr << "Refine ball " << v->point();
      Point_3 new_point(v->point().point(), v->point().weight() / 16);
      Ball_context* context = v->info().context();
      const Polyline_id id = context->id;
      std::cerr << ", id= " << id << "\n";
      const Tr_vertex_handle pred = context->pred;
      const Tr_vertex_handle succ = context->succ;
      const Polyline::const_iterator it = context->iterator;
      // const FT position = context->position;
      
      c2t3.triangulation().remove(v);
      // Tr_vertex_handle new_v = insert(new_point, id);
      // Ball_context* new_context = new_v->info().context();
      // new_context->is_corner = false;
      // new_context->pred = pred;
      // new_context->succ = succ;
      // new_context->iterator = it;
      // new_context->position = position;

      Polyline::const_iterator pred_it, succ_it;
      FT pred_pos, succ_pos;
      boost::tie(pred_it, pred_pos) = get_it_and_pos(pred, id);
      boost::tie(succ_it, succ_pos) = get_it_and_pos(succ, id);
      if(pred_it != Polyline::const_iterator() && 
         succ_it != Polyline::const_iterator()) 
      {
        std::cerr << "OK!\n";
        if(pred_it > succ_it) {
          std::swap(pred_it, succ_it);
          std::swap(pred_pos, succ_pos);
        }
        cover(pred_it, succ_it,
              id,
              CGAL_NTS sqrt(new_point.weight()),
              pred, succ,
              pred_pos, succ_pos);
      }
      else {
        std::cerr << "ERROR!\n";
      }
    }
  }

  std::pair<Polyline::const_iterator, FT>
  get_it_and_pos(Tr_vertex_handle v, Polyline_id id)
  {
    Ball_context* context = v->info().context();
    if(!context->is_corner) {
      if(context->id == id)
        return std::make_pair(context->iterator, context->position);
    }
    else {
      typedef Tr_corner_vertices::const_iterator const_iterator;
      typedef std::pair<const_iterator, const_iterator> Range;
      Range range_a = tr_corner_vertices.equal_range(v);
      for(const_iterator it = range_a.first; it != range_a.second; ++it)
      {
        using CGAL::cpp0x::get;
        if(get<0>(it->second) == id)
          return std::make_pair(get<1>(it->second), 0);
      }
    }
    return std::make_pair(Polyline::const_iterator(), 0);
  }

  Insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size, 
                 bool refine_balls)
    : c2t3(c2t3), pMesh(pMesh), size(size)
  {
    Tr& tr = c2t3.triangulation();

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

    typedef std::vector<std::pair<Tr_vertex_handle,
      Tr_vertex_handle> > Polylines_endpoints;
    Polylines_endpoints polylines_endpoints;

    // Then insert the protecting balls.
    for(Polyline_id i = 0; i < polylines.size(); ++i)
    {
      Tr_vertex_handle va = insert(Point_3(*polylines[i].begin(), 0), i);
      Tr_vertex_handle vb = insert(Point_3(*(--polylines[i].end()), 0), i);
      polylines_endpoints.push_back(std::make_pair(va, vb));
      CGAL_assertion(!polyline_is_a_cycle[i] ||
                     va == vb);
    }

    typedef std::vector<std::pair<FT, FT> > Polylines_endpoints_sq_radii;
    std::map<Tr_vertex_handle, FT> sq_radii;
    Polylines_endpoints_sq_radii polylines_endpoints_sq_radii;
    
    FT size2 = CGAL::square(size);
    for(Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin(), 
          end = tr.finite_vertices_end(); vit != end; ++vit)
    {
      const Point_3& a = vit->point();
      std::vector<Tr_vertex_handle> incident_vertices;
      incident_vertices.reserve(32);
      tr.incident_vertices(vit, std::back_inserter(incident_vertices));
      std::vector<Tr_vertex_handle>::const_iterator vit2 = incident_vertices.begin();
      FT min_sq_distance = CGAL::squared_distance(a, (*vit2)->point());
      for(std::vector<Tr_vertex_handle>::const_iterator end = incident_vertices.end(); 
          vit2 != end; ++vit2)
      {
        FT sq_distance = CGAL::squared_distance(a, (*vit2)->point());
        if(sq_distance < min_sq_distance) {
          min_sq_distance = sq_distance;
        }
      }
      sq_radii[vit] = (std::min)(min_sq_distance, size2);
    }

    for(Polyline_id i = 0; i < polylines.size(); ++i) {
      polylines_endpoints_sq_radii.push_back(std::make_pair(sq_radii[polylines_endpoints[i].first],
                                                            sq_radii[polylines_endpoints[i].second]));
    }

    polylines_endpoints.clear();
    tr.clear();

    // Then insert the protecting balls.
    for(Polyline_id i = 0; i < polylines.size(); ++i)
    {
      Tr_vertex_handle va = insert(Point_3(*polylines[i].begin(),
                                           polylines_endpoints_sq_radii[i].first), i);
      Tr_vertex_handle vb = insert(Point_3(*(polylines[i].end()-1),
                                           polylines_endpoints_sq_radii[i].second), i);
      if(!polyline_is_a_cycle[i]) {
        va->info().context()->is_corner = true;
        vb->info().context()->is_corner = true;
      }
      else {
        CGAL_assertion(va == vb);
      }

      cover(polylines[i].begin(),
            polylines[i].end()-1,
            i,
            CGAL_NTS sqrt((std::min)(va->point().weight(),
                                     vb->point().weight())),
            va,
            vb,
            0,
            0);
    }
    polylines_endpoints_sq_radii.clear();

    unsigned counter = 0;
    for(unsigned i = 0; i < polyline_is_a_cycle.size(); ++i)
      if(polyline_is_a_cycle[i]) ++counter;
    std::cerr << "Number of cycles: " << counter << std::endl;
    std::cerr << "Number of polylines: " << polylines.size() - counter << "\n";
    std::cerr << "Number of protecting balls: " 
              << c2t3.triangulation().number_of_vertices() << std::endl;

    counter = 0;
    for(Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin(), 
          end = tr.finite_vertices_end(); vit != end; ++vit) {
      if(vit->info().context()->is_corner)
        ++counter;
    }
    std::cerr << "Number of vertex balls: " << counter << std::endl;

    std::cerr << "Number of hidden balls: " << hidden_balls.size() << "\n";
    for(Hidden_balls::const_iterator it = hidden_balls.begin(),
          end = hidden_balls.end(); it != end; ++it)
    {
      std::cerr << "Hidden ball " << *it << ", hidden by the ball "
                << c2t3.triangulation().nearest_power_vertex(*it)->point()
                << "\n";
    }
    if(refine_balls) {
      separate_balls();
    }
  }
};

#include "Polyhedron_demo_remeshing_plugin_protection_cgal_code.h"

void insert_spheres(C2t3& c2t3, Polyhedron* pMesh, const FT size, 
                    bool refine_balls) {

  Insert_spheres go(c2t3, pMesh, size, refine_balls);
}


