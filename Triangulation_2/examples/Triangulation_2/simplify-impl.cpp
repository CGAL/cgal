

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constraint_hierarchy_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/Timer.h>
#include <list>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

typedef CGAL::Timer Timer;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_predicates_tag                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
typedef CGAL::Polyline_constrained_triangulation_2<CDT>       PCT;
typedef PCT::Point                                    Point;
typedef PCT::Constraint_id                            Cid;
typedef PCT::Vertices_in_constraint_iterator          Vertices_in_constraint_iterator;
typedef PCT::Constraint_iterator                      Constraint_iterator;
typedef PCT::Vertex_handle                            Vertex_handle;
typedef PCT::Vertex_circulator                         Vertex_circulator;
typedef PCT::Face_handle                              Face_handle;
typedef PCT::Edge                                     Edge;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost<double> Cost;

 
  struct Compare_cost 
  { 
    bool operator() ( Vertices_in_constraint_iterator const& x, 
                      Vertices_in_constraint_iterator const& y ) const 
    { 
      return x->cost < y->cost; 
    } 
  } ;
  
  struct Id_map : public boost::put_get_helper<std::size_t, Id_map>
  { 
    typedef boost::readable_property_map_tag category;
    typedef std::size_t                      value_type;
    typedef value_type                       reference;
    typedef Vertices_in_constraint_iterator  key_type;
    
    reference operator[] ( key_type const& x ) const { return x->id ; }
  } ;
  
typedef CGAL::Modifiable_priority_queue<Vertices_in_constraint_iterator,Compare_cost,Id_map> MPQ ;
  

Vertices_in_constraint_iterator
decrement(Vertices_in_constraint_iterator it)
{
  do{ 
    --it;
  } while(it->removed);
 return it;
}

Vertices_in_constraint_iterator
increment(Vertices_in_constraint_iterator it)
{
 do{ 
   ++it;
 } while(it->removed);
 return it;
}


void
dump(const PCT& pct)
{
  Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
  for(; cit!=e; ++cit){
    Cid cid = *cit;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        it++){
      std::cerr << it->point; 
      if(it->fixed) std::cerr << " f";
      if(it->removed) std::cerr << " r";
      if(! it->removed && it->point != it->vertex->point()){ std::cerr << "   " << it->vertex->point();}
      std::cerr << std::endl;
    }
  }
}


void 
write_poly(const PCT& pct, Cid cid)
{
  std::cout << "POLYGON" << std::endl;
  for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
      it != pct.vertices_in_constraint_end(cid);
      it++){
    
    if(! it->removed){
      std::cout << it->point; 
    }
    std::cout << std::endl;
  }
}

void
write_poly(const PCT& pct)
{
  Constraint_iterator b = pct.constraints_begin(), e = pct.constraints_end();
  for(; b!=e; ++b){
    Cid cid = *b;
    write_poly(pct,cid);
  }
}

// Fix the leftmost, rightmost, topmost and bottommost vertex
void
fix_lrbt_vertices(const PCT& pct)
{
  Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
  for(; cit!=e; ++cit){
    Cid cid = *cit;
    Vertices_in_constraint_iterator l,r,b,t;
    l = r = b = t = pct.vertices_in_constraint_begin(cid);
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        it++){ 
      if(it->point.x() < l->point.x()) l = it; 
      if(it->point.x() > r->point.x()) r = it; 
      if(it->point.y() < b->point.y()) b = it; 
      if(it->point.y() > t->point.y()) t = it; 
    }
    l->fixed = r->fixed = t->fixed = b->fixed = true;
  }
}

// For all polyline constraints we compute the cost of all non fixed and not removed vertices
void
initialize_costs(PCT& pct, MPQ& mpq)
{
  Cost cost;
  Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
  for(; cit!=e; ++cit){
    Cid cid = *cit;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      if(! it->fixed  && ! it->removed){
        Vertices_in_constraint_iterator u = decrement(it);
        Vertices_in_constraint_iterator w = increment(it);

        boost::optional<double> dist = cost(pct, u, it, w);
        if(dist){
          it->cost = *dist;
          mpq.push(it);
        } else {
          it->cost = (std::numeric_limits<double>::max)();
        } 
      }
    }
  }
}


bool
is_removable(const PCT& pct, Vertices_in_constraint_iterator it)
{
  typedef PCT::Geom_traits Geom_traits;
  if( it->removed || it->fixed){
    return false;
  }
  
  Vertex_handle vh = it->vertex;
  Vertices_in_constraint_iterator u = decrement(it);
  Vertex_handle uh = u->vertex;
  Vertices_in_constraint_iterator w = increment(it);
  Vertex_handle wh = w->vertex;

  Geom_traits::Orientation_2 orientation_2 = pct.geom_traits().orientation_2_object();
  CGAL::Orientation o = orientation_2(uh->point(), vh->point(), wh->point());
  if(o == CGAL::COLLINEAR){
    return true;
  }
  if(o == CGAL::LEFT_TURN){
    std::swap(uh,wh);
  }

  // uh, vh, wh perform a right turn 
  const Point& up = uh->point();
  const Point& wp = wh->point();
  Vertex_circulator circ = pct.incident_vertices(vh);
  while(circ != uh){
    ++circ;
  }
  ++circ;
  while(circ != wh){
    o = orientation_2(up, circ->point(), wp);
    if(orientation_2(up, wp, circ->point()) != CGAL::RIGHT_TURN){
      return false;
    }
    ++circ;
  }
  return true;
}

// If several polylines pass through the same vertex we have to mark the node
// in the polyline  constraint as fixed
// At the same time ix the vertices with the smallest and largest x and y coordinates
void
fix_shared_vertices(const PCT& pct)
{ 
  CGAL::Unique_hash_map<Vertex_handle,Vertices_in_constraint_iterator> M;

  Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
  for(; cit!=e; ++cit){
    Cid cid = *cit;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      if(! M.is_defined(it->vertex)){
        M[it->vertex]=it;
      } else {
        it->fixed = true;
      } 
    }
  }

}


int
initialize_indices(PCT& pct)
{
  int id = 0;
  Constraint_iterator b = pct.constraints_begin(), e = pct.constraints_end();
  for(; b!=e; ++b){
    Cid cid = *b;
    for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        ++it){
      it->id = id++;
    }
  }
  return id;
}

bool
simplify(PCT& pct,MPQ& mpq)
{
  Vertices_in_constraint_iterator v = mpq.top();
  mpq.pop();
  if(is_removable(pct,v)){
    Vertices_in_constraint_iterator u = decrement(v), w = increment(v);
    pct.simplify(u,v,w);
    Cost cost;
    if(! u->fixed){
      Vertices_in_constraint_iterator uu = decrement(u);
      boost::optional<double> dist = cost(pct, uu,u,w);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        u->cost = *dist;
        mpq.update(u, true);
      }
    }
    
    if(! w->fixed){
      Vertices_in_constraint_iterator ww = increment(w);
      boost::optional<double> dist = cost(pct, u,w,ww);
      if(! dist){
        std::cerr << "undefined cost not handled yet" << std::endl;
      } else {
        w->cost = *dist;
        mpq.update(w, true);
      }
    }
    return true;
  }
  return false;
}


void 
simplify(PCT& pct, int n)
{  
  int m = initialize_indices(pct);
  Compare_cost cc;
  Id_map idm;
  MPQ mpq(m, cc, idm);

  initialize_costs(pct, mpq);
  while(! mpq.empty() && n--){
    simplify(pct,mpq);
  } 
}

int
main( )
{
  Timer timer;
  int n=0, m=0;
  PCT pct;
  std::vector<Point> poly ;
#if 0  
  Point p;
  std::cin >> n;
  for(int i=0;i<n;i++){
    std::cin >> p;
    poly.push_back(p);
  }
  pct.insert_polyline(poly.begin(), poly.end());

  fix_lrbt_vertices(pct);
  fix_shared_vertices(pct);
  simplify(pct, n/2);

  write_poly(pct);

#else
  std::string line;
  std::string POLYGON("POLYGON");
  while(std::getline(std::cin, line)){
    std::string::size_type idx = line.find(POLYGON);
    if(idx == 0){
      if(!poly.empty()){
        pct.insert_polyline(poly.begin(), poly.end());
        m+= poly.size();
        poly.clear();
        n++;
      }
    } else {
      double x,y ;
      std::istringstream ss(line);
      ss >> x >> y ;
      poly.push_back( Point(x,y) );
    }
  }
  if(!poly.empty()){
    pct.insert_polyline(poly.begin(), poly.end());
    m+= poly.size();
    poly.clear();
    n++;
  }
  std::cerr << n << " polylines with in total " << m << " vertices" << std::endl;
  timer.start();
  fix_lrbt_vertices(pct);
  fix_shared_vertices(pct);
  simplify(pct, m/2);
  timer.stop();
  write_poly(pct);
#endif

  std::cerr << timer.time() << "sec." << std::endl;

  return 0;
}


