#ifndef CGAL_MESH_H
#define CGAL_MESH_H
#include <list>
#include <map>
#include <cmath>
#include <iostream.h>
#include <fstream.h>
//#include <qwidget.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>
//#include <CGAL/Triangulation_vertex_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
//#include <CGAL/predicate_classes_2.h>
//#include <CGAL/Triangulation_circulators_2.h>

//CGAL_BEGIN_NAMESPACE
using namespace CGAL;

#define INVALID_EDGE Edge(NULL, -1);

// #define UPDATE {viewer->repaint(); qApp->processEvents(); float ff; for(float f=0.0; f<1.2e6; f++) ff += sin(f);}
//#define UPDATE

// template <class Gt, class Tds> class Triangulation_finite_vertex_circulator :
// 	public Triangulation_vertex_circulator_2<Gt, Tds>
// {

// };

template <class Gt, class Tds, class Mtraits = void>
class Mesh:
public Constrained_Delaunay_triangulation_2 <Gt, Tds>
{

public:
  //  QWidget *viewer;

  typedef Gt      Triangulation_geom_traits;
  typedef Tds     Triangulation_data_structure;
  typedef typename Gt::FT FT;
  typedef FT      Length;
  typedef FT      Square_length;

  typedef CGAL::Triangulation_2<Gt, Tds>                Tr;
  typedef CGAL::Constrained_triangulation_2<Gt, Tds>   Ct;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds>     CDt;
  //  typedef typename CDt::Constraint             Constraint;
  typedef typename CDt::Vertex                 Vertex;
  typedef typename CDt::Edge                   Edge;
  typedef typename CDt::Face_handle            Face_handle;
  typedef typename CDt::Vertex_handle          Vertex_handle;
  typedef typename CDt::Finite_faces_iterator  Finite_faces_iterator;
  
  typedef std::pair<Vertex_handle,Vertex_handle> Constrained_edge;

  typedef typename Gt::Point_2                        Point;
  typedef std::pair<Point,Point>             Constraint;
  typedef std::list<Constraint>              List_constraints;

  typedef Mtraits  Mesh_traits;
  
  enum Cluster_status {REDUCED,NON_REDUCED};
  struct Cluster {
    Cluster_status status;
    Length length;
    FT alpha_min;
    std::map<Vertex_handle, Length> vertices;
  };

  void write(ostream &f);
  void read(istream &f);
  void reset() {
  cluster_map.clear();
  clear();
}

private:
  int current_time;
  std::multimap<FT, Face_handle>    Bad_faces;
  std::list<Constrained_edge>       c_edge_queue;
  std::map<Vertex_handle, int>      insertion_time;
  std::map<Vertex_handle, FT>       insertion_radius_map;
  std::multimap<Vertex_handle, Cluster>  cluster_map;

public:
  //INSERTION-REMOVAL
  void refine_mesh(/*QWidget *v */);
  void bounds(FT &xmin, FT &ymin, 
	      FT &xmax, FT &ymax,
	      FT &xcenter, FT &ycenter);
  void bounding_box();


    
  Mesh(const Geom_traits& gt=Geom_traits()):CDt(gt){};

    Mesh(List_constraints& lc, const Geom_traits& gt=Geom_traits()):CDt(gt)
    {
      typename List_constraints::iterator lcit=lc.begin();
      for( ; lcit != lc.end(); ++lcit)
	{
	  insert( (*lcit).first, (*lcit).second);
	}
      CGAL_triangulation_postcondition(is_valid());
    }

    template <class InputIterator>
    Mesh(InputIterator first, InputIterator last, const Geom_traits&
	 gt=Gt()):CDt(gt)
    {
      while(first != last){
	insert((*first).first, (*first).second);
	++first;
      }
      CGAL_triangulation_postcondition(is_valid());
    }


 
private:
  void refine_face(Face_handle f);
  void refine_edge(Vertex_handle va, Vertex_handle vb);
  void create_clusters();
  void create_clusters_of_vertex(Vertex_handle v);
  void show_clusters();
  Vertex_handle insert_middle(Vertex_handle va, Vertex_handle vb);
  void cut_cluster(Vertex_handle va, Vertex_handle vb);
  void cut_reduced_cluster(Vertex_handle va, Vertex_handle vb);
  void cut_cluster_edge(Vertex_handle va, Vertex_handle vb);
  Vertex_handle insert_in_c_edge(Vertex_handle va, Vertex_handle vb, Point p);
  void fill_edge_queue();
  void fill_facette_map();
  void process_edge_queue();
  void process_facette_map();
  void update_c_edge_queue(Vertex_handle va,
			   Vertex_handle vb,
			   Vertex_handle vm);
  void update_facette_map(Vertex_handle v);
  void update_cluster(Vertex_handle va, Vertex_handle vb, Vertex_handle vm);
  inline Edge edge_between(Vertex_handle va, Vertex_handle vb);

 

  //CHECK
  bool is_encroached(Vertex_handle va, Vertex_handle vb, Point p);
  bool is_encroached(Vertex_handle va, Vertex_handle vb);
  bool min_insertion_radius(Vertex_handle v, Cluster &c);
  bool is_bad(Face_handle f);
  bool is_small_angle(Vertex_handle vleft, 
		      Vertex_handle vmiddle, 
		      Vertex_handle vright);
  bool is_cluster_reduced(const Cluster&); //look cluster status
  bool is_cluster(Vertex_handle va, Vertex_handle vb);
 
  FT shortest_edge_of_cluster(Vertex_handle v, Cluster &cluster);
  void check_cluster_status( Cluster&); 



  // HELPING functions
  Square_length shortest_edge(Face_handle f);
  FT angle(Vertex_handle vleft, Vertex_handle vmiddle, Vertex_handle vright);
  FT circumradius_to_shortest_edge_ratio(Face_handle f);
  FT insertion_radius(Vertex_handle v);
  Vertex_handle nearest_incident_vertex(Vertex_handle v);
  bool find_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c);
  inline Vertex_circulator incr(Vertex_circulator &c) {
  c++;
  if(Vertex_handle(c) == infinite_vertex()) {
    c++;
  }
  return c;
}

inline Vertex_circulator incr_constraint(Vertex_handle v, Vertex_circulator &c)
{
  Vertex_circulator cbegin = c;
  incr(c);
  while(!edge_between(v, c).first->is_constrained(edge_between(v, c).second) 
	&& c != cbegin ) {
    incr(c);
  }
  return c;
}
	
inline Vertex_circulator succ(Vertex_circulator c) { return incr(c); }
inline Vertex_circulator succ_constraint(Vertex_handle v, 
					 Vertex_circulator c) 
{
  return incr_constraint(v, c); 
}
inline bool get_conflicting_edges(const Point &p,
							  list<Constrained_edge> &cel)
{ //cel = constrained edges list
  Face_handle start;
  int li;
  Locate_type lt;
  Face_handle fh = locate(p,lt,li, start);
  switch(lt) {
  case OUTSIDE_AFFINE_HULL:
  case VERTEX:
    return false;
  case FACE:
  case EDGE:
  case OUTSIDE_CONVEX_HULL:
    propagate_conflicting_edges(p,fh,0,cel);
    propagate_conflicting_edges(p,fh,1,cel);
    propagate_conflicting_edges(p,fh,2,cel);
    return true;    
  }
  CGAL_triangulation_assertion(false);
  return false;
}
//private:


	inline void propagate_conflicting_edges (const Point  &p,
				    Face_handle fh, 
				    int i, 
				    list<Constrained_edge> &cel)
  {
    if(is_infinite(make_pair(fh,i))) return;
    Face_handle fn = fh->neighbor(i);
    Edge e(fh,i);
    if(!is_encroached(e.first->vertex(cw(e.second)), 
		      e.first->vertex(ccw(e.second)), p)) {
      return;
    }
    if ( fh->is_constrained(i)) {// || ! test_conflict(p,fn)) {
      cel.push_back(Constrained_edge(e.first->vertex(cw(e.second)), 
				e.first->vertex(ccw(e.second))));//Edge(fn, fn->index(fh));
      return;
    }
    
    int j = fn->index(fh);
    propagate_conflicting_edges(p,fn,ccw(j),cel);
    propagate_conflicting_edges(p,fn,cw(j),cel);
    return;
  }

bool is_infinite(Face_handle fh) {
  Vertex_handle va, vb, vc;
  va=fh->vertex(0);
  vb=fh->vertex(1);
  vc=fh->vertex(2);
  if(va==infinite_vertex() || vb==infinite_vertex() || 
     vc==infinite_vertex()) return true;
  //check whether the face still exists...
  if(!is_face(va, vb, vc)) return true;
  return false;
}

bool is_infinite(Edge e) {
  Vertex_handle va, vb;
  va=e.first->vertex(cw(e.second));
  vb=e.first->vertex(ccw(e.second));
  if(va==infinite_vertex() || vb==infinite_vertex()) return true;
  //check whether the face still exists...
  if(!is_edge(va, vb)) return true;
  return false;
}

};

template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
find_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c)
{
    // check if vb is in any cluster of va 
  pair<multimap<Vertex_handle, Cluster>::iterator,
    multimap<Vertex_handle, Cluster>::iterator> range = 
    cluster_map.equal_range(va);
  for(multimap<Vertex_handle, Cluster>::iterator it = range.first;
      it != range.second; it++)
    {
      Cluster &cl = it->second;
      if(cl.vertices.find(vb)!=cl.vertices.end()) {
	c = it->second;
	return true;
      }
    }
  return false;
}

//the function that writes a file
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
write(ostream &f)
{
  int nedges = 0;
  Edge_iterator eit = edges_begin();
  while(eit != edges_end()) {
    if((*eit).first->is_constrained((*eit).second))
      nedges++;
    eit++;
  }
  f<<nedges<<endl;
  eit=edges_begin();
  while(eit!=edges_end()) {
    if((*eit).first->is_constrained((*eit).second)) {
      f<<(*eit).first->vertex(cw((*eit).second))->point().x()<<" ";
      f<<(*eit).first->vertex(cw((*eit).second))->point().y()<<" ";
      f<<(*eit).first->vertex(ccw((*eit).second))->point().x()<<" ";
      f<<(*eit).first->vertex(ccw((*eit).second))->point().y()<<endl;
    }
    eit++;
  }
}


//the function that reads a file
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
read(istream &f)
{
  int nedges = 0;
  clear();
  f>>nedges;
  for(int n=0; n<nedges; n++) {
    FT x1, y1, x2, y2;
    f>>x1>>y1>>x2>>y2;
    insert(Point(x1, y1), Point(x2, y2));
  }
}

template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
fill_edge_queue()
{
  Edge_iterator ei = edges_begin();
  while(ei != edges_end()) {
    Vertex_handle va = (*ei).first->vertex(cw((*ei).second));
    Vertex_handle vb = (*ei).first->vertex(ccw((*ei).second));
    if((*ei).first->is_constrained((*ei).second) && 
       is_encroached(va, vb)&&
       va!=infinite_vertex()&&
       vb!=infinite_vertex()) {
      c_edge_queue.push_back(make_pair(va, vb));
    }
    ei++;
  }
}

//it is necessarry for process_facette_map
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
fill_facette_map()
{
  Face_iterator fit = faces_begin();
  while(fit != faces_end())
    {
      if(!is_infinite(fit)) {
	if( is_bad(fit))
	  {
	    Bad_faces.insert(make_pair(
	      circumradius_to_shortest_edge_ratio(fit), fit));
	  }
      }
      fit++;
    }
}


//is used by process_edge_queue
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
process_edge_queue()
{
  while(! c_edge_queue.empty() )
    {
      Constrained_edge ce;
      ce=c_edge_queue.front();
      c_edge_queue.pop_front();
      // if ce is not an edge, choose another. should be more robust
      while(!is_edge(ce.first, ce.second) && !c_edge_queue.empty()) {
				ce=c_edge_queue.front();
				c_edge_queue.pop_front();
      }
      Face_handle fh;
      int i;
      if ( is_edge( ce.first, ce.second, fh,i))
			{
	  		if(fh->is_constrained(i) &&
	     		is_encroached(ce.first, ce.second))
	    	{
		  //Vertex_handle va = ce.first;
		  //Vertex_handle vb = ce.second;
	      	refine_edge(ce.first, ce.second);
	      	//c_edge_queue.pop_front();
	    	}
			}
    }
}




template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
process_facette_map()
{
  if(Bad_faces.empty())
    {
    }
  while( !Bad_faces.empty())
    {
      Face_handle Bf = (*(Bad_faces.begin())).second;
      while(is_infinite(Bf)) {
	Bad_faces.erase(Bad_faces.begin());
	if(Bad_faces.empty()) return;
	Bf = (*(Bad_faces.begin())).second;
      }
      Bad_faces.erase(Bad_faces.begin());
      Vertex_handle va, vb, vc;
      va = Bf->vertex(0);
      vb = Bf->vertex(1);
      vc = Bf->vertex(2);
      if( is_face(va,vb,vc )&&is_bad(Bf))
	{
	  refine_face(Bf);
	}
    }
}


//this function split all the segments that are encroached
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
refine_edge(Vertex_handle va, Vertex_handle vb)
{
  // UPDATE;
  Face_handle f;
  int i;
  is_edge(va, vb, f, i);
  if( !is_cluster(va,vb) && 
      !is_cluster(vb,va) && 
      f->is_constrained(i)) {
    Vertex_handle vm = insert_middle(va,vb);
    update_c_edge_queue(va, vb, vm);
    update_facette_map(vm);
    return;
  }
  else if (is_cluster(va,vb) && is_cluster(vb,va)) {
    Vertex_handle vm = insert_middle(va,vb);
    update_c_edge_queue(va, vb, vm);
    update_facette_map(vm);
    if(vm.is_null()) {
			exit(-1);
		}
    update_cluster(va,vb,vm);
    update_cluster(vb,va,vm);
    return;
  }
  Vertex_handle vaa = is_cluster(va,vb) ? va : vb;
  Vertex_handle vbb = is_cluster(va,vb) ? vb : va;
  cut_cluster_edge(vaa,vbb);
}

 //split all the bad faces
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
refine_face(Face_handle f)
{
  Point pc; 
  Vertex_handle v;
  //check if the face still exists
  pc = circumcenter(f);
  list<Constrained_edge> conflicts;
  list<Constrained_edge>::iterator out_conflicts=c_edge_queue.begin();
  if(get_conflicting_edges(pc, c_edge_queue))
    {
      if(!c_edge_queue.empty()){
	while(! c_edge_queue.empty() )
	  {
	    Constrained_edge ce;
	    ce=c_edge_queue.front();
	    c_edge_queue.pop_front();
	    while(!is_edge(ce.first, ce.second) && !c_edge_queue.empty()) {
	      ce=c_edge_queue.front();
	      c_edge_queue.pop_front();
	    }
	    Face_handle fh;
	    int i;
	    if ( is_edge( ce.first, ce.second, fh,i))
	      {
		if(	fh->is_constrained(i) &&
			is_encroached(ce.first, ce.second))
		  {
		    Vertex_handle va = ce.first;
		    Vertex_handle vb = ce.second;
		    Square_length rg = shortest_edge(f);
		    Cluster cluster;
		    Vertex_handle vaa = is_cluster(va,vb) ? va : vb;
		    Vertex_handle vbb = is_cluster(va,vb) ? vb : va;
		    find_cluster(vaa, vbb, cluster);
		    Length rmin = shortest_edge_of_cluster(vaa, cluster);
		    if(rmin >rg) {
		      refine_edge(ce.first, ce.second);
		    } else {
		      exit(-1);
		    }
		      		//c_edge_queue.pop_front();
		  }
	      }
	  }  
      }
      else
	{
	  Face_handle start;
	  int li;
	  Locate_type lt;
	  Face_handle fh = locate(pc,lt,li, start);
	  if(lt!=OUTSIDE_CONVEX_HULL) {
	    v = insert(pc);
	    update_facette_map(v);
	  }
	}
  }
  // UPDATE; 
  current_time++; 
  insertion_time[v] = current_time; 
  //flip_around(v);  
}


template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
create_clusters()
{
	Vertex_iterator vit = vertices_begin();
	while(vit != vertices_end()) {
		create_clusters_of_vertex(vit);
		vit++;
	}
}

template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
create_clusters_of_vertex(Vertex_handle v)
{
  // prerequisite: at least too vertices must exist in the triangulation
  Vertex_circulator vcirc = incident_vertices(v);
  
  incr_constraint(v, vcirc);
  Vertex_circulator vbegin = vcirc;
  incr_constraint(v, vcirc);
  do {
    if(!is_small_angle(vcirc, v, succ_constraint(v, vcirc))) {
      break;
    }
    incr_constraint(v, vcirc);
  } while(vcirc!=vbegin);
  // now vcirc points to:
  //  * either the point JUST BEFORE the first large angle
  //  * or vbegin if such an angle has not been encountered

  // begin with the second case: a single cluster exists
  if(vcirc==vbegin) {
    FT min_angle = 1000.0; // a BIG value for an angle
    Cluster cl;
    
    do{
      if(is_small_angle(vcirc, v, succ_constraint(v, vcirc))) {
	cl.vertices.insert(make_pair(vcirc, squared_distance(v->point(), 
							     vcirc->point())));
	min_angle = min(min_angle, angle(vcirc, v, succ_constraint(v, vcirc)));
      }
      incr_constraint(v, vcirc);
    }while(vcirc!=vbegin);

    cl.alpha_min = min_angle;
    cluster_map.insert(make_pair(v, cl));
  } else {
    // the first case
    vbegin = vcirc;
    do{
      Cluster cl;
      cl.alpha_min = 1000;
      do {
	cl.vertices.insert(make_pair(vcirc, squared_distance(v->point(),
							     vcirc->point())));
	cl.alpha_min = min(cl.alpha_min, angle(vcirc, v,
					       succ_constraint(v, vcirc)));
	incr_constraint(v, vcirc);
      } while(is_small_angle(vcirc, v, succ_constraint(v, vcirc)));
      cluster_map.insert(make_pair(v, cl));
    }  while(vcirc != vbegin);
  }
}



template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
show_clusters()
{
  multimap<Vertex_handle, Cluster>::iterator cmit = cluster_map.begin();
  while(cmit != cluster_map.end()) {
    Vertex_handle v = (*cmit).first;
    map<Vertex_handle, Length>::iterator cit = (*cmit).second.vertices.begin();
    while(cit != (*cmit).second.vertices.end()) {
      Vertex_handle v1 = (*cit).first;
      cit++;
    }
    cmit++;
  }
}

//refine the cluster
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
cut_cluster(Vertex_handle va, Vertex_handle vb)
{
  Cluster c;
  find_cluster(va, vb, c);
  check_cluster_status(c);
  if( !is_cluster_reduced(c))
    {
      cut_cluster_edge(va, vb);
    }
  else
    {
      cut_reduced_cluster(va, vb);
    }
}


//refine the reduced_cluster
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
cut_reduced_cluster(Vertex_handle va, Vertex_handle vb)
{
  Cluster c;
  find_cluster(va, vb, c);
  Vertex_handle vm;
  if(is_cluster_reduced(c))
    {
      // vm = insert_in_c_edge(va, vb, va->point() + (vb->point() - va->point())/2.0);
    }
  update_c_edge_queue(va, vb, vm);
  update_facette_map(vm);
}


//refine the cluster edges
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
cut_cluster_edge(Vertex_handle va, Vertex_handle vb)
{
  //seulement les constraints
  Cluster c;
  find_cluster(va, vb, c);
  Square_length l2 = 2.0; //shortest_edge_of_cluster(va, c)/1.001;
  Square_length L2 = squared_distance(va->point(), vb->point());
  if(L2 > l2)
    {
      FT rapport = pow(2.0, rint(0.5*log(L2/l2)/log(2.0)))*::sqrt(l2/(L2*4));
      Point pc=va->point()+(vb->point()-va->point())*rapport;
      Vertex_handle vc = insert_in_c_edge(va,vb,pc);
      flip_around(vc);
      update_c_edge_queue(va, vb, vc);
      update_facette_map(vc);
      update_cluster(va, vb, vc);
 		}
}


template <class Gt, class Tds, class Mtraits>
Mesh<Gt, Tds, Mtraits>::Vertex_handle Mesh<Gt, Tds, Mtraits>::
insert_middle(Vertex_handle va, Vertex_handle vb)
{
  Point midpoint;
  midpoint = va->point() + (vb->point() - va->point())/2.0 ;
 return insert_in_c_edge(va,vb, midpoint);
}


//insert in constraint edge the middle
template <class Gt, class Tds, class Mtraits>
Mesh<Gt, Tds, Mtraits>::Vertex_handle Mesh<Gt, Tds, Mtraits>::
insert_in_c_edge(Vertex_handle va, Vertex_handle vb, Point p)
{
  Face_handle f;
  int i;
  is_edge(va, vb, f, i);
  Vertex_handle v = Ct::special_insert_in_edge(p, f, i);
  flip_around(v) ;
  update_c_edge_queue(va, vb, v);
  return v;
}




//update the encroached segments list
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
update_c_edge_queue(Vertex_handle va, Vertex_handle vb, Vertex_handle vm)
{
  //c_edge_queue.pop_back();
  Face_circulator fc = incident_faces(vm);
  Face_circulator fcbegin = fc;
  if(fc.is_empty()) return;
  do {
    int i;
    for(i=0; i<3; i++) {
      Edge e=Edge(fc, i);
      Vertex_handle v1=e.first->vertex(cw(i));
      Vertex_handle v2=e.first->vertex(ccw(i));
      if(is_encroached(v1, v2)&&
	 e.first->is_constrained(i)&&
	 v1!=infinite_vertex()&&v2!=infinite_vertex()) 
	{
	  c_edge_queue.push_back(Constrained_edge(v1, v2));
	}
    }
    fc++;
  } while(fc != fcbegin);
  if(is_encroached(va, vm)) {
    c_edge_queue.push_back(Constrained_edge(va, vm));
  }

  if(is_encroached(vb, vm)) {
    c_edge_queue.push_back(Constrained_edge(vb, vm));
  }
 

}

template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
update_facette_map(Vertex_handle v)
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fcbegin = fc;
  do {
    if(!is_infinite(fc)) {
      if(is_bad(fc)) {
	Bad_faces.insert(make_pair(
	  circumradius_to_shortest_edge_ratio(fc), fc));
      }
    }
    fc++;
  } while(fc!=fcbegin);
}


//ok
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
update_cluster(Vertex_handle va, Vertex_handle vb, Vertex_handle vm)
{
  multimap<Vertex_handle, Cluster>::iterator it_va_cluster =
    cluster_map.find(va);
  int n_va_cluster = cluster_map.count(va);
  for(int i=0; i<n_va_cluster; i++)
    {
      map<Vertex_handle, Length> &vertices = (*it_va_cluster).second.vertices;
      map<Vertex_handle, Length>::iterator it_vertices = vertices.find(vb);
      if(it_vertices != vertices.end())
	{
	  vertices.erase(vb);
	  if(va.is_null()) {}
	  if(vm.is_null()) {}
				
	  vertices[vm] = (squared_distance(va->point(), vm->point()));
	  break;
	}
      it_va_cluster++;
    }

  multimap<Vertex_handle, Cluster>::iterator it_vb_cluster =
    cluster_map.find(vb);
  int n_vb_cluster = cluster_map.count(vb);
  for(int i=0; i<n_vb_cluster; i++)
    {
      map<Vertex_handle, Length> &vertices = (*it_vb_cluster).second.vertices;
      map<Vertex_handle, Length>::iterator it_vertices = vertices.find(va);
      if(it_vertices != vertices.end())
	{
	  vertices.erase(va);
	  vertices[vm] = (squared_distance(vb->point(), vm->point()));
	  break;
	}
      it_vb_cluster++;
    }
}

template <class Gt, class Tds, class Mtraits>
inline Mesh<Gt, Tds, Mtraits>::Edge 
Mesh<Gt, Tds, Mtraits>::edge_between(Vertex_handle va, Vertex_handle vb) {
  Edge_circulator ec = va->incident_edges();
  Edge_circulator ecbegin = ec;
  do {
    Edge e = (*ec);
    Face_handle f = e.first;
    int iedge = e.second;
    Vertex_handle v1 = f->vertex(cw(iedge));
    Vertex_handle v2 = f->vertex(ccw(iedge));
    if((v1 == va && v2 == vb) || (v2 == va && v1 == vb)) {
      return e;
    }
    ec++;
  } while(ec != ecbegin);
  cerr<<"invalid edge"<<endl;
  return INVALID_EDGE;
}




//CHECK


template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_encroached(Vertex_handle va, Vertex_handle vb, Point p)
{
   Point pm=midpoint(va->point(), vb->point());
   if(2.0*::sqrt(squared_distance(pm, p)) < ::sqrt(squared_distance(va->point(), vb->point())))
     {
       return true;
     }
   return false;
}


template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_encroached(Vertex_handle va, Vertex_handle vb)
{
  Vertex_iterator vi=vertices_begin();
  while(vi!=vertices_end())
    {
      if(is_encroached(va, vb, vi->point())&&
	 va!=Vertex_handle(vi)&&
	 vb!=Vertex_handle(vi))
	{
	  return true;
	}
      vi++;
    }
  return false;

}


template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
min_insertion_radius(Vertex_handle v, Cluster &c)
{

}


//the measure of faces quality
template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_bad(Face_handle f)
{
  FT quality = circumradius_to_shortest_edge_ratio(f);
    if((quality >1) || (quality < 0.5))
      //set_a:(>1.0 || <0.4)
      //set_b:(>1.0 || <0.5)
      //set_c:(>1.1 || <0.5)
      //set_d:(>1.1 || <0.4)
      //set_e:(>0.9 || <0.5)
      //set_f:(>0.9 || <0.4)//~NOK
      //set_g:(>0.8 || <0.4)//NOK
      //set_h:(>0.8 || <0.5)
      {
	return TRUE;
      }
    else
      {
	return FALSE;
      }
    // }
}


//
template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_small_angle(Vertex_handle vleft,
	       Vertex_handle vmiddle,
	       Vertex_handle vright)
{
  FT cos_alpha;
  cos_alpha = angle(vleft, vmiddle, vright);
  if(cos_alpha > 1/2)
    {
      return TRUE; //the same cluster
    }
  else
    {
      return FALSE; //another cluster
    }
}


//
template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_cluster_reduced(const Cluster& c)
{
  return c.status == REDUCED;
}



template <class Gt, class Tds, class Mtraits>
bool Mesh<Gt, Tds, Mtraits>::
is_cluster(Vertex_handle va, Vertex_handle vb)
{
  // check if vb is in any cluster of va 
  pair<multimap<Vertex_handle, Cluster>::iterator,
    multimap<Vertex_handle, Cluster>::iterator> range = 
    cluster_map.equal_range(va);
  for(multimap<Vertex_handle, Cluster>::iterator it = range.first;
      it != range.second; it++)
    {
      Cluster &cl = it->second;
      if(cl.vertices.find(vb)!=cl.vertices.end()) {
	return true;
      }
    }
  return false;
}


template <class Gt, class Tds, class Mtraits>
typename Gt::FT Mesh<Gt, Tds, Mtraits>::
shortest_edge_of_cluster(Vertex_handle v, Cluster &cluster)
{ 
  map<Vertex_handle, Length>::iterator vit = cluster.vertices.begin();
  FT min_edge = squared_distance(((*vit).first)->point(), (*v).point());
  vit++;
  while(vit != cluster.vertices.end())
    {
      FT temp = squared_distance(((*vit).first)->point(), (*v).point());
      if(temp < min_edge)
	{
	  min_edge = temp;
	}
     vit++;
    }
  return min_edge;
}

// look if all edges have the same lenght
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>:: 
check_cluster_status( Cluster& cluster)
{
  Length initl;
  map<Vertex_handle, Length>::iterator vi = cluster.vertices.begin();
  initl = (*(cluster.vertices.find((*vi).first))).second;
  vi++;
  for(;vi != cluster.vertices.end(); vi++)
    {
      if(( (*(cluster.vertices.find((*vi).first))).second)  != initl)
	{
	  cluster.status = NON_REDUCED;

	  return;
	}
      vi++;  
    }
  cluster.status = REDUCED;
}









template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
bounds(FT &xmin, FT &ymin, 
       FT &xmax, FT &ymax,
       FT &xcenter, FT &ycenter)
{
  Vertex_iterator vi=vertices_begin();
  xmin=xmax=xcenter=vi->point().x();
  ymin=ymax=ycenter=vi->point().y();
  vi++;
  while(vi != vertices_end())
  {
    xcenter+=vi->point().x();
    ycenter+=vi->point().y();
    if(vi->point().x() < xmin) xmin=vi->point().x();
    if(vi->point().x() > xmax) xmax=vi->point().x();
    if(vi->point().y() < ymin) ymin=vi->point().y();
    if(vi->point().y() > ymax) ymax=vi->point().y();
    vi++;
  }
  xcenter /= number_of_vertices();
  ycenter /= number_of_vertices();
}


//bounding box
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
bounding_box()
{
 
  FT span;
  FT xmin, xmax, ymin, ymax, xcenter, ycenter;
  bounds(xmin, ymin, xmax, ymax, xcenter, ycenter);
  span = max((xmax-xmin), (ymax-ymin));
  Point bb1(xcenter - 1.5*span, ycenter - 1.5*span);
  Point bb2(xcenter + 1.5*span, ycenter - 1.5*span);
  Point bb3(xcenter + 1.5*span, ycenter + 1.5*span);
  Point bb4(xcenter - 1.5*span, ycenter + 1.5*span);
  insert(bb1);
  insert(bb2);
  insert(bb3);
  insert(bb4);
  insert(bb1, bb2);
  insert(bb2, bb3);
  insert(bb3, bb4);
  insert(bb4, bb1);

}


//the angle that are between 2 edges from the triangulation
template <class Gt, class Tds, class Mtraits>
typename Gt::FT Mesh<Gt, Tds, Mtraits>::
angle(Vertex_handle vleft, Vertex_handle vmiddle, Vertex_handle vright)
{
  Point 
    pa = vleft->point(),
    pb = vmiddle->point(),
    pc = vright->point();
  FT a, b, c;
  a = ::sqrt(squared_distance(pb, pc));
  b = ::sqrt(squared_distance(pc, pa));
  c = ::sqrt(squared_distance(pa, pb));
  FT cos_alpha = ((a*a+c*c-b*b)/(2*a*c));
  return cos_alpha;
}


//the shortest edge that are in a triangle
template <class Gt, class Tds, class Mtraits>
typename Gt::FT Mesh<Gt, Tds, Mtraits>::
shortest_edge(Face_handle f)
{
  Point 
    pa = (f->vertex(0))->point(),
    pb = (f->vertex(1))->point(),
    pc = (f->vertex(2))->point();
  FT a, b, c;
  a =(squared_distance(pb, pc));
  b =(squared_distance(pc, pa));
  c =(squared_distance(pa, pb));
  return (min(a, min(b, c))); // regarder la documentation du CGAL sur
  // MIN (Developer Manual)
}


//the triangle quality is represented by the
//circumradius_to_shortest_edge_ratio value
template <class Gt, class Tds, class Mtraits>
typename Gt::FT Mesh<Gt, Tds, Mtraits>::
circumradius_to_shortest_edge_ratio(Face_handle f)
{
  Point p;
  p = circumcenter(f);
  Point A = (f->vertex(0))->point();
  FT radius =::sqrt(squared_distance(p, A));
  FT sh_edge =::sqrt(shortest_edge(f));
  return (radius/sh_edge);
}


//insertion radius: the definition
template <class Gt, class Tds, class Mtraits>
typename Gt::FT Mesh<Gt, Tds, Mtraits>::
insertion_radius(Vertex_handle v)
{
  Vertex_handle v1 = nearest_incident_vertex(v);
  return (squared_distance(v->point(), v1->point()));
}


//this function must compute the vertex that are so close to our vertex.
template <class Gt, class Tds, class Mtraits>
Mesh<Gt, Tds, Mtraits>::Vertex_handle Mesh<Gt, Tds, Mtraits>::
nearest_incident_vertex(Vertex_handle v)
{
  Vertex_handle vbegin, vcurrent, vnearest;
  Vertex_circulator circ = v->incident_vertices();
  vbegin = circ;
  vnearest = vbegin;
  circ++;
  FT dist = (squared_distance(v->point(), vbegin->point()));
  while((vcurrent = circ) != vbegin) {
    FT d1 = (squared_distance(v->point(), vcurrent->point()));
    if(d1 < dist) {
      vnearest = vcurrent;
    }
    circ++;
  }
  return vnearest;
}


//the mesh refine function 
template <class Gt, class Tds, class Mtraits>
void Mesh<Gt, Tds, Mtraits>::
refine_mesh(/* QWidget *w */)
{
  //  viewer = w;
  // bounding_box();
  create_clusters();
  show_clusters();

  fill_edge_queue();
  fill_facette_map();
  process_edge_queue();
  process_facette_map();
  
}

//CGAL_END_NAMESPACE


#endif

