#ifndef CGAL_MESH_H
#define CGAL_MESH_H
#include <CGAL/basic.h>
#include <list>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

CGAL_BEGIN_NAMESPACE

template <class Tr, class Mtraits = void>
class Mesh: public Tr
{

public:
  typedef Tr Triangulation;
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Geom_traits::FT FT;
  typedef FT      Length;
  typedef FT      Square_length;

  typedef typename Tr::Vertex                 Vertex;
  typedef typename Tr::Edge                   Edge;
  typedef typename Tr::Edge_iterator          Edge_iterator;
  typedef typename Tr::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Tr::Edge_circulator        Edge_circulator;
  typedef typename Tr::Face_handle            Face_handle;
  typedef typename Tr::Face_iterator          Face_iterator;
  typedef typename Tr::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Tr::Face_circulator        Face_circulator;
  typedef typename Tr::Vertex_handle          Vertex_handle;
  typedef typename Tr::Vertex_iterator        Vertex_iterator;
  typedef typename Tr::Vertex_circulator      Vertex_circulator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Tr::Locate_type            Locate_type;
  
  typedef std::pair<Vertex_handle,Vertex_handle> Constrained_edge;

  typedef typename Tr::Point                  Point;

  typedef typename Tr::Constraint             Constraint;
  typedef typename Tr::List_constraints       List_constraints;

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
  // Bad_faces: list of bad finite faces
  // # warning: some faces could be recycled during insertion in the
  //  triangulation
  // # TODO: create a const iterator that give real bad faces
  std::multimap<FT, Face_handle>    Bad_faces;

  // c_edge_queue: list of encroached constrained edges
  //  warning: some edges could be destroyed
  // # TODO: idem const iterator
  std::list<Constrained_edge>       c_edge_queue;

  //  std::map<Vertex_handle, int>      insertion_time;
  //  std::map<Vertex_handle, FT>       insertion_radius_map;
  std::multimap<Vertex_handle, Cluster>  cluster_map;

public:
  //INSERTION-REMOVAL
  void refine_mesh();
  // TODO: refine_mesh_step(), that do a step of refinment
  void bounds(FT &xmin, FT &ymin, 
	      FT &xmax, FT &ymax,
	      FT &xcenter, FT &ycenter);
  void bounding_box();

  Mesh():Tr(){};

    Mesh(List_constraints& lc, const Geom_traits& gt=Geom_traits()):Tr(gt)
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
						  gt=Geom_traits()):Tr(gt)
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
  //  void show_clusters();
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
  //  bool min_insertion_radius(Vertex_handle v, Cluster &c);
  bool is_bad(Face_handle f);
  bool is_small_angle(Vertex_handle vleft, 
		      Vertex_handle vmiddle, 
		      Vertex_handle vright);
  bool is_cluster_reduced(const Cluster&); //look cluster status
  bool is_cluster(Vertex_handle va, Vertex_handle vb);
 
  FT shortest_edge_of_cluster(Vertex_handle v, Cluster &cluster);
  void check_cluster_status( Cluster&); 



  // HELPING functions
  Square_length shortest_edge_squared_lenght(Face_handle f);
  FT cosinus_of_angle(Vertex_handle vleft, Vertex_handle vmiddle, Vertex_handle vright);
  FT aspect_ratio(Face_handle f); // r^2/l^2
  //  FT insertion_radius(Vertex_handle v);
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
  // TODO: create a finite_contraints_circulator
  Vertex_circulator cbegin = c;
  incr(c);
  while(!edge_between(v, c).first->is_constrained(edge_between(v, c).second) 
	&& c != cbegin ) {
    incr(c);
  }
  return c;
}
	
inline Vertex_circulator succ_constraint(Vertex_handle v, 
					 Vertex_circulator c) 
{
  return incr_constraint(v, c); 
}
inline bool get_conflicting_edges(const Point &p,
				  list<Constrained_edge> &cel)
{ //cel = constrained edges list
  int li;
  Locate_type lt;
  Face_handle fh = locate(p,lt,li);
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

}; // end of Mesh

template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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

template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
fill_edge_queue()
{
  for(Finite_edges_iterator ei = finite_edges_begin();
      ei != finite_edges_end();
      ei++)
    {
      Vertex_handle va = (*ei).first->vertex(cw((*ei).second));
      Vertex_handle vb = (*ei).first->vertex(ccw((*ei).second));
      if((*ei).first->is_constrained((*ei).second) && 
	 is_encroached(va, vb))
	c_edge_queue.push_back(make_pair(va, vb));
    }
}

//it is necessarry for process_facette_map
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
fill_facette_map()
{
  for(Finite_faces_iterator fit = finite_faces_begin();
      fit != finite_faces_end();
      fit++)
    {
      if( is_bad(fit))
	Bad_faces.insert(make_pair(
	   aspect_ratio(fit), fit));
    }
}


//is used by process_edge_queue
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
process_edge_queue()
{
  while(! c_edge_queue.empty() )
    {
      Constrained_edge ce;
      do {// if ce is not an edge, choose another. should be more robust
	ce=c_edge_queue.front();
	c_edge_queue.pop_front();
      }
      while(!is_edge(ce.first, ce.second) && !c_edge_queue.empty());

      Face_handle fh;
      int i;
      if ( is_edge( ce.first, ce.second, fh,i))
	{
	  if(fh->is_constrained(i) &&
	     is_encroached(ce.first, ce.second))
	    refine_edge(ce.first, ce.second);
	}
    }
}

template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
process_facette_map()
{
  while( !Bad_faces.empty())
    {
      Face_handle Bf = (*(Bad_faces.begin())).second;
//       while(!is_face(Bf)) {
// 	Bad_faces.erase(Bad_faces.begin());
// 	if(Bad_faces.empty()) return;
// 	Bf = (*(Bad_faces.begin())).second;
//       }
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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
    // cluster at both ends
    Vertex_handle vm = insert_middle(va,vb);
    update_c_edge_queue(va, vb, vm);
    update_facette_map(vm);

    CGAL_assertion(!vm.is_null());

    update_cluster(va,vb,vm);
    update_cluster(vb,va,vm);
    return;
  }
  Vertex_handle vaa = is_cluster(va,vb) ? va : vb;
  Vertex_handle vbb = is_cluster(va,vb) ? vb : va;
  cut_cluster_edge(vaa,vbb);
}

 //split all the bad faces
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
refine_face(Face_handle f)
{
  Point pc; 
  Vertex_handle v;
  // //check if the face still exists
  pc = circumcenter(f);
  //  list<Constrained_edge> conflicts;
  //  list<Constrained_edge>::iterator out_conflicts=c_edge_queue.begin();
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
		    Square_length rg = shortest_edge_squared_lenght(f);
		    Cluster cluster;
		    Vertex_handle vaa = is_cluster(va,vb) ? va : vb;
		    Vertex_handle vbb = is_cluster(va,vb) ? vb : va;
		    find_cluster(vaa, vbb, cluster);
		    Length rmin = shortest_edge_of_cluster(vaa,
							   cluster);


		    // PROBLEM!!
		    CGAL_assertion(rmin <= rg);

		    refine_edge(ce.first, ce.second);

		  }
		      		//c_edge_queue.pop_front();
	      }
	  }  
      }
      else
	{
 	  int li;
 	  Locate_type lt;
 	  /*Face_handle fh = */locate(pc,lt,li);
	  if(lt!=OUTSIDE_CONVEX_HULL) {
	    v = insert(pc);
	    update_facette_map(v);
	  }
	}
  }
  // UPDATE; 
  //  current_time++; 
  //  insertion_time[v] = current_time; 
  //flip_around(v);  
}


template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
create_clusters()
{
  for(Vertex_iterator vit = vertices_begin();
      vit != vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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
	min_angle = min(min_angle, cosinus_of_angle(vcirc, v, succ_constraint(v, vcirc)));
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
	cl.alpha_min = min(cl.alpha_min, cosinus_of_angle(vcirc, v,
					       succ_constraint(v, vcirc)));
	incr_constraint(v, vcirc);
      } while(is_small_angle(vcirc, v, succ_constraint(v, vcirc)));
      cluster_map.insert(make_pair(v, cl));
    }  while(vcirc != vbegin);
  }
}



// template <class Tr, class Mtraits>
// void Mesh<Tr, Mtraits>::
// show_clusters()
// {
//   multimap<Vertex_handle, Cluster>::iterator cmit = cluster_map.begin();
//   while(cmit != cluster_map.end()) {
//     Vertex_handle v = (*cmit).first;
//     map<Vertex_handle, Length>::iterator cit = (*cmit).second.vertices.begin();
//     while(cit != (*cmit).second.vertices.end()) {
//       Vertex_handle v1 = (*cit).first;
//       cit++;
//     }
//     cmit++;
//   }
// }

//refine the cluster
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
cut_reduced_cluster(Vertex_handle va, Vertex_handle vb)
{
  Cluster c;
  find_cluster(va, vb, c);
  Vertex_handle vm;
  if(is_cluster_reduced(c))
    {
      // WHY IS IT COMMENTED??
      // vm = insert_in_c_edge(va, vb, va->point() + (vb->point() - va->point())/2.0);
    }
  update_c_edge_queue(va, vb, vm);
  update_facette_map(vm);
}


//refine the cluster edges
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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
      //      flip_around(vc);
      update_c_edge_queue(va, vb, vc);
      update_facette_map(vc);
      update_cluster(va, vb, vc);
 		}
}


template <class Tr, class Mtraits>
Mesh<Tr, Mtraits>::Vertex_handle Mesh<Tr, Mtraits>::
insert_middle(Vertex_handle va, Vertex_handle vb)
{
  Point midpoint;
  midpoint = va->point() + (vb->point() - va->point())/2.0 ;
 return insert_in_c_edge(va,vb, midpoint);
}


//insert in constraint edge the middle
template <class Tr, class Mtraits>
Mesh<Tr, Mtraits>::Vertex_handle Mesh<Tr, Mtraits>::
insert_in_c_edge(Vertex_handle va, Vertex_handle vb, Point p)
{
  Face_handle f;
  int i;
  is_edge(va, vb, f, i);
  Vertex_handle v = special_insert_in_edge(p, f, i);
  // WARNING: special_insert_in_edge is not robust!
  // We should deconstrained the constrained edge, 
  // 
  //  flip_around(v) ; // no longer necessar, because
  //  special_insert_in_edge does it for Constrained_Delaunay_triangulation_2
  update_c_edge_queue(va, vb, v);
  return v;
}




//update the encroached segments list
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
update_c_edge_queue(Vertex_handle va, Vertex_handle vb, Vertex_handle vm)
{
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

template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
update_facette_map(Vertex_handle v)
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fcbegin = fc;
  do {
    if(!is_infinite(fc)) {
      if(is_bad(fc)) {
	Bad_faces.insert(make_pair(
	  aspect_ratio(fc), fc));
      }
    }
    fc++;
  } while(fc!=fcbegin);
}


//ok
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
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

	  CGAL_assertion(!va.is_null() && !vm.is_null());
				
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

// == is_edge ???????????????
// used by: incr_constraint
template <class Tr, class Mtraits>
inline Mesh<Tr, Mtraits>::Edge 
Mesh<Tr, Mtraits>::edge_between(Vertex_handle va, Vertex_handle vb) {
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
  CGAL_assertion(false); // invalid edge
}

//CHECK


// TO IMPLEMENT WITH A SINGLE SCALAR PRODUCT!
// ->traits
template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
is_encroached(Vertex_handle va, Vertex_handle vb, Point p)
{
   Point pm=midpoint(va->point(), vb->point());
   if(2.0*::sqrt(squared_distance(pm, p)) < ::sqrt(squared_distance(va->point(), vb->point())))
     {
       return true;
     }
   return false;
}

// NOT ALL VERTICES, ONLY THE TWO NEIGHBORS
template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
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

// ?????????????
// ->traits
//the measure of faces quality
template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
is_bad(Face_handle f)
{
  FT quality = aspect_ratio(f); // recall: r^2/l^2
  //    if((quality >1) || (quality < 0.5))
  if (quality > 3) // $ B=\sqrt{3} $
      //set_a:(>1.0 || <0.4)
      //set_b:(>1.0 || <0.5)
      //set_c:(>1.1 || <0.5)
      //set_d:(>1.1 || <0.4)
      //set_e:(>0.9 || <0.5)
      //set_f:(>0.9 || <0.4)//~NOK
      //set_g:(>0.8 || <0.4)//NOK
      //set_h:(>0.8 || <0.5)
      {
	return true;
      }
    else
      {
	return false;
      }
    // }
}


// -> traits?
template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
is_small_angle(Vertex_handle vleft,
	       Vertex_handle vmiddle,
	       Vertex_handle vright)
{
  FT cos_alpha = cosinus_of_angle(vleft, vmiddle, vright);
  if(cos_alpha > 1/2)
    {
      return true; //the same cluster
    }
  else
    {
      return false; //another cluster
    }
}


// # used by: cut_cluster, cut_reduced_cluster
template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
is_cluster_reduced(const Cluster& c)
{
  return c.status == REDUCED;
}



template <class Tr, class Mtraits>
bool Mesh<Tr, Mtraits>::
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


template <class Tr, class Mtraits>
typename Mesh<Tr, Mtraits>::FT Mesh<Tr, Mtraits>::
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>:: 
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










// ->traits?
//the angle that are between 2 edges from the triangulation
template <class Tr, class Mtraits>
typename Mesh<Tr, Mtraits>::FT Mesh<Tr, Mtraits>::
cosinus_of_angle(Vertex_handle vleft, Vertex_handle vmiddle, Vertex_handle vright)
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

// ->traits?
//the shortest edge that are in a triangle
// # used by: refine_face, aspect_ratio
template <class Tr, class Mtraits>
typename Mesh<Tr, Mtraits>::FT Mesh<Tr, Mtraits>::
shortest_edge_squared_lenght(Face_handle f)
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

// ->traits
//the triangle quality is represented by the
//aspect_ratio value
// # used by: fill_facette_map, update_facette_map, is_bad
template <class Tr, class Mtraits>
typename Mesh<Tr, Mtraits>::FT Mesh<Tr, Mtraits>::
aspect_ratio(Face_handle f)
{
  Point p;
  p = circumcenter(f);
  Point A = (f->vertex(0))->point();
  FT radius = squared_distance(p, A);
  FT sh_edge = shortest_edge_squared_lenght(f);
  return (radius/sh_edge);
}


// //insertion radius: the definition
// template <class Tr, class Mtraits>
// typename Mesh<Tr, Mtraits>::FT Mesh<Tr, Mtraits>::
// insertion_radius(Vertex_handle v)
// {
//   Vertex_handle v1 = nearest_incident_vertex(v);
//   return (squared_distance(v->point(), v1->point()));
// }


//this function must compute the vertex that are so close to our vertex.
template <class Tr, class Mtraits>
Mesh<Tr, Mtraits>::Vertex_handle Mesh<Tr, Mtraits>::
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
template <class Tr, class Mtraits>
void Mesh<Tr, Mtraits>::
refine_mesh()
{
  //  viewer = w;
  //  bounding_box();
  create_clusters();

  fill_edge_queue();
  fill_facette_map();
  process_edge_queue();
  process_facette_map();
  
}

CGAL_END_NAMESPACE


#endif

