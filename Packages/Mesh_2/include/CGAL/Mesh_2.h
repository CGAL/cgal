#ifndef CGAL_MESH_H
#define CGAL_MESH_H
#include <CGAL/basic.h>
#include <list>
#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Filtred_container.h>
#include <CGAL/Filtred_circulator.h>

CGAL_BEGIN_NAMESPACE

// auxiliary classes
// # TODO: use this Filtred_iterator usable by Filtred_container, so
// that it erase bad elements. It will need a pointer to the container.
template <class In, class Pred>
class Filtred_iterator_from_container : public In
{
  Pred test;
  In _end;
public:
  typedef Pred Predicate;
  typedef In InputIterator;
  typedef Filtred_iterator_from_container<In,Pred> Self;

  Filtred_iterator_from_container(In current, In end, Pred p=Pred())
    : test(p), _end(end) 
    {
      while(!(*this==_end) & !test(*this))
	this->In::operator++();
    };

  Self& operator++() {
    do {
      this->In::operator++();
    } while(!(*this==_end) & !test(*(*this)));
    return *this;
  }

  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
};


// Tr is a Delaunay constrained triangulation (with intersections or not)
template <class Tr>
class Mesh_2: public Tr
{
public:
  typedef Tr Triangulation;
  typedef Mesh_2<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Geom_traits::FT FT;
  typedef FT      Squared_length;

  typedef typename Tr::Vertex                 Vertex;
  typedef typename Tr::Edge                   Edge;
  typedef typename Tr::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Tr::Edge_circulator        Edge_circulator;
  typedef typename Tr::Face_handle            Face_handle;
  typedef typename Tr::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Tr::All_faces_iterator  All_faces_iterator;
  typedef typename Tr::Face_circulator        Face_circulator;
  typedef typename Tr::Vertex_handle          Vertex_handle;
  typedef typename Tr::Vertex_circulator      Vertex_circulator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Tr::Locate_type            Locate_type;
  
  typedef std::pair<Vertex_handle,Vertex_handle>
                                              Constrained_edge;

  typedef typename Tr::Point                  Point;

  typedef typename Tr::List_constraints       List_constraints;

public:
  void write(std::ostream &f);
  void read(std::istream &f);
  void reset() {
    cluster_map.clear();
    c_edge_queue.clear();
    Bad_faces.clear();
    clear();
  }

private:

  // traits type
  typedef typename Geom_traits::Vector_2 Vector_2;
  typedef typename Geom_traits::Construct_translated_point_2
      Construct_translated_point_2;
  typedef typename Geom_traits::Compute_squared_distance_2
      Compute_squared_distance_2;
  typedef typename Geom_traits::Angle_2 Angle_2;
  typedef typename Geom_traits::Construct_vector_2
      Construct_vector_2;
  typedef typename Geom_traits::Construct_scaled_vector_2
      Construct_scaled_vector_2;
  typedef typename Geom_traits::Construct_midpoint_2
      Construct_midpoint_2;
  typedef typename Geom_traits::Orientation_2 Orientation_2;
  typedef typename Geom_traits::Compute_squared_minimum_sine_2 
      Compute_squared_minimum_sine_2;
  typedef typename Geom_traits::Is_bad Is_bad;

  struct Cluster {
    bool is_reduced ;

    // smallest_angle gives the two vertices defining the
    // smallest angle in the cluster
    std::pair<Vertex_handle, Vertex_handle> smallest_angle;

    FT rmin; // WARNING: rmin has no meaning if is_reduced=false!!!
    Squared_length minimum_squared_length;

    // the following map tells what segments are in the cluster and
    // they already have been splitted once
    typedef std::map<Vertex_handle, bool> Vertices_map;
    Vertices_map vertices;

    inline bool reduced() const {
      return is_reduced;
    }

    inline bool reduced(const Vertex_handle v) {
      return vertices[v];
    }
  };

  class Is_this_edge_constrained {
    Self* _m;
    Vertex_handle _v;
  public:
    Is_this_edge_constrained(Self* m, Vertex_handle v)
      : _m(m), _v(v) {}

    Is_this_edge_constrained(const Is_this_edge_constrained& other)
      : _m(other._m), _v(other._v) {}

    Is_this_edge_constrained& 
    operator=(const Is_this_edge_constrained& other)
      {
	_m = other._m;
	_v = other._v;
	return *this;
      }

    bool operator()(Vertex& v2) const 
      {
	if(_m->is_infinite(v2.handle()))
	  return false;
	Face_handle fh;
	int i;
	_m->is_edge(_v,v2.handle(),fh,i);
	return fh->is_constrained(i);
      }
  };

  typedef Filtred_circulator<Vertex_circulator,
    Is_this_edge_constrained> Constrained_vertex_circulator;

  // Bad_faces: list of bad finite faces
  // warning: some faces could be recycled during insertion in the
  //  triangulation, that's why I use a wrapper around the map
  typedef CGAL::Threetuple<Vertex_handle> Threevertices;
  class Is_really_bad {
    Self& _m;
  public:
    Is_really_bad(Self& m) : _m(m) {};
    inline
    bool operator()(const std::pair<FT, Threevertices>& p) const
      {
	const Threevertices& t = p.second;
	const Vertex_handle&
	  va = t.e0,
	  vb = t.e1,
	  vc = t.e2;
	Face_handle f;
	return( _m.is_face(va,vb,vc,f) && _m.is_bad(f));
      }
  };

  const Is_really_bad test_is_bad;
  typedef std::multimap<double, Threevertices> 
                                              Bad_faces_container_primal_type;
  CGAL::Filtred_container<Bad_faces_container_primal_type, 
    Is_really_bad> Bad_faces;

  // c_edge_queue: list of encroached constrained edges
  //  warning: some edges could be destroyed, use the same wrapper
  class Is_really_a_contrained_edge {
    const Self& _m;
  public:
    explicit Is_really_a_contrained_edge(const Self& m) : _m(m) {};
    inline
    bool operator()(const Constrained_edge& ce) const
      {
	Face_handle fh;
	int i;
	return _m.is_edge(ce.first, ce.second, fh,i) &&
	  fh->is_constrained(i);
      }
  };

  const Is_really_a_contrained_edge test_is_encroached;
  typedef std::list<Constrained_edge> List_of_constraints;
  CGAL::Filtred_container<List_of_constraints,
    Is_really_a_contrained_edge>       c_edge_queue;

  typedef std::multimap<Vertex_handle, Cluster> Cluster_map_type;
  Cluster_map_type cluster_map;
  // each vertex can have several clusters

public:
  void refine();
  void conform();
  void init();
  void mark_convex_hull();

  // It is an iterator of points
  template <class It> void init(It begin, It end);
  template <class It> void mark_facets(It begin, It end);

  void process_one_edge();
  void process_one_face();
  bool refine_step();
  //CHECK
  bool is_encroached(const Vertex_handle va, 
		     const Vertex_handle vb,
		     Point p) const;
  bool is_encroached(const Vertex_handle va, 
		     const Vertex_handle vb) const;

  inline 
  bool is_bad(const Face_handle f) const;
  
  inline
  double aspect_ratio(const Face_handle f) const;

  Mesh_2() : Tr(), test_is_bad(*this), Bad_faces(test_is_bad),
    test_is_encroached(*this), c_edge_queue(test_is_encroached) {};

  Mesh_2(List_constraints& lc, const Geom_traits& gt = Geom_traits())
    : Tr(gt), test_is_bad(*this), Bad_faces(test_is_bad),
    test_is_encroached(*this), c_edge_queue(test_is_encroached)
    {
      typename List_constraints::iterator lcit = lc.begin();
      for( ; lcit != lc.end(); ++lcit)
	{
	  insert( (*lcit).first, (*lcit).second);
	}
      CGAL_triangulation_postcondition(is_valid());
    }
  
  template <class InputIterator>
  Mesh_2(InputIterator first, InputIterator last, const Geom_traits&
       gt = Geom_traits()) : Tr(gt), test_is_bad(*this),
	 Bad_faces(test_is_bad), test_is_encroached(*this), 
	 c_edge_queue(test_is_encroached)
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
  void split_face(const Face_handle& f, const Point& circum_center);
  void create_clusters();
  void create_clusters_of_vertex(Vertex_handle v);
  void construct_cluster(Vertex_handle v,
			 Constrained_vertex_circulator begin,
			 const Constrained_vertex_circulator& end,
			 Cluster c = Cluster());
  Vertex_handle insert_middle(Face_handle f, int i);
  void cut_cluster_edge(Vertex_handle va, Vertex_handle vb, Cluster&
			c);
  Vertex_handle insert_in_c_edge(Vertex_handle va, Vertex_handle vb, Point p);
  void fill_edge_queue();
  void fill_facette_map();
  void update_c_edge_queue(Vertex_handle va,
			   Vertex_handle vb,
			   Vertex_handle vm);
  void update_facette_map(Vertex_handle v);

  // update_cluster update the cluster of [va,vb], putting vm instead
  // of vb. If reduction=false, the edge [va,vm] is not set reduced.
  void update_cluster(Cluster& c, Vertex_handle va, Vertex_handle vb,
		      Vertex_handle vm, bool reduction = true);


  inline Edge edge_between(Vertex_handle va, Vertex_handle vb);

  bool is_small_angle(const Point& pleft,
		      const Point& pmiddle, 
		      const Point& pright) const;

  // HELPING functions
  Squared_length shortest_edge_squared_length(Face_handle f);
  FT squared_cosine_of_angle_times_4(const Point& pleft,
				     const Point& pmiddle,
				     const Point& pright) const;

  // get_cluster returns the cluster of [va,vb] in c and return true
  // if it is in a cluster. If erase=true, the cluster is remove from
  // the cluster map.
  bool get_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c,
		   bool erase = false);

}; // end of Mesh_2

// # used by refine_face and cut_cluster_edge
template <class Tr>
bool Mesh_2<Tr>::
get_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c, bool erase)
{
  typedef Cluster_map_type::iterator Iterator;
  typedef std::pair<Iterator, Iterator> Range;
  Range range = cluster_map.equal_range(va);
  for(Iterator it = range.first; it != range.second; it++)
    {
      Cluster &cl = it->second;
      if(cl.vertices.find(vb)!=cl.vertices.end()) {
	c = it->second;
	if(erase)
	  cluster_map.erase(it);
	return true;
      }
    }
  return false;
}

//the function that writes a file
template <class Tr>
void Mesh_2<Tr>::
write(std::ostream &f)
{
  int nedges = 0;
  Finite_edges_iterator eit = finite_edges_begin();
  while(eit != finite_edges_end()) {
    if((*eit).first->is_constrained((*eit).second))
      nedges++;
    eit++;
  }
  f<<nedges<<std::endl;
  eit = finite_edges_begin();
  while(eit!=finite_edges_end()) {
    if((*eit).first->is_constrained((*eit).second)) {
      f<<(*eit).first->vertex(cw((*eit).second))->point().x()<<" ";
      f<<(*eit).first->vertex(cw((*eit).second))->point().y()<<" ";
      f<<(*eit).first->vertex(ccw((*eit).second))->point().x()<<" ";
      f<<(*eit).first->vertex(ccw((*eit).second))->point().y()<<std::endl;
    }
    eit++;
  }
}

//the function that reads a file
template <class Tr>
void Mesh_2<Tr>::
read(std::istream &f)
{
  int nedges = 0;
  clear();
  f>>nedges;
  for(int n = 0; n<nedges; n++) {
    FT x1, y1, x2, y2;
    f>>x1>>y1>>x2>>y2;
    insert(Point(x1, y1), Point(x2, y2));
  }
}

template <class Tr>
void Mesh_2<Tr>::
fill_edge_queue()
{
  for(Finite_edges_iterator ei = finite_edges_begin();
      ei != finite_edges_end();
      ++ei)
    {
      Vertex_handle va = (*ei).first->vertex(cw((*ei).second));
      Vertex_handle vb = (*ei).first->vertex(ccw((*ei).second));
      if((*ei).first->is_constrained((*ei).second) && 
	 is_encroached(va, vb))
	{
	  c_edge_queue.push_back(std::make_pair(va, vb));
	}
    }
}

template <class Tr>
inline
double Mesh_2<Tr>::
aspect_ratio(const Face_handle fh) const
{
  Compute_squared_minimum_sine_2 squared_sine = 
    geom_traits().compute_squared_minimum_sine_2_object();
  const Vertex_handle&
    va = fh->vertex(0),
    vb = fh->vertex(1),
    vc = fh->vertex(2);
  return squared_sine(va->point(), vb->point(), vc->point());
}
  
//it is necessarry for process_facette_map
template <class Tr>
void Mesh_2<Tr>::
fill_facette_map()
{
  for(Finite_faces_iterator fit = finite_faces_begin();
      fit != finite_faces_end();
      ++fit)
    {
      if( is_bad(fit))
	{
	  const Vertex_handle&
	    va = fit->vertex(0),
	    vb = fit->vertex(1),
	    vc = fit->vertex(2);
	  Bad_faces.insert(std::make_pair(aspect_ratio(fit),
				     Threevertices(va,vb,vc)));
	}
    }
}

//this function split all the segments that are encroached
template <class Tr>
void Mesh_2<Tr>::
refine_edge(Vertex_handle va, Vertex_handle vb)
{
  Face_handle f;
  int i;
  is_edge(va, vb, f, i); // get the edge (f,i)
  CGAL_assertion(f->is_constrained(i));
  
  Cluster c,c2;

  if( get_cluster(va,vb,c,true) )
    if( get_cluster(vb,va,c2,true) )
      { // both ends are clusters
	Vertex_handle vm = insert_middle(f,i);
	update_cluster(c,va,vb,vm,false);
	update_cluster(c2,vb,va,vm,false);
      }
    else
      // va only is a cluster
      cut_cluster_edge(va,vb,c);
  else
    if( get_cluster(vb,va,c,true) )
      // vb only is a cluster
      cut_cluster_edge(vb,va,c);
    else
      // no cluster
      insert_middle(f,i);
};

 //split all the bad faces
template <class Tr>
void Mesh_2<Tr>::
refine_face(Face_handle f)
{
  Point pc = circumcenter(f);

  typedef std::list<Edge> List_of_edges;
  typedef std::list<Face_handle> List_of_face_handles;
  List_of_edges zone_of_pc_boundary;
  List_of_face_handles zone_of_pc;

  // find conflicts around pc (starting from f as hint)
  get_conflicts_and_boundary(pc, 
			    std::back_inserter(zone_of_pc), 
			    std::back_inserter(zone_of_pc_boundary), 
			    f);
  // For the moment, we don't use the zone_of_pc.
  // It will be used when we will destroyed old bad faces in Bad_faces

  bool split_the_face = true;
  bool keep_the_face_bad = false;

  for(List_of_edges::iterator it = zone_of_pc_boundary.begin();
      it!= zone_of_pc_boundary.end();
      it++)
    {
      const Face_handle& fh = it->first;
      const int& i = it->second;
      const Vertex_handle&
	va = fh->vertex(cw(i)),
	vb = fh->vertex(ccw(i));
      if(fh->is_constrained(i) && is_encroached(va,vb,pc))
	{
	  split_the_face = false;
	  Cluster c,c2;
	  bool 
	    is_cluster_at_va = get_cluster(va,vb,c),
	    is_cluster_at_vb = get_cluster(vb,va,c2);
	  if( ( is_cluster_at_va &&  is_cluster_at_vb) || 
	      (!is_cluster_at_va && !is_cluster_at_vb) )
	    {
	      // two clusters or no cluster
	      c_edge_queue.push_back(Constrained_edge(va,vb));
	      keep_the_face_bad = true;
	    }
	  else
	    {
	      // only one cluster: c or c2
	      if(is_cluster_at_vb)
		c = c2;
// What Shewchuck says:
// - If the cluster is not reduced (all segments don't have the same
// length as [va,vb]), then split the edge
// - Else, let rmin be the minimum insertion radius introduced by the
// potential split, let T be the triangle whose circumcenter
// encroaches [va,vb] and let rg be the length of the shortest edge
// of T. If rmin >= rg, then split the edge.

	      if( !c.reduced() || 
		  c.rmin >= shortest_edge_squared_length(f) )
		{
		  c_edge_queue.push_back(Constrained_edge(va,vb));
		  keep_the_face_bad = true;
		}
	    }
	}
    }; // after here c_edge_queue contains edges encroached by pc

  const Vertex_handle&
    va = f->vertex(0),
    vb = f->vertex(1),
    vc = f->vertex(2);

  if(split_the_face)
    {
      int li;
      Locate_type lt;
      locate(pc,lt,li,f);
      if(lt!=OUTSIDE_CONVEX_HULL)
	split_face(f, pc);
    }
  else
    if(keep_the_face_bad)
      {
	Bad_faces.insert(std::make_pair(aspect_ratio(f),
				   Threevertices(va,vb,vc)));
      }
}

// # used by refine_face
template <class Tr>
inline
void Mesh_2<Tr>::
split_face(const Face_handle& f, const Point& circum_center)
{
  Vertex_handle v = insert(circum_center,f);
  update_facette_map(v);
}

template <class Tr>
void Mesh_2<Tr>::
create_clusters()
{
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr>
void Mesh_2<Tr>::
create_clusters_of_vertex(Vertex_handle v)
{
  Is_this_edge_constrained test(this, v);
  Constrained_vertex_circulator begin(incident_vertices(v),test);
  // This circulator represents all constrained edges around the
  // vertex v. An edge [v,v'] is represented by the vertex v'.

  if(begin == 0) return; // if there is only one vertex

  Constrained_vertex_circulator
    current(begin), next(begin), cluster_begin(begin);
  ++next; // next is always just after current.
  if(current == next) return;

  bool in_a_cluster = false;
  do
    {
      if(is_small_angle(current->point(), v->point(), next->point()))
	{
	  if(!in_a_cluster)
	    {
	      // at this point, current is the beginning of a cluster
	      in_a_cluster = true;
	      cluster_begin = current;
	    }
	}
      else
	if(in_a_cluster)
	  {
	    // at this point, current is the end of a cluster and
	    // cluster_begin is its beginning
 	    construct_cluster(v, cluster_begin, current);
	    in_a_cluster = false;
	  }
      ++next;
      ++current;
    } while( current!=begin );
  if(in_a_cluster)
    {
      Cluster c;
      if(get_cluster(v, begin, c, true)) 
	// get the cluster and erase it from the clusters map
	construct_cluster(v, cluster_begin, begin, c);
      else
	construct_cluster(v, cluster_begin, current);
    }
}

template <class Tr>
void Mesh_2<Tr>::
construct_cluster(Vertex_handle v,
		  Constrained_vertex_circulator begin,
		  const Constrained_vertex_circulator& end,
		  Cluster c)
{
  if(c.vertices.empty())
    {
      c.is_reduced = false;
      // c.rmin is not initialized because
      // is_reduced=false!
      c.minimum_squared_length = 
	squared_distance(v->point(), begin->point());
      Constrained_vertex_circulator second(begin);
      ++second;
      c.smallest_angle.first = begin;
      c.smallest_angle.second = second;
    }

  bool all_edges_in_cluster=false;
  if(begin==end)
    all_edges_in_cluster=true;

  Point& vp = v->point();
  
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();

  FT greatest_cosine = 
    squared_cosine_of_angle_times_4(c.smallest_angle.first->point(),
				    v->point(),
				    c.smallest_angle.second->point());

  Constrained_vertex_circulator next(begin);
  ++next;
  do
    {
      c.vertices[begin] = false;
      Squared_length l = squared_distance(vp,
					begin->point());
      c.minimum_squared_length = 
	std::min(l,c.minimum_squared_length);
      
      if(all_edges_in_cluster || begin!=end)
	{
	  FT cosine = 
	    squared_cosine_of_angle_times_4(begin->point(),
					    v->point(),
					    next->point());
	  if(cosine>greatest_cosine)
	    {
	      greatest_cosine = cosine;
	      c.smallest_angle.first = begin;
	      c.smallest_angle.second = next;
	    }
	}
    }
  while(next++,begin++!=end);
  cluster_map.insert(std::make_pair(v,c));
}

template <class Tr>
void Mesh_2<Tr>::
cut_cluster_edge(Vertex_handle va, Vertex_handle vb, Cluster& c)
{
  Construct_vector_2 vector =
    geom_traits().construct_vector_2_object();
  Construct_scaled_vector_2 scaled_vector =
    geom_traits().construct_scaled_vector_2_object();
  Compute_squared_distance_2 squared_distance =
    geom_traits().compute_squared_distance_2_object();
  Construct_midpoint_2 midpoint = 
    geom_traits().construct_midpoint_2_object();
  Construct_translated_point_2 translate =
    geom_traits().construct_translated_point_2_object();

  Vertex_handle vc;

  if(c.reduced(vb))
    {
      Face_handle fh;
      int i;
      is_edge(va,vb,fh,i);
      vc = insert_middle(fh,i);
    }
  else
    {
      const Point&
	a = va->point(),
	b = vb->point(),
	m = midpoint(a, b);

      Vector_2 v = vector(a,m);
      v = scaled_vector(v,CGAL_NTS sqrt(c.minimum_squared_length /
				      squared_distance(a,b)));
      Point 
	i = translate(a,v),
	i2(i);
	
      do {
	i = translate(a,v);
	v = scaled_vector(v,FT(2));
	i2 = translate(a,v);
      }	while(squared_distance(a,i2) <= squared_distance(a,m));
      if( squared_distance(i,m) > squared_distance(m,i2) )
	i = i2;
      //here i is the best point for splitting
      Face_handle fh;
      int index;
      is_edge(va,vb,fh,index);
      
      vc = special_insert_in_edge(i, fh, index);
    }
  update_c_edge_queue(va, vb, vc);
  update_facette_map(vc);
  update_cluster(c, va, vb, vc);
}


template <class Tr>
Mesh_2<Tr>::Vertex_handle
Mesh_2<Tr>::
insert_middle(Face_handle f, int i)
{
  Construct_midpoint_2
    midpoint = geom_traits().construct_midpoint_2_object();

  Vertex_handle
    va = f->vertex(cw(i)),
    vb = f->vertex(ccw(i));

  Point mp = midpoint(va->point(), vb->point());

  Vertex_handle vm = special_insert_in_edge(mp, f, i);
  // WARNING: special_insert_in_edge is not robust!
  // We should deconstrained the constrained edge, 
  update_c_edge_queue(va, vb, vm);
  update_facette_map(vm);
  return vm;
}

//update the encroached segments list
template <class Tr>
void Mesh_2<Tr>::
update_c_edge_queue(Vertex_handle va, Vertex_handle vb, Vertex_handle vm)
{
  Face_circulator fc = incident_faces(vm), fcbegin(fc);

  do {
    for(int i = 0; i<3; i++) {
      Vertex_handle v1 = fc->vertex(cw(i));
      Vertex_handle v2 = fc->vertex(ccw(i));
      if( fc->is_constrained(i) && !is_infinite(v1) &&
	  !is_infinite(v2) && is_encroached(v1, v2) )
	c_edge_queue.push_back(Constrained_edge(v1, v2));
    }
    ++fc;
  } while(fc != fcbegin);

  if(is_encroached(va, vm)) {
    c_edge_queue.push_back(Constrained_edge(va, vm));
  }

  if(is_encroached(vb, vm)) {
    c_edge_queue.push_back(Constrained_edge(vb, vm));
  }
}

template <class Tr>
void Mesh_2<Tr>::
update_facette_map(Vertex_handle v)
{
  Face_circulator fc = v->incident_faces(), fcbegin(fc);
  do {
    if(!is_infinite(fc)) {
      if(is_bad(fc)) {
	const Vertex_handle&
	  va = fc->vertex(0),
	  vb = fc->vertex(1),
	  vc = fc->vertex(2);
	
	Bad_faces.insert(std::make_pair(aspect_ratio(fc),
				   Threevertices(va,vb,vc)));
      }
    }
    fc++;
  } while(fc!=fcbegin);
}

template <class Tr>
void Mesh_2<Tr>::
update_cluster(Cluster& c, Vertex_handle va,Vertex_handle vb,
	       Vertex_handle vm, bool reduction)
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();
  c.vertices.erase(vb);
  c.vertices[vm] = reduction;
  
  if(vb==c.smallest_angle.first)
    c.smallest_angle.first = vm;
  if(vb==c.smallest_angle.second)
    c.smallest_angle.second = vm;

  FT l = squared_distance(va->point(),vm->point());
  if(l<c.minimum_squared_length)
    c.minimum_squared_length = l;

  if(!c.reduced())
    {
      typename Cluster::Vertices_map::iterator it = c.vertices.begin();
      while(it!=c.vertices.end() && c.reduced(it->first))
	++it; // TODO: use std::find and an object class
      if(it==c.vertices.end())
	c.is_reduced = true;
    }

  if(c.reduced())
    c.rmin = squared_distance(c.smallest_angle.first->point(),
			      c.smallest_angle.second->point())/FT(4);
  cluster_map.insert(std::make_pair(va,c));
}

template <class Tr>
bool Mesh_2<Tr>::
is_encroached(const Vertex_handle va, const Vertex_handle vb,
	      const Point p) const
{
  Angle_2 angle = geom_traits().angle_2_object();

  return angle(va->point(), p, vb->point())==OBTUSE;
}

// WARNING, TODO: NOT ALL VERTICES, ONLY THE TWO NEIGHBORS
// TODO: si les graines sont utilisées, il faut tester uniquement avec 
// les vertex dans les faces marquées.
template <class Tr>
bool Mesh_2<Tr>::
is_encroached(const Vertex_handle va, const Vertex_handle vb) const
{
  Finite_vertices_iterator vi = finite_vertices_begin();
  while(vi!=finite_vertices_end())
    {
      if(is_encroached(va, vb, vi->point()) &&
	 va!=Vertex_handle(vi) &&
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
// the measure of faces quality
// # We can add here other contraints, such as a bound on the size
template <class Tr>
inline
bool Mesh_2<Tr>::
is_bad(const Face_handle f) const
{
  const Point&
    a = f->vertex(0)->point(),
    b = f->vertex(1)->point(),
    c = f->vertex(2)->point();

  return geom_traits().is_bad_object()(a,b,c);
}


// -> traits?
template <class Tr>
bool Mesh_2<Tr>::
is_small_angle(const Point& pleft,
	       const Point& pmiddle,
	       const Point& pright) const
{
  Angle_2 angle = geom_traits().angle_2_object();
  Orientation_2 orient = geom_traits().orientation_2_object();
  
  if( angle(pleft, pmiddle, pright)==OBTUSE )
    return false;
  if( orient(pmiddle,pleft,pright)==RIGHTTURN)
    return false;

  FT cos_alpha = squared_cosine_of_angle_times_4(pleft, pmiddle,
						 pright);

  if(cos_alpha > 1)
    {
      return true; //the same cluster
    }
  else
    {
      return false; //another cluster
    }
}

// # used by: is_small_angle, create_clusters_of_vertex
// # compute 4 times the square of the cosine of the angle (ab,ac)
// # WARNING, TODO: this is not exact with doubles and can lead to crashes!!
template <class Tr>
typename Mesh_2<Tr>::FT Mesh_2<Tr>::
squared_cosine_of_angle_times_4(const Point& pb, const Point& pa,
				const Point& pc) const
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();

  const FT
    a = squared_distance(pb, pc),
    b = squared_distance(pa, pb),
    c = squared_distance(pa, pc);

  const FT num = a-(b+c);

  return (num*num)/(b*c);
}

// ->traits?
//the shortest edge that are in a triangle
// # used by: refine_face, aspect_ratio
template <class Tr>
typename Mesh_2<Tr>::FT Mesh_2<Tr>::
shortest_edge_squared_length(Face_handle f)
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();
  Point 
    pa = (f->vertex(0))->point(),
    pb = (f->vertex(1))->point(),
    pc = (f->vertex(2))->point();
  FT a, b, c;
  a = squared_distance(pb, pc);
  b = squared_distance(pc, pa);
  c = squared_distance(pa, pb);
  return (min(a, min(b, c)));
}

//the mesh refine function 
template <class Tr>
void Mesh_2<Tr>::
refine()
{
  init();
  while(! (c_edge_queue.empty() && Bad_faces.empty()) )
    {
      conform();
      if ( !Bad_faces.empty() )
	process_one_face();
    }
}

template <class Tr>
inline
void Mesh_2<Tr>::
conform()
{
  while( !c_edge_queue.empty() )
    {
      process_one_edge();	
    };
}

template <class Tr>
inline
void Mesh_2<Tr>::
init()
{
  cluster_map.clear();
  c_edge_queue.clear();
  Bad_faces.clear();
  
  create_clusters();
  fill_edge_queue();
  fill_facette_map();
  mark_convex_hull();
}

template <class Tr>
template <class It>
inline
void Mesh_2<Tr>::
init(It begin, It end)
{
  init();
  mark_facets(begin, end);
}

template <class Tr>
inline
void Mesh_2<Tr>::
process_one_edge()
{
  Constrained_edge ce = c_edge_queue.front();
  c_edge_queue.pop_front();

  refine_edge(ce.first, ce.second);
}

template <class Tr>
inline
void Mesh_2<Tr>::
process_one_face()
{
  const Threevertices& t = Bad_faces.front().second;
  const Vertex_handle&
    va = t.e0,
    vb = t.e1,
    vc = t.e2;

  Face_handle f;
  is_face(va,vb,vc,f);
  Bad_faces.pop_front();
  refine_face(f);
}

template <class Tr>
inline
bool Mesh_2<Tr>::
refine_step()
{
  if( !c_edge_queue.empty() )
    process_one_edge();
  else
    if ( !Bad_faces.empty() )
      process_one_face();
    else
      return false;
  return true;
}

template <class Tr>
void Mesh_2<Tr>::
mark_convex_hull()
{
  for(Finite_faces_iterator fit=finite_faces_begin();
      fit!=finite_faces_end();
      ++fit)
    fit->set_marked(true);
  infinite_face()->set_marked(false);
}

template <class Tr>
template <class It>
void Mesh_2<Tr>::
mark_facets(It begin, It end)
{
  for(All_faces_iterator it=all_faces_begin();
      it!=all_faces_end();
      ++it)
    it->set_marked(false);
  
  for(It it=begin; it!=end; ++it)
    {
      std::queue<Face_handle> face_queue;
      Face_handle fh=locate(*it);
      if(fh!=NULL)
      {
	face_queue.push(fh);
	fh->set_marked(true);
      }
      while( !face_queue.empty() )
	{
	  Face_handle fh = face_queue.front();
	  face_queue.pop();
	  for(int i=0;i<3;i++)
	    {
	      const Face_handle& nb = fh->neighbor(i);
	      if( !fh->is_constrained(i) && !nb->is_marked() )
		{
		  nb->set_marked(true);
		  face_queue.push(nb);
		}
	    }
	}
    }
  infinite_face()->set_marked(false);
}

CGAL_END_NAMESPACE


#endif
