#ifndef CGAL_MESH_H
#define CGAL_MESH_H
#include <CGAL/basic.h>
#include <list>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Threetuple.h>

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

template <class Circ, class Pred>
class Filtred_circulator : public Circ
{
  bool is_null;
  Pred test;

public:
  typedef Filtred_circulator<Circ,Pred> Self;


  Filtred_circulator(const Pred p=Pred()): is_null(true), test(p) {};

  Filtred_circulator(const Self& c): Circ(c), is_null(c.is_null),
    test(c.test) {};

  Self& operator=(const Self& c)
    {
      this->Circ::operator=(c);
      is_null=c.is_null;
      test=c.test;
      return *this;
    }

  Filtred_circulator(const Circ& c, const Pred p=Pred())
    : Circ(c), is_null(false), test(p)
    {
      if(test(**this))
	is_null=false;
      else
	{
	  Self end(*this);
	  do { 
	    this->Circ::operator++();
	  } while( !test(**this) && (*this)!=end );
	  if((*this)==end)
	    is_null=true;
	}
    };

  bool operator==( CGAL_NULL_TYPE p) const {
    return is_null;
  }

  bool operator!=( CGAL_NULL_TYPE p) const {
    return !is_null;
  }

  bool operator==(const Self& c) const 
    {
      return is_null==c.is_null && this->Circ::operator==(c);
    }

  bool operator!=( const Self& c) const { return !(*this == c); }

  Self& operator++() {
    CGAL_assertion(!is_null);
    do {
      this->Circ::operator++();
    } while( !test(**this) );
    return *this;
  }

  Self  operator++(int) {
    Self tmp= *this;
    ++*this;
    return tmp;
  }
};

template <class Cont, class Pred>
class Filtred_container : public Cont
{
  Pred test;
public:
  Filtred_container(Pred p) : Cont(), test(p) {};
  Filtred_container(Cont& c, Pred p=Pred()) : Cont(c), test(p) {};

  typedef typename Cont::reference reference;
  typedef typename Cont::iterator iterator;

  inline
  reference front()
    {
      iterator r=begin();
      while(!test(*r))
	{
	  erase(r);
	  r=begin();
	}
      return *r;
    }

  inline
  bool empty()
    {
      if(Cont::empty())
	return true;
      else
	{
	  while(!Cont::empty())
	    {
	      iterator r=begin();
	    if(!test(*r))
	      pop_front();
	    else
	      return false;
	    }
	  return true;
	}
    }

  inline
  void pop_front()
    {
      erase(begin());
    }
};

// Tr is a Delaunay constrained triangulation
// Mtraits is a mesh trait
template <class Tr, class Mtraits = void>
class Mesh_2: public Tr
{
public:
  typedef Tr Triangulation;
  typedef Mtraits  Mesh_2_traits;
  typedef Mesh_2<Triangulation, Mesh_2_traits> Self;
  
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;
  typedef typename Geom_traits::FT FT;
  //  typedef FT      Length;
  typedef FT      Squared_length;

  typedef typename Tr::Vertex                 Vertex;
  typedef typename Tr::Edge                   Edge;
  //  typedef typename Tr::Edge_iterator          Edge_iterator;
  typedef typename Tr::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Tr::Edge_circulator        Edge_circulator;
  typedef typename Tr::Face_handle            Face_handle;
//   typedef typename Tr::Face_iterator          Face_iterator;
  typedef typename Tr::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Tr::Face_circulator        Face_circulator;
  typedef typename Tr::Vertex_handle          Vertex_handle;
  //  typedef typename Tr::Vertex_iterator        Vertex_iterator;
  typedef typename Tr::Vertex_circulator      Vertex_circulator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Tr::Locate_type            Locate_type;
  
  typedef std::pair<Vertex_handle,Vertex_handle>
                                              Constrained_edge;

  typedef typename Tr::Point                  Point;

//   typedef typename Tr::Constraint             Constraint;
  typedef typename Tr::List_constraints       List_constraints;

  void write(ostream &f);
  void read(istream &f);
  void reset() {
    cluster_map.clear();
    clear();
  }

private:
  typedef typename Geom_traits::Compute_squared_distance_2
      Compute_squared_distance_2;
  typedef typename Geom_traits::Angle_2 Angle_2;
    
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
	_m=other._m;
	_v=other._v;
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
	const Threevertices& t=p.second;
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
  class Is_really_an_encroached_edge {
    const Self& _m;
  public:
    explicit Is_really_an_encroached_edge(const Self& m) : _m(m) {};
    inline
    bool operator()(const Constrained_edge& ce) const
      {
	Face_handle fh;
	int i;
	return _m.is_edge( ce.first, ce.second, fh,i) &&
	  fh->is_constrained(i) && _m.is_encroached(ce.first, ce.second);
      }
  };

  const Is_really_an_encroached_edge test_is_encroached;
  typedef std::list<Constrained_edge> List_of_constraints;
  CGAL::Filtred_container<List_of_constraints,
    Is_really_an_encroached_edge>       c_edge_queue;

  std::multimap<Vertex_handle, Cluster>  cluster_map;
  // each vertex can have several clusters

public:
  //INSERTION-REMOVAL
  void refine_mesh();
  // TODO: refine_mesh_step(), that do a step of refinment

  //CHECK
  bool is_encroached(const Vertex_handle va, 
		     const Vertex_handle vb,
		     Point p) const;
  bool is_encroached(const Vertex_handle va, 
		     const Vertex_handle vb) const;

  bool is_bad(const Face_handle f);
  FT aspect_ratio(Face_handle f); // r^2/l^2

  Mesh_2() : Tr(), test_is_bad(*this), Bad_faces(test_is_bad),
    test_is_encroached(*this), c_edge_queue(test_is_encroached) {};

  Mesh_2(List_constraints& lc, const Geom_traits& gt=Geom_traits())
    : Tr(gt), test_is_bad(*this), Bad_faces(test_is_bad),
    test_is_encroached(*this), c_edge_queue(test_is_encroached)
    {
      typename List_constraints::iterator lcit=lc.begin();
      for( ; lcit != lc.end(); ++lcit)
	{
	  insert( (*lcit).first, (*lcit).second);
	}
      CGAL_triangulation_postcondition(is_valid());
    }
  
  template <class InputIterator>
  Mesh_2(InputIterator first, InputIterator last, const Geom_traits&
       gt=Geom_traits()) : Tr(gt), test_is_bad(*this),
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
  Vertex_handle insert_middle(Face_handle f, int i);
  void cut_cluster_edge(Vertex_handle va, Vertex_handle vb, Cluster&
			c);
  // TODO: rewrite cut_cluster_edge
  Vertex_handle insert_in_c_edge(Vertex_handle va, Vertex_handle vb, Point p);
  void fill_edge_queue();
  void fill_facette_map();
  void update_c_edge_queue(Vertex_handle va,
			   Vertex_handle vb,
			   Vertex_handle vm);
  void update_facette_map(Vertex_handle v);
  void update_cluster(Cluster& c, Vertex_handle va, Vertex_handle vb,
		      Vertex_handle vm);
  inline Edge edge_between(Vertex_handle va, Vertex_handle vb);

  bool is_small_angle(Vertex_handle vleft, 
		      Vertex_handle vmiddle, 
		      Vertex_handle vright);

  // HELPING functions
  Squared_length shortest_edge_squared_length(Face_handle f);
  FT squared_cosine_of_angle_times_4(Vertex_handle va,
					Vertex_handle vb,
					Vertex_handle vc);
  Vertex_handle nearest_incident_vertex(Vertex_handle v);
  bool get_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c);

}; // end of Mesh_2

// # used by refine_face and cut_cluster_edge
template <class Tr, class Mtraits>
bool Mesh_2<Tr, Mtraits>::
get_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c)
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
void Mesh_2<Tr, Mtraits>::
write(ostream &f)
{
  int nedges = 0;
  Finite_edges_iterator eit = finite_edges_begin();
  while(eit != finite_edges_end()) {
    if((*eit).first->is_constrained((*eit).second))
      nedges++;
    eit++;
  }
  f<<nedges<<endl;
  eit=finite_edges_begin();
  while(eit!=finite_edges_end()) {
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
void Mesh_2<Tr, Mtraits>::
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
void Mesh_2<Tr, Mtraits>::
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
	  c_edge_queue.push_back(make_pair(va, vb));
	}
    }
}

//it is necessarry for process_facette_map
template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
fill_facette_map()
{
  for(Finite_faces_iterator fit = finite_faces_begin();
      fit != finite_faces_end();
      ++fit)
    {
      if( is_bad(fit))
	{
	  const Vertex_handle&
	    va=fit->vertex(0),
	    vb=fit->vertex(1),
	    vc=fit->vertex(2);
	  Bad_faces.insert(make_pair(aspect_ratio(fit),
				     Threevertices(va,vb,vc)));
	}
    }
}

//this function split all the segments that are encroached
template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
refine_edge(Vertex_handle va, Vertex_handle vb)
{
  Face_handle f;
  int i;
  is_edge(va, vb, f, i); // get the edge (f,i)
  CGAL_assertion(f->is_constrained(i));
  
  Cluster c,c2;

  if( get_cluster(va,vb,c) )
    if( get_cluster(vb,va,c2) )
      { // both ends are clusters
	Vertex_handle vm = insert_middle(f,i);
	update_cluster(c,va,vb,vm);
	update_cluster(c2,vb,va,vm);
      }
    else
      // va only is a cluster
      cut_cluster_edge(va,vb,c);
  else
    if( get_cluster(vb,va,c) )
      // vb only is a cluster
      cut_cluster_edge(vb,va,c);
    else
      // no cluster
      insert_middle(f,i);
};

 //split all the bad faces
template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
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

  for(List_of_edges::iterator it=zone_of_pc_boundary.begin();
      it!=zone_of_pc_boundary.end();
      it++)
    {
      const Vertex_handle&
	va=it->first->vertex(cw(it->second)),
	vb=it->first->vertex(ccw(it->second));
      if(is_encroached(va,vb,pc))
	{
	  Cluster c;
	  if(get_cluster(va,vb,c))
	    {
	      if( !c.is_reduced || 
		  c.rmin >= shortest_edge_squared_length(f) )
		c_edge_queue.push_back(Constrained_edge(va,vb));
	    }
	  else
	    c_edge_queue.push_back(Constrained_edge(va,vb));
	}
    }; // after here c_edge_queue contains edges encroached by pc

  if(c_edge_queue.empty())
    {
      int li;
      Locate_type lt;
      locate(pc,lt,li,f);
      if(lt!=OUTSIDE_CONVEX_HULL)
	split_face(f, pc);
    }
}

// # used by refine_face
template <class Tr, class Mtraits>
inline
void Mesh_2<Tr, Mtraits>::
split_face(const Face_handle& f, const Point& circum_center)
{
  Vertex_handle v = insert(circum_center,f);
  update_facette_map(v);
}

template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
create_clusters()
{
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
create_clusters_of_vertex(Vertex_handle v)
{
  Is_this_edge_constrained test(this, v);
  Constrained_vertex_circulator begin(incident_vertices(v),test);
  // This circulator represents all constrained edges around the
  // vertex v. An edge [v,v'] is represented by the vertex v'.

  if(begin == 0) return; // if there is only one vertex

  Constrained_vertex_circulator
    current(begin), next(begin), cluster_edge(begin);
  ++next; // next is always just after current.

  bool in_a_cluster=false;
  while(next!=begin)
    {
      if(is_small_angle(current, v, next))
	{
	  if(!in_a_cluster)
	    { 
	      // at this point, current is the beginning of a cluster
	      in_a_cluster=true;
	      cluster_edge=current;
	    }
	}
      else
	if(in_a_cluster)
	  {
	    // at this point, current is the end of a cluster and
	    // cluster_edge is its beginning
	    in_a_cluster=false;
	    Cluster new_cluster;
	    Point& vp = v->point();

	    Compute_squared_distance_2 squared_length = 
	      geom_traits().compute_squared_distance_2_object();

	    new_cluster.is_reduced=false;
	    // new_cluster.rmin is not initialized because is_reduced=false!
	    new_cluster.minimum_squared_length = 
	      squared_length(vp, cluster_edge->point());

	    Constrained_vertex_circulator
	      next_cluster_edge(cluster_edge);
	    ++next_cluster_edge;
	    FT greatest_cosine = 
	      squared_cosine_of_angle_times_4(v,
					      cluster_edge,
					      next_cluster_edge);
	    new_cluster.smallest_angle.first=cluster_edge;
	    new_cluster.smallest_angle.second=next_cluster_edge;
	    do
	      {
		new_cluster.vertices[cluster_edge]=false;
		Squared_length l=squared_length(vp,
					       cluster_edge->point());
		new_cluster.minimum_squared_length=
		  std::min(l,new_cluster.minimum_squared_length);
		new_cluster.is_reduced=false;

		if(cluster_edge!=next)
		  {
		    FT
		      cosine=squared_cosine_of_angle_times_4(v, cluster_edge, next_cluster_edge);
		    if(l>greatest_cosine)
		      {
			greatest_cosine=cosine;
			new_cluster.smallest_angle.first=cluster_edge;
			new_cluster.smallest_angle.second=next_cluster_edge;
		      }
		  }
		++next_cluster_edge;
		++cluster_edge;
	      }
	    while(cluster_edge!=next);
	  }
      ++next;
      ++current;
    }
}

//refine the cluster edges
// # used by: refine_edge
template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
cut_cluster_edge(Vertex_handle va, Vertex_handle vb, Cluster& c)
{
  // What Shewchuck says:
  // - If the cluster is not reduced (all segments don't have the same
  // length as [va,vb]), then split the edge (at midpoint)
  // - Else, let rmin be the minimum insertion radius introduced by the
  // potential split, let T be the triangle whose circumcenter
  // encroaches [va,vb] and let rg be the length of the shortest edge
  // of T. If rmin >= rg, then split the edge.

  Squared_length l2 = 2.0; //shortest_edge_of_cluster(va, c)/1.001;
  Squared_length L2 = squared_distance(va->point(), vb->point());

  // # WARNING, TODO: what will append if a edge already has a power
  // of two length??

  if(L2 > l2)
    {
      FT rapport = pow(2.0, rint(0.5*log(L2/l2)/log(2.0)))*::sqrt(l2/(L2*4));
      Point pc=va->point()+(vb->point()-va->point())*rapport;

      Face_handle fh;
      int i;
      is_edge(va,vb,fh,i);

      Vertex_handle vc = special_insert_in_edge(pc, fh, i);
      update_c_edge_queue(va, vb, vc);
      update_facette_map(vc);
      update_cluster(c, va, vb, vc);
    }
}


template <class Tr, class Mtraits>
Mesh_2<Tr, Mtraits>::Vertex_handle
Mesh_2<Tr, Mtraits>::
insert_middle(Face_handle f, int i)
{
  Vertex_handle
    va=f->vertex(cw(i)),
    vb=f->vertex(ccw(i));

  Point mp = midpoint(va->point(), vb->point());

  Vertex_handle vm = special_insert_in_edge(mp, f, i);
  // WARNING: special_insert_in_edge is not robust!
  // We should deconstrained the constrained edge, 
  update_c_edge_queue(va, vb, vm);
  update_facette_map(vm);
  return vm;
}

//update the encroached segments list
// # TODO: rewrite this!!
template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
update_c_edge_queue(Vertex_handle va, Vertex_handle vb, Vertex_handle vm)
{
  Face_circulator fc = incident_faces(vm), fcbegin(fc);

  do {
    for(int i=0; i<3; i++) {
      Vertex_handle v1=fc->vertex(cw(i));
      Vertex_handle v2=fc->vertex(ccw(i));
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

template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
update_facette_map(Vertex_handle v)
{
  Face_circulator fc = v->incident_faces(), fcbegin(fc);
  do {
    if(!is_infinite(fc)) {
      if(is_bad(fc)) {
	const Vertex_handle&
	  va=fc->vertex(0),
	  vb=fc->vertex(1),
	  vc=fc->vertex(2);
	
	Bad_faces.insert(make_pair(aspect_ratio(fc),
				   Threevertices(va,vb,vc)));
      }
    }
    fc++;
  } while(fc!=fcbegin);
}

template <class Tr, class Mtraits>
void Mesh_2<Tr, Mtraits>::
update_cluster(Cluster& c, Vertex_handle va,Vertex_handle vb,
	       Vertex_handle vm)
{
  Compute_squared_distance_2 squared_length = 
    geom_traits().compute_squared_distance_2_object();
  c.vertices.erase(vb);
  c.vertices[vm]=true;
  
  if(vb==c.smallest_angle.first)
    c.smallest_angle.first=vm;
  if(vb==c.smallest_angle.second)
    c.smallest_angle.second=vm;

  FT l=squared_length(va->point(),vm->point());
  if(l<c.minimum_squared_length)
    c.minimum_squared_length=l;

  if(c.is_reduced)
    c.is_reduced=false; // c *was* reduced!
  else
    {
      typename Cluster::Vertices_map::iterator it=c.vertices.begin();
      while(it!=c.vertices.end() &&
	    squared_length(va->point(),it->first->point()) ==
	    c.minimum_squared_length)
	++it; // TODO: use std::find and an object class
      if(it==c.vertices.end())
	c.is_reduced=true;
      else
	c.is_reduced=false;
    }

  if(c.is_reduced)
    c.rmin=squared_length(c.smallest_angle.first->point(),
			c.smallest_angle.second->point())/FT(2);
}

template <class Tr, class Mtraits>
bool Mesh_2<Tr, Mtraits>::
is_encroached(const Vertex_handle va, const Vertex_handle vb,
	      const Point p) const
{
  Angle_2 angle=geom_traits().angle_2_object();

  return angle(va->point(), p, vb->point())==OBTUSE;
}

// NOT ALL VERTICES, ONLY THE TWO NEIGHBORS
template <class Tr, class Mtraits>
bool Mesh_2<Tr, Mtraits>::
is_encroached(const Vertex_handle va, const Vertex_handle vb) const
{
  Finite_vertices_iterator vi=finite_vertices_begin();
  while(vi!=finite_vertices_end())
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
// the measure of faces quality
// # We can add here other contraints, such as a bound on the size
template <class Tr, class Mtraits>
bool Mesh_2<Tr, Mtraits>::
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
bool Mesh_2<Tr, Mtraits>::
is_small_angle(Vertex_handle vleft,
	       Vertex_handle vmiddle,
	       Vertex_handle vright)
{
  FT cos_alpha = squared_cosine_of_angle_times_4(vmiddle, vleft, vright);
  if(cos_alpha > 1)
    {
      return true; //the same cluster
    }
  else
    {
      return false; //another cluster
    }
}

// template <class Tr, class Mtraits>
// typename Mesh_2<Tr, Mtraits>::FT Mesh_2<Tr, Mtraits>::
// shortest_edge_of_cluster(Vertex_handle v, Cluster &cluster)
// { 
//   map<Vertex_handle, Length>::iterator vit = cluster.vertices.begin();
//   FT min_edge = squared_distance(((*vit).first)->point(), (*v).point());
//   vit++;
//   while(vit != cluster.vertices.end())
//     {
//       FT temp = squared_distance(((*vit).first)->point(), (*v).point());
//       if(temp < min_edge)
// 	{
// 	  min_edge = temp;
// 	}
//      vit++;
//     }
//   return min_edge;
// }

// look if all edges have the same length
// # A reduced cluster is a cluster wherein all edges have the same
// length
// template <class Tr, class Mtraits>
// void Mesh_2<Tr, Mtraits>:: 
// check_cluster_status( Cluster& cluster)
// {
//   Length initl;
//   map<Vertex_handle, Length>::iterator vi = cluster.vertices.begin();
//   initl = (*(cluster.vertices.find((*vi).first))).second;
//   vi++;
//   for(;vi != cluster.vertices.end(); vi++)
//     {
//       if(( (*(cluster.vertices.find((*vi).first))).second)  != initl)
// 	{
// 	  cluster.status = NON_REDUCED;

// 	  return;
// 	}
//       vi++;  
//     }
//   cluster.status = REDUCED;
// }

// # used by: is_small_angle
// compute 4 times the square of the cosine of the angle (ab,ac)
template <class Tr, class Mtraits>
typename Mesh_2<Tr, Mtraits>::FT Mesh_2<Tr, Mtraits>::
squared_cosine_of_angle_times_4(Vertex_handle va, Vertex_handle vb,
				   Vertex_handle vc)
{
  Compute_squared_distance_2 squared_length = 
    geom_traits().compute_squared_distance_2_object();  

  Point 
    pa = va->point(),
    pb = vb->point(),
    pc = vc->point();
  FT
    a = squared_distance(pb, pc),
    b = squared_distance(pc, pa),
    c = squared_distance(pa, pb);
  return (a-(b+c))/(b*c);
}

// ->traits?
//the shortest edge that are in a triangle
// # used by: refine_face, aspect_ratio
template <class Tr, class Mtraits>
typename Mesh_2<Tr, Mtraits>::FT Mesh_2<Tr, Mtraits>::
shortest_edge_squared_length(Face_handle f)
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();
  Point 
    pa = (f->vertex(0))->point(),
    pb = (f->vertex(1))->point(),
    pc = (f->vertex(2))->point();
  FT a, b, c;
  a =squared_distance(pb, pc);
  b =squared_distance(pc, pa);
  c =squared_distance(pa, pb);
  return (min(a, min(b, c)));
}

// ->traits
//the triangle quality is represented by the
//aspect_ratio value
// # used by: fill_facette_map, update_facette_map, is_bad
template <class Tr, class Mtraits>
typename Mesh_2<Tr, Mtraits>::FT Mesh_2<Tr, Mtraits>::
aspect_ratio(Face_handle f)
{
  Point p;
  p = circumcenter(f);
  Point A = (f->vertex(0))->point();
  FT radius = squared_distance(p, A);
  FT sh_edge = shortest_edge_squared_length(f);
  return (radius/sh_edge);
}

//this function must compute the vertex that are so close to our vertex.
template <class Tr, class Mtraits>
Mesh_2<Tr, Mtraits>::Vertex_handle Mesh_2<Tr, Mtraits>::
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
void Mesh_2<Tr, Mtraits>::
refine_mesh()
{
  create_clusters();

  fill_edge_queue();
  fill_facette_map();
  while(! (c_edge_queue.empty() && Bad_faces.empty()) )
    {
      while( !c_edge_queue.empty() )
	{
	  Constrained_edge ce=c_edge_queue.front();
	  c_edge_queue.pop_front();

	  std::cerr << "Encroached edge: "
		    << ce.first->point() << " , "
		    << ce.second->point() << std::endl;

	  refine_edge(ce.first, ce.second);
	};
      while( !Bad_faces.empty() )
	{
	  const Threevertices& t = Bad_faces.front().second;
	  const Vertex_handle&
	    va = t.e0,
	    vb = t.e1,
	    vc = t.e2;

	  std::cerr << "Bad face: "
		    << va->point() << " , "
		    << vb->point() << " , "
		    << vc->point() << std::endl;

	  Face_handle f;
	  is_face(va,vb,vc,f);
	  Bad_faces.pop_front();
	  refine_face(f);
	}
    }
}

CGAL_END_NAMESPACE


#endif

