#ifndef CGAL_CONFORM_2_H
#define CGAL_CONFORM_2_H
#include <CGAL/Conform_2.h>
#include <list>
#include <map>
#include <iterator>
#include <functional>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/utility.h>
#include <CGAL/iterator.h>
#include <CGAL/Filtered_container.h>
#include <CGAL/Filter_circulator.h>

#ifdef CGAL_USE_BOOST
#  include <boost/iterator_adaptors.hpp>
#else
#  include <CGAL/Iterator_project.h>
#endif

CGAL_BEGIN_NAMESPACE

/**
   - Tr is a Delaunay constrained triangulation (with intersections or
   not).
*/
template <class Tr>
class Conform_triangulation_2: public Tr
{
public:
  // -- public typedef --
  typedef Tr Triangulation;
  typedef Conform_triangulation_2<Tr> Conform; // for use by nested types
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  typedef typename Geom_traits::FT FT;
  typedef FT      Squared_length;

  // -- types inherited from the templated base class --
  // must be redefined for use in the implementation
  typedef typename Tr::Vertex                 Vertex;
  typedef typename Tr::Vertex_handle          Vertex_handle;
  typedef typename Tr::Vertex_circulator      Vertex_circulator;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Tr::Finite_edges_iterator  Finite_edges_iterator;

  typedef typename Tr::Face_handle            Face_handle;
  typedef typename Tr::Face_circulator        Face_circulator;
  typedef typename Tr::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Tr::All_faces_iterator     All_faces_iterator;

  typedef typename Tr::Locate_type            Locate_type;
  
  typedef typename Tr::Point                  Point;

protected:
  // typedefs for private members types
  typedef std::pair<Vertex_handle,Vertex_handle> Constrained_edge;

  /**
     Is_edge_constrained: function object class. 
     Take a Conform_triangulation_2 object m and a Vertex_handle v in
     constructor. Is_edge_constrained(Vertex_handle v2) tells if
     [v,v2] is a constrained edge in m. */
  class Is_edge_constrained {
    Conform* _m;
    Vertex_handle _v;
  public:
    Is_edge_constrained(Conform* m, Vertex_handle v)
      : _m(m), _v(v) {}

    bool operator()(/*const*/ Vertex& v2) const 
      {
	if(_m->is_infinite(v2.handle()))
	  return false;
	Face_handle fh;
	int i;
	_m->is_edge(_v,v2.handle(),fh,i);
	return fh->is_constrained(i);
      }
  };

  class Is_really_a_constrained_edge {
    const Conform& _m;
  public:
    explicit Is_really_a_constrained_edge(const Conform& m) : _m(m) {};
    bool operator()(const Constrained_edge& ce) const
      {
	Face_handle fh;
	int i;
	return _m.is_edge(ce.first, ce.second, fh,i) &&
	  fh->is_constrained(i);
      }
  };

  typedef Filtred_circulator<Vertex_circulator, Is_edge_constrained>
    Constrained_vertex_circulator;

  typedef std::list<Constrained_edge> List_of_constraints;
  typedef CGAL::Filtered_container<List_of_constraints, 
                                  Is_really_a_constrained_edge>
    Constrained_edges_queue;

  // Cluster register several informations about clusters.
  // A cluster is a is a set of vertices v_i incident to one vertice
  // v_0, so that angles between segments [v_0, v_i] is less than 60°.
  struct Cluster {
    bool reduced ; // Is the cluster reduced

    // smallest_angle gives the two vertices defining the
    // smallest angle in the cluster
    std::pair<Vertex_handle, Vertex_handle> smallest_angle;

    FT rmin; // WARNING: rmin has no meaning if reduced=false!!!
    Squared_length minimum_squared_length;

    // The following map tells what vertices are in the cluster and if 
    // the corresponding segment has been splitted once.
    typedef std::map<Vertex_handle, bool> Vertices_map;
    Vertices_map vertices;

    bool is_reduced() const {
      return reduced;
    }

    bool is_reduced(const Vertex_handle v) {
      return vertices[v];
    }
  };

  typedef std::multimap<Vertex_handle, Cluster> Cluster_map;

private:
  template <class Pair>
  struct Pair_get_first: public std::unary_function<Pair,
						    typename Pair::first_type>
  {
    typedef typename Pair::first_type result;
    const result& operator()(const Pair& p) const
    {
      return p.first;
    }
  };

  typedef typename Cluster::Vertices_map Cluster_vertices_map;

#ifdef CGAL_USE_BOOST
public:
  typedef typename boost::projection_iterator_generator<
    Pair_get_first<typename Cluster_map::value_type>,
    typename Cluster_map::const_iterator>::type
  Cluster_vertices_iterator;

  typedef typename boost::projection_iterator_generator<
    Pair_get_first<typename Cluster_vertices_map::value_type>,
    typename Cluster_vertices_map::const_iterator>::type
  Vertices_in_cluster_iterator;
#else
  typedef CGAL::Iterator_project<typename Cluster_map::const_iterator,
    Pair_get_first<typename Cluster_map::value_type> >
  Cluster_vertices_iterator;

  typedef CGAL::Iterator_project<
    typename Cluster_vertices_map::const_iterator,
    Pair_get_first<typename Cluster_vertices_map::value_type> >
  Vertices_in_cluster_iterator;
#endif // CGAL_USE_BOOST

  // -- conform criteria --
public:
  struct Is_locally_gabriel_conform
  {
    bool operator()(const Conform& ct,
		    const Face_handle& fh,
		    const int i) const
      {
	typedef typename Geom_traits::Angle_2 Angle_2;
	
	const Angle_2 angle = ct.geom_traits().angle_2_object();
	
	const Point& a = fh->vertex(ct.cw(i))->point();
	const Point& b = fh->vertex(ct.ccw(i))->point();
	const Point& c = fh->vertex(i)->point();
	const Point& d = fh->mirror_vertex(i)->point();
	
	return( angle(a, c, b) != OBTUSE &&
		angle(a, d, b) != OBTUSE );
      }
    bool operator()(const Conform& ct,
		    const Face_handle& fh,
		    const int i,
		    const Point& p) const
      {
	typedef typename Geom_traits::Angle_2 Angle_2;
	
	const Angle_2 angle = ct.geom_traits().angle_2_object();
	
	const Point& a = fh->vertex(ct.cw(i))->point();
	const Point& b = fh->vertex(ct.ccw(i))->point();
	
	return( angle(a, p, b) != OBTUSE );
      }

  };
  
  struct Is_locally_delaunay_conform
  {
    bool operator()(const Conform& ct,
		    const Face_handle& fh,
		    const int i) const
      {
	typedef typename Geom_traits::Side_of_oriented_circle_2
	  Side_of_oriented_circle_2;
	
	const Side_of_oriented_circle_2 in_circle =
	  ct.geom_traits.side_of_oriented_circle_2_object();
	
	const Point& a = fh->vertex(ct.cw(i))->point();
	const Point& b = fh->vertex(ct.ccw(i))->point();
	const Point& c = fh->vertex(i)->point();
	const Point& d = fh->mirror_vertex(i)->point();
	
	return( in_circle(c, b, a, d) == ON_POSITIVE_SIDE );
      }
  };

public:
  // --- CONSTRUCTORS ---

  /** default constructor */
  explicit
  Conform_triangulation_2(const Geom_traits& gt = Geom_traits());

  // --- SURCHARGED INSERTION-DELETION FONCTIONS ---
  // TODO!

  // --- ASSIGNEMENT ---
  // TODO!

  // --- ACCESS FUNCTIONS ---
  int number_of_constrained_edges() const;

  int number_of_clusters_vertices() const
    {
      return cluster_map.size();
    }

  Cluster_vertices_iterator clusters_vertices_begin() const
  {
    return Cluster_vertices_iterator(cluster_map.begin());
  }

  Cluster_vertices_iterator clusters_vertices_end() const
  {
    return Cluster_vertices_iterator(cluster_map.end());
  }

  unsigned int number_of_clusters_at_vertex(Vertex_handle vh)
  {
    typedef typename Cluster_map::iterator Iterator;
    typedef std::pair<Iterator, Iterator> Range;
    Range range = cluster_map.equal_range(vh);
    return std::distance(range.first, range.second);
  }

  // returns the sequence of vertices bellonging to the n-th cluster of vh
  std::pair<Vertices_in_cluster_iterator, Vertices_in_cluster_iterator>
  vertices_in_cluster_sequence(Vertex_handle vh, unsigned int n)
  {
    typedef typename Cluster_map::iterator Iterator;
    typedef std::pair<Iterator, Iterator> Range;
    typedef typename Range::first_type Clusters_iterator;
    typedef Pair_get_first<typename Cluster_vertices_map::value_type>
      Get_first;

    Range range = cluster_map.equal_range(vh);
    Iterator first = range.first;
    std::advance(first, n);
    const Cluster& c = first->second;

    return
      std::make_pair(Vertices_in_cluster_iterator(c.vertices.begin()),
		     Vertices_in_cluster_iterator(c.vertices.end()));
  }

  // --- HELPING FUNCTION ---
  void clear();

  // --- CONFORMING FUNCTIONS ---

  // Conform edges
  void gabriel_conform();
  void delaunay_conform();

  // -- other conform policies --
  // TODO: should be protected and not public
  template <class Is_locally_conform>
  void conform(const Is_locally_conform& is_loc_conf)
  {
    while( !edges_to_be_conformed.empty() )
    {
      process_one_edge(is_loc_conf);
    };
  }

  // --- STEP BY STEP FUNCTIONS ---

  /**
     init(...): Initialize the data structures 
     (The call of this function is REQUIRED before any step by step
     operation).
  */
  template <class Is_locally_conform>
  void init(const Is_locally_conform&);

  /** Execute on step of the algorithm.
      init() should have been called before.
  */
  template <class Is_locally_conform>
  bool refine_step(const Is_locally_conform&);

  /** Tells if all constrained edges are conformed. */
  bool is_conformed()
    // This function cannot be "const" because, as edges_to_be_conformed is
    // filtred, its empty() method is not const.
  { return edges_to_be_conformed.empty(); }

  // --- PROTECTED TYPES ---
protected:
  typedef std::list<Face_handle> List_of_face_handles;

  // -- traits type --
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

private:
  // PRIVATE MEMBER DATAS

  // edges_to_be_conformed: list of encroached constrained edges
  //  warning: some edges could be destroyed, use the same wrapper
  // is_really_a_constrained_edge: tester to filter edges_to_be_conformed
  const Is_really_a_constrained_edge is_really_a_constrained_edge;
  Constrained_edges_queue edges_to_be_conformed;

  // Cluster_map: multimap Vertex_handle -> Cluster
  // each vertex can have several clusters
  Cluster_map cluster_map;

  // tell if init has been called since the last insertion of a
  // constraint
  bool initialized;

  // --- PRIVATE MEMBER FUNCTIONS ---
private: 

  // -- auxiliary functions to handle clusters --

  // for all vertices, call create_clusters_of_vertex
  void create_clusters();

  // compute clusters of the vertex v, using the auxiliary function
  // construct_cluster
  void create_clusters_of_vertex(const Vertex_handle v);

  // add the sequence [begin, end] to the cluster c and add it to the
  // clusters of the vertex v
  void construct_cluster(const Vertex_handle v,
			 Constrained_vertex_circulator begin,
			 const Constrained_vertex_circulator& end,
			 Cluster c = Cluster());

  // see more functions about clusters in protected member functions


  // -- functions that maintain the queue of encroached edges --

  // scan all constrained edges and put them in the queue if they are
  // encroached
  template <class Is_locally_conform>
  void fill_edge_queue(const Is_locally_conform&);

  // update the queue with edges incident to vm
  template <class Is_locally_conform>
  void update_edges_to_be_conformed(const Vertex_handle va,
				    const Vertex_handle vb,
				    const Vertex_handle vm,
				    const Is_locally_conform&);

  // -- inlined functions that compose the refinement process --

  // take one edge in the queue and call refine_edge
  template <class Is_locally_conform>
  void process_one_edge(const Is_locally_conform&);

  // handle the encroached edge, call cut_cluster_edge+update_cluster
  // or insert_middle
  template <class Is_locally_conform>
  void refine_edge(const Vertex_handle va, const Vertex_handle vb,
		   const Is_locally_conform&);

  // -- auxiliary functions that return a boolean --

  // tell if [va,vb] is encroached, by looking for the two neighbors
  // This function takes care of markers.
//   bool is_encroached(const Vertex_handle va, 
// 		     const Vertex_handle vb) const;

//   // tell if [va,vb] is encroached by p
//   bool is_encroached(const Vertex_handle va, 
// 		     const Vertex_handle vb,
// 		     const Point& p) const;

  // tell if the angle <pleft, pmiddle, pright> is less than 60°
  // Uses squared_cosine_of_angle_times_4 and used by
  // create_clusters_of_vertex
  bool is_small_angle(const Point& pleft,
		      const Point& pmiddle, 
		      const Point& pright) const;

  // -- auxiliary functions that are called to split an edge or --
  // -- a face                                                  --

  // cut [va,vb] knowing that it is in the cluster c
  Vertex_handle cut_cluster_edge(const Vertex_handle va,
				 const Vertex_handle vb,
				 Cluster& c);

  // insert the midpoint of the edge (f,i)
  Vertex_handle insert_middle(Face_handle f, const int i);


  // -- functions that really insert points --

  // virtual function that inserts the point p in the edge
  // (fh,edge_index)
  virtual
  Vertex_handle virtual_insert_in_the_edge(Face_handle fh,
				   const int edge_index,
				   const Point& p);
  
  // -- helping computing functions -- 

  // return the squared cosine of the angle <pleft, pmiddle, pright>
  // times 4
  FT squared_cosine_of_angle_times_4(const Point& pleft,
				     const Point& pmiddle,
				     const Point& pright) const;

  // --- PROTECTED MEMBER FUNCTIONS ---
protected:
  // -- functions to manage clusters from derived classes --

  // update the cluster of [va,vb], putting vm instead of vb. If
  // reduction=false, the edge [va,vm] is not set reduced.
  void update_cluster(Cluster& c, const Vertex_handle va,
		      const Vertex_handle vb, const Vertex_handle vm,
		      bool reduction = true);

  // get_cluster returns the cluster of [va,vb] in c and return true
  // if it is in a cluster. If erase=true, the cluster is remove from
  // the cluster map
  bool get_cluster(const Vertex_handle va, const Vertex_handle vb,
		   Cluster &c, bool erase = false);


  void add_contrained_edge_to_be_conform(const Vertex_handle& va,
					 const Vertex_handle& vb)
  {
    edges_to_be_conformed.push_back(Constrained_edge(va,vb));
  }
}; // end of Conform_triangulation_2

// --- CONSTRUCTORS ---

template <class Tr>
Conform_triangulation_2<Tr>::
Conform_triangulation_2(const Geom_traits& gt)
  : Tr(gt), is_really_a_constrained_edge(*this),
    edges_to_be_conformed(is_really_a_constrained_edge),
    cluster_map(),
    initialized(false)
{};

// --- ACCESS FUNCTIONS ---

template <class Tr>
int Conform_triangulation_2<Tr>::
number_of_constrained_edges() const
{
  int nedges = 0;
  for(Finite_edges_iterator eit = finite_edges_begin();
      eit != finite_edges_end();
      ++eit)
    if((*eit).first->is_constrained((*eit).second))
      ++nedges;
  return nedges;
}


// --- HELPING FUNCTIONS ---

template <class Tr>
void Conform_triangulation_2<Tr>::
clear() 
{
  cluster_map.clear();
  edges_to_be_conformed.clear();
  Triangulation::clear();
}

// --- CONFORMING FUNCTIONS ---

template <class Tr>
inline
void Conform_triangulation_2<Tr>::
gabriel_conform()
{
  if(!initialized) init(Is_locally_gabriel_conform());
  while( !edges_to_be_conformed.empty() )
    {
      process_one_edge(Is_locally_gabriel_conform());	
    };
}

template <class Tr>
inline
void Conform_triangulation_2<Tr>::
delaunay_conform()
{
  if(!initialized) init(Is_locally_delaunay_conform());
  while( !edges_to_be_conformed.empty() )
    {
      process_one_edge(Is_locally_delaunay_conform());
    };
}

// --- STEP BY STEP FUNCTIONS ---
template <class Tr>
template <class Is_locally_conform>
inline
void Conform_triangulation_2<Tr>::
init(const Is_locally_conform& is_loc_conf)
{
  cluster_map.clear();
  edges_to_be_conformed.clear();
  create_clusters();
  fill_edge_queue(is_loc_conf);
  initialized = true;
}

template <class Tr>
template <class Is_locally_conform>
inline
bool Conform_triangulation_2<Tr>::
refine_step(const Is_locally_conform& is_loc_conf)
{
  if( !edges_to_be_conformed.empty() )
    process_one_edge(is_loc_conf);
  else
    return false;
  return true;
}



// --- PRIVATE MEMBER FUNCTIONS ---

template <class Tr>
void Conform_triangulation_2<Tr>::
create_clusters()
{
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr>
void Conform_triangulation_2<Tr>::
create_clusters_of_vertex(const Vertex_handle v)
{
  Is_edge_constrained test(this, v);
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
      Face_handle f;
      int i;
      is_edge(v, next, f, i); // put in f the face on the right side 
                              // of (v,next)
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
void Conform_triangulation_2<Tr>::
construct_cluster(Vertex_handle v,
		  Constrained_vertex_circulator begin,
		  const Constrained_vertex_circulator& end,
		  Cluster c)
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();

  if(c.vertices.empty())
    {
      c.reduced = false;
      // c.rmin is not initialized because
      // reduced=false!
      c.minimum_squared_length = 
	squared_distance(v->point(), begin->point());
      Constrained_vertex_circulator second(begin);
      ++second;
      c.smallest_angle.first = begin;
      c.smallest_angle.second = second;
    }

  bool all_edges_in_cluster=false; // tell if all incident edges are
  // in the cluster
  if(begin==end)
    all_edges_in_cluster=true;

  const Point& vp = v->point();
  
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
template <class Is_locally_conform>
void Conform_triangulation_2<Tr>::
fill_edge_queue(const Is_locally_conform& is_locally_conform)
{
  for(Finite_edges_iterator ei = finite_edges_begin();
      ei != finite_edges_end();
      ++ei)
    {
      if((*ei).first->is_constrained((*ei).second) && 
	 !is_locally_conform(*this, (*ei).first, (*ei).second) )
	{
	  const Vertex_handle& va = (*ei).first->vertex(cw((*ei).second));
	  const Vertex_handle& vb = (*ei).first->vertex(ccw((*ei).second));
	  edges_to_be_conformed.push_back(std::make_pair(va, vb));
	}
    }
}

//update the encroached segments list
// TODO: perhaps we should remove destroyed edges too
// TODO: rewrite this function one day
template <class Tr>
template <class Is_locally_conform>
void Conform_triangulation_2<Tr>::
update_edges_to_be_conformed(Vertex_handle va,
			     Vertex_handle vb,
			     Vertex_handle vm,
			     const Is_locally_conform& is_locally_conform)
{
  Face_circulator fc = incident_faces(vm), fcbegin(fc);
  if( fc == 0 ) return;

  do {
    for(int i = 0; i<3; i++) {
      if( fc->is_constrained(i) && !is_infinite(fc,cw(i)) &&
	  !is_infinite(fc,ccw(i)) &&
	  !is_locally_conform(*this, fc->handle(), i) )
	{
	  const Vertex_handle& v1 = fc->vertex(ccw(i));
	  const Vertex_handle& v2 = fc->vertex(cw(i));
	  edges_to_be_conformed.push_back(Constrained_edge(v1, v2));
	}
    }
    ++fc;
  } while(fc != fcbegin);

  Face_handle fh;
  int index;
  is_edge(va, vm, fh, index);
  if(!is_locally_conform(*this, fh, index)) {
    edges_to_be_conformed.push_back(Constrained_edge(va, vm));
  }
  is_edge(vb, vm, fh, index);
  if(!is_locally_conform(*this, fh, index)) {
    edges_to_be_conformed.push_back(Constrained_edge(vb, vm));
  }
}

template <class Tr>
template <class Is_locally_conform>
inline
void Conform_triangulation_2<Tr>::
process_one_edge(const Is_locally_conform& is_loc_conf)
{
  Constrained_edge ce = edges_to_be_conformed.front();
  edges_to_be_conformed.pop_front();

  refine_edge(ce.first, ce.second, is_loc_conf);
}

//this function split all the segments that are encroached
template <class Tr>
template <class Is_locally_conform>
void Conform_triangulation_2<Tr>::
refine_edge(Vertex_handle va, Vertex_handle vb,
	    const Is_locally_conform& is_loc_conf)
{
  Face_handle f;
  int i;
  is_edge(va, vb, f, i); // get the edge (f,i)
  CGAL_assertion(f->is_constrained(i));
  
  Cluster c,c2;
  Vertex_handle vm;

  if( get_cluster(va,vb,c,true) )
    if( get_cluster(vb,va,c2,true) )
      { // both ends are clusters
	vm = insert_middle(f,i);
	update_cluster(c,va,vb,vm,false);
	update_cluster(c2,vb,va,vm,false);
      }
    else
      // va only is a cluster
      vm = cut_cluster_edge(va,vb,c);
  else
    if( get_cluster(vb,va,c,true) )
      // vb only is a cluster
      vm = cut_cluster_edge(vb,va,c);
    else
      // no cluster
      vm = insert_middle(f,i);
  update_edges_to_be_conformed(va, vb, vm, is_loc_conf);
};

// template <class Tr>
// inline
// bool Conform_triangulation_2<Tr>::
// is_encroached(const Vertex_handle va, const Vertex_handle vb) const
// {
//   Face_handle fh;
//   int i;
//   is_edge(va,vb,fh,i);

//   const Point& candidat_1 = fh->vertex(i)->point();
//   const Point& candidat_2 = fh->mirror_vertex(i)->point();

//   return ( (/* fh->is_marked() && */
// 	    is_encroached(va, vb, candidat_1) ) ||
// 	   (/* fh->neighbor(i)->is_marked() && */
// 	     is_encroached(va, vb, candidat_2) )
// 	   );
// }

// -> traits?
// TODO, FIXME: not robust!
template <class Tr>
bool Conform_triangulation_2<Tr>::
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

template <class Tr>
typename Conform_triangulation_2<Tr>::Vertex_handle
Conform_triangulation_2<Tr>::
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

  if(c.is_reduced(vb))
    {
      Face_handle fh;
      int i;
      is_edge(va,vb,fh,i);
      vc = insert_middle(fh,i);
    }
  else
    {
      const Point
	& a = va->point(),
	& b = vb->point(),
	& m = midpoint(a, b);



      Vector_2 v = vector(a,m);
      v = scaled_vector(v,CGAL_NTS sqrt(c.minimum_squared_length /
				      squared_distance(a,b)));

      Point i = translate(a,v), i2(i);

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

      vc = virtual_insert_in_the_edge(fh, index, i);
    }
  update_cluster(c, va, vb, vc);
  return vc;
}

template <class Tr>
typename Conform_triangulation_2<Tr>::Vertex_handle
Conform_triangulation_2<Tr>::
insert_middle(Face_handle f, int i)
{
  Construct_midpoint_2
    midpoint = geom_traits().construct_midpoint_2_object();

  const Vertex_handle
    & va = f->vertex(cw(i)),
    & vb = f->vertex(ccw(i));

  const Point& mp = midpoint(va->point(), vb->point());

  Vertex_handle vm = virtual_insert_in_the_edge(f, i, mp);

  return vm;
}

template <class Tr>
inline 
typename Conform_triangulation_2<Tr>::Vertex_handle
Conform_triangulation_2<Tr>::
virtual_insert_in_the_edge(Face_handle fh, int edge_index, const Point& p)
  // insert the point p in the edge (fh, edge_index). It updates seeds 
  // too.
{
  List_of_face_handles zone_of_p;

  // deconstrain the edge before finding the conflicts
  fh->set_constraint(edge_index,false);
  fh->neighbor(edge_index)->set_constraint(fh->mirror_index(edge_index),
					   false);

  get_conflicts_and_boundary(p, 
			     std::back_inserter(zone_of_p), 
			     Emptyset_iterator(), fh);
  
  // reconstrain the edge
  fh->set_constraint(edge_index,true);
  fh->neighbor(edge_index)->set_constraint(fh->mirror_index(edge_index),true);

  Vertex_handle vp = insert(p, Triangulation::EDGE, fh, edge_index);
  // TODO, WARNING: this is not robust!
  // We should deconstrained the constrained edge, insert the two
  // subconstraints and re-constrain them

  return vp;
}

// # used by: is_small_angle, create_clusters_of_vertex
// # compute 4 times the square of the cosine of the angle (ab,ac)
// # WARNING, TODO: this is not exact with doubles and can lead to crashes!!
template <class Tr>
typename Conform_triangulation_2<Tr>::FT
Conform_triangulation_2<Tr>::
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


// --- PROTECTED MEMBER FUNCTIONS ---

template <class Tr>
void Conform_triangulation_2<Tr>::
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

  if(!c.is_reduced())
    {
      typename Cluster::Vertices_map::iterator it = c.vertices.begin();
      while(it!=c.vertices.end() && c.is_reduced(it->first))
	++it; // TODO: use std::find and an object class
      if(it==c.vertices.end())
	c.reduced = true;
    }

  if(c.is_reduced())
    c.rmin = squared_distance(c.smallest_angle.first->point(),
			      c.smallest_angle.second->point())/FT(4);
  cluster_map.insert(std::make_pair(va,c));
}

// # used by refine_face and cut_cluster_edge
template <class Tr>
bool Conform_triangulation_2<Tr>::
get_cluster(Vertex_handle va, Vertex_handle vb, Cluster &c, bool erase)
{
  typedef typename Cluster_map::iterator Iterator;
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


// --- GLOBAL FUNCTIONS ---
template <class Tr>
void
gabriel_conform(Tr& t)
{
  typedef Conform_triangulation_2<Tr> Conform;

  Conform conform;
  conform.swap(t);
  conform.gabriel_conform();
  t.swap(conform);
}

template <class Tr>
void
delaunay_conform(Tr& t)
{
  typedef Conform_triangulation_2<Tr> Conform;

  Conform conform;
  conform.swap(t);
  conform.delaunay_conform();
  t.swap(conform);
}

CGAL_END_NAMESPACE


#endif //CGAL_CONFORM_2_H
