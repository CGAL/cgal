// Copyright (c) 2001-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_H
#define CGAL_MESH_2_H
#include <CGAL/Conforming_Delaunay_triangulation_2.h>
#include <CGAL/Double_map.h>

namespace CGAL {

/**
   Tr is a Delaunay constrained triangulation (with intersections or not)
*/
template <class Tr,
	  class Extras = 
	    Conforming_Delaunay_triangulation_2_default_extras<Tr> >
class Delaunay_mesh_2: public Conforming_Delaunay_triangulation_2<Tr, Extras>
{
public:
  // --- public typedefs ---
  typedef Conforming_Delaunay_triangulation_2<Tr, Extras> Conform;

  typedef Conform Base;
  typedef Delaunay_mesh_2<Tr, Extras> Self;

  /** \name Types inherited from the templated base class */
  //@{
  typedef typename Base::Geom_traits Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef FT      Squared_length;

  typedef typename Base::Vertex_handle        Vertex_handle;
  typedef typename Base::Face_handle          Face_handle;
 
  typedef typename Base::Face_circulator        Face_circulator;
  typedef typename Base::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Base::All_faces_iterator     All_faces_iterator;
  typedef typename Base::Point                  Point;
  //@}

  /** \name Types needed to access member datas
   (Inherited from Conforming_DT.) */
  //@{
  typedef typename Base::Seeds Seeds;
  typedef typename Base::Seeds_const_iterator Seeds_const_iterator;
  typedef Seeds_const_iterator Seeds_iterator;
  //@}

  /** \name Traits types */
  //@{
  typedef typename Geom_traits::Quality Quality;
  //@}

  /** \name CONSTRUCTORS */
  //@{
  
  explicit Delaunay_mesh_2(const Geom_traits& gt = Geom_traits(),
			   const Extras& extras = Extras()):
    Base(gt, extras) {}

  //@}

  /** \name ACCESS FUNCTION */
  //@{
  bool is_bad(const Face_handle fh, Quality& q) const;
  bool is_bad(const Face_handle fh) const ;

  Seeds_const_iterator seeds_begin() const;
  Seeds_const_iterator seeds_end() const;

  //@}

  /** \name HELPING FUNCTION */
  void clear();

  /** \name MARKING FUNCTIONS */
  /** The value type of InputIterator should be Point, and represents
      seeds. Connected components of seeds are marked with the value of 
      "mark". Other components are marked with !mark. The connected
      component of infinite faces is always marked with false.
  */
  template <class InputIterator>
  void set_seeds(InputIterator b, InputIterator e,
		 const bool mark = false,
		 const bool do_it_now = false)
  {
    seeds.clear();
    std::copy(b, e, std::back_inserter(seeds));
    seeds_mark=mark;
    if(do_it_now) mark_facets();
  }

  void clear_seeds()
  {
    seeds.clear();
    seeds_mark = false;
  }

  /** Procedure that forces facets to be marked immediatelly. \em Not
      documented. Call the base classes'one. */
  void mark_facets();

  /** \name MESHING FUNCTIONS */

  // Perform meshing. 
  void refine_mesh();

  /** \name REMESHING FUNCTIONS */

  // Set the geom_traits nut DO NOT recalculate the list of bad faces (must
  // call set_bad_faces of calculate_bad_faces bellow)
  void set_geom_traits(const Geom_traits& gt,
		       bool recalculate_bad_faces = true);

  // Set the geom_traits and add the sequence [begin, end[ to the list
  // of bad faces.
  // Fh_it is a iterator of Face_Handle.
  // Use this overriden function if the list of bad faces can be
  // computed easily without testing all faces.
  template <class Fh_it>
  void set_bad_faces(Fh_it begin, Fh_it end)
  {
    bad_faces.clear();
    for(Fh_it pfit=begin; pfit!=end; ++pfit)
      push_in_bad_faces(*pfit, Quality());
  }

  /** \name STEP BY STEP FUNCTIONS */

  /**
     init(): Initialize the data structures 
     (The call of this function is REQUIRED before any step by step
     operation).
  */
  void init();

  /** Execute one step of the algorithm.
      Needs init() see above */
  bool step_by_step_refine_mesh();

private:
  /** \name PRIVATE TYPES */
  typedef CGAL::Triple<Vertex_handle,
                       Vertex_handle,
                       Vertex_handle> Threevertices; 

  typedef std::list<typename Base::Edge> List_of_edges;
  typedef std::list<Face_handle> List_of_face_handles;

  /** \name traits type */
  typedef typename Geom_traits::Is_bad Is_bad;
  typedef typename Geom_traits::Compute_squared_distance_2
      Compute_squared_distance_2;
//   typedef typename Geom_traits::Construct_circumcenter_2
//       Construct_circumcenter_2;

  /** \name typedefs for private members types */
  typedef CGAL::Double_map<Face_handle, Quality> Bad_faces;

private:
  /** \name PRIVATE MEMBER DATAS */

  // bad_faces: list of bad finite faces
  // warning: some faces could be recycled during insertion in the
  //  triangulation, that's why we need to be able to remoce faces
  //  from the map.
  Bad_faces bad_faces;

  Seeds seeds;
  bool seeds_mark;

private: 
  /** \name PRIVATE MEMBER FUNCTIONS */

  // -- functions that maintain the map of bad faces

  // auxiliary functions called to put a new face in the map, two
  // forms
  void push_in_bad_faces(Face_handle fh, const Quality& q);
 
  // scan all faces and put them if needed in the map
  void fill_facet_map();

  // update the map with faces incident to the vertex v
  void compute_new_bad_faces(Vertex_handle v);

  /** \name inlined functions that compose the refinement process */

  // take one face in the queue and call refine_face
  void process_one_face();

  /** handle one face; call split_face or put in the edges_to_be_conformed the
     list of edges that would be encroached by the circum_center of f
     This function uses Shewchuk's terminator criteria. 
     \todo This function calls get_conflicts_and_boundary and should
     pass the result to split_face. */
  void refine_face(Face_handle f, const Quality& q);

  /** \name functions that really insert points */

  // split the face f by inserting its circum center circum_center
  void split_face(const Face_handle& f, const Point& circum_center);

  // overrideen functions that inserts the point p in the edge
  // (fh,edge_index)
  Vertex_handle virtual_insert_in_the_edge(Face_handle fh,
					   const int edge_index,
					   const Point& p);

  /** \name helping computing functions */ 

  // return the squared length of the triangle corresponding to the
  // face f
  Squared_length shortest_edge_squared_length(Face_handle f);


  /** \name debugging functions */

  bool is_bad_faces_valid();

}; // end of Delaunay_mesh_2

// --- ACCESS FUNCTIONS ---

// ?????????????
// ->traits
// the measure of faces quality
// # We can add here other contraints, such as a bound on the size
template <class Tr, class Extras>
inline
bool Delaunay_mesh_2<Tr, Extras>::
is_bad(const Face_handle f, typename Delaunay_mesh_2<Tr, Extras>::Quality& q) const
{
  const Point
    & a = f->vertex(0)->point(),
    & b = f->vertex(1)->point(),
    & c = f->vertex(2)->point();

  return geom_traits().is_bad_object()(a,b,c,q);
}

template <class Tr, class Extras>
inline
bool Delaunay_mesh_2<Tr, Extras>::
is_bad(const Face_handle f) const
{ 
  Quality q;
  return is_bad(f, q);
}

template <class Tr, class Extras>
inline
typename Delaunay_mesh_2<Tr, Extras>::Seeds_const_iterator
Delaunay_mesh_2<Tr, Extras>::
seeds_begin() const
{
  return seeds.begin();
}

template <class Tr, class Extras>
inline
typename Delaunay_mesh_2<Tr, Extras>::Seeds_const_iterator
Delaunay_mesh_2<Tr, Extras>::
seeds_end() const
{
  return seeds.end();
}



// --- HELPING FUNCTIONS ---

template <class Tr, class Extras>
void Delaunay_mesh_2<Tr, Extras>::
clear() 
{
  bad_faces.clear();
  seeds.clear();
  Base::clear();
}

// --- MARKING FUNCTIONS ---

template <class Tr, class Extras>
void Delaunay_mesh_2<Tr, Extras>::
mark_facets()
{
  Conform::mark_facets(seeds.begin(), seeds.end(), seeds_mark);
}

// --- MESHING FUNCTIONS ---

//the mesh refine function 
template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
refine_mesh()
{
  if(get_initialized() != this->GABRIEL) init();

  while(!Conform::is_conforming_done() || !bad_faces.empty() )
    {
      Conform::make_conforming_Gabriel();
      if ( !bad_faces.empty() )
	process_one_face();
    }
}

// --- REMESHING FUNCTIONS ---

template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
set_geom_traits(const Geom_traits& gt,
		bool recalculate_bad_faces/* = true */)
{
  this->_gt = gt;
  if (recalculate_bad_faces) fill_facet_map();
}

// --- STEP BY STEP FUNCTIONS ---

template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
init()
{
  bad_faces.clear();
  mark_facets(); // facets should be marked before the call to Base::init()

  Base::init_Gabriel();
  // handles clusters and edges

  fill_facet_map();
}

template <class Tr, class Extras>
inline
bool Delaunay_mesh_2<Tr, Extras>::
step_by_step_refine_mesh()
{
  if( !step_by_step_conforming_Gabriel() )
    if ( !bad_faces.empty() )
      process_one_face();
    else
      return false;
  return true;
}

// --- PRIVATE MEMBER FUNCTIONS ---

template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
push_in_bad_faces(Face_handle fh, const Quality& q)
{
  CGAL_assertion(fh->is_marked());
  bad_faces.insert(fh, q);
}

//it is necessarry for process_facet_map
template <class Tr, class Extras>
void Delaunay_mesh_2<Tr, Extras>::
fill_facet_map()
{
  for(Finite_faces_iterator fit = finite_faces_begin();
      fit != finite_faces_end();
      ++fit)
    {
      Quality q;
      if( is_bad(fit, q) && fit->is_marked() )
	push_in_bad_faces(fit, q);
    }
}

template <class Tr, class Extras>
void Delaunay_mesh_2<Tr, Extras>::
compute_new_bad_faces(Vertex_handle v)
{
  Face_circulator fc = v->incident_faces(), fcbegin(fc);
  do {
    Quality q;
    if(!is_infinite(fc))
      if( is_bad(fc, q) && fc->is_marked() )
	push_in_bad_faces(fc, q);
    fc++;
  } while(fc!=fcbegin);
}

template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
process_one_face()
{
  CGAL_assertion(is_bad_faces_valid());

  Face_handle f = bad_faces.front()->second;
  Quality q = bad_faces.front()->first;
  bad_faces.pop_front();
  refine_face(f, q);
}

//split all the bad faces
template <class Tr, class Extras>
void Delaunay_mesh_2<Tr, Extras>::
refine_face(const Face_handle f, const Quality& q)
{
  typename Conform::Is_locally_conforming_Gabriel is_gabriel_conform;

//   Construct_circumcenter_2 circumcenter =
//     geom_traits().construct_circumcenter_2_object();

  const Point& pc = circumcenter(f);

  List_of_edges zone_of_pc_boundary;
  List_of_face_handles zone_of_pc;

  // find conflicts around pc (starting from f as hint)
  get_conflicts_and_boundary(pc, 
			    std::back_inserter(zone_of_pc), 
			    std::back_inserter(zone_of_pc_boundary), 
			    f);
  // For the moment, we don't use the zone_of_pc.
  // It will be used when we will destroyed old bad faces in bad_faces

  bool split_the_face = true;
  bool keep_the_face_bad = false;

  for(typename List_of_edges::iterator it = zone_of_pc_boundary.begin();
      it!= zone_of_pc_boundary.end();
      it++)
    {
      const Face_handle& fh = it->first;
      const int& i = it->second;
      if(fh->is_constrained(i) && !is_gabriel_conform(*this, fh, i, pc))
        {
          const Vertex_handle& va = fh->vertex(cw(i));
          const Vertex_handle& vb = fh->vertex(ccw(i));
          split_the_face = false;
          typename Base::Cluster c,c2;
          bool
            is_cluster_at_va = get_cluster(va,vb,c),
            is_cluster_at_vb = get_cluster(vb,va,c2);
          if( ( is_cluster_at_va &&  is_cluster_at_vb) ||
              (!is_cluster_at_va && !is_cluster_at_vb) )
            {
              // two clusters or no cluster
              add_contrained_edge_to_be_conform(va,vb);
              keep_the_face_bad = true;
            }
          else
            {
              // only one cluster: c or c2
              if(is_cluster_at_vb)
                c = c2;
// What Shewchuk says:
// - If the cluster is not reduced (all segments don't have the same
// length as [va,vb]), then split the edge
// - Else, let rmin be the minimum insertion radius introduced by the
// potential split, let T be the triangle whose circumcenter
// encroaches [va,vb] and let rg be the length of the shortest edge
// of T. If rmin >= rg, then split the edge.

	      if( !c.is_reduced() || 
		  c.rmin >= shortest_edge_squared_length(f) )
		{
		  add_contrained_edge_to_be_conform(va,vb);
		  keep_the_face_bad = true;
		}
	    }
	}
    }; // after here edges encroached by pc are in the list of edges to
       // be conformed.

  if(split_the_face)
    {
      CGAL_assertion(f->is_marked());
      split_face(f, pc);
    }
  else
    if(keep_the_face_bad)
      push_in_bad_faces(f, q);
}

// # used by refine_face
template <class Tr, class Extras>
inline
void Delaunay_mesh_2<Tr, Extras>::
split_face(const Face_handle& f, const Point& circum_center)
{
  CGAL_assertion(f->is_marked());

  List_of_face_handles zone_of_cc;
  List_of_edges zone_of_cc_boundary;

  get_conflicts_and_boundary(circum_center, 
			     std::back_inserter(zone_of_cc),
			     std::back_inserter(zone_of_cc_boundary),
			     f);
  CGAL_assertion(is_bad_faces_valid());
  for(typename List_of_face_handles::iterator fh_it = zone_of_cc.begin();
      fh_it != zone_of_cc.end();
      ++fh_it)
    bad_faces.erase(*fh_it);

  extras().signal_before_inserted_vertex_in_face(static_cast<const Tr&>(*this),
						 f,
						 zone_of_cc_boundary.begin(),
						 zone_of_cc_boundary.end(),
						 zone_of_cc.begin(),
						 zone_of_cc.end(),
						 circum_center);

  // insert the point in the triangulation with star_hole
  Vertex_handle v = star_hole(circum_center,
			      zone_of_cc_boundary.begin(),
			      zone_of_cc_boundary.end(),
			      zone_of_cc.begin(),
			      zone_of_cc.end());

  extras().signal_after_inserted_vertex_in_face(static_cast<const Tr&>(*this),
					      v);

  Face_circulator fc = incident_faces(v), fcbegin(fc);
  do {
    fc->set_marked(true);
  } while (++fc != fcbegin);

  compute_new_bad_faces(v);
}

template <class Tr, class Extras>
inline 
typename Delaunay_mesh_2<Tr, Extras>::Vertex_handle
Delaunay_mesh_2<Tr, Extras>::
virtual_insert_in_the_edge(Face_handle fh, int edge_index, const Point& p)
  // insert the point p in the edge (fh, edge_index). It updates seeds 
  // too.
{
  const Vertex_handle& va = fh->vertex( cw(edge_index));
  const Vertex_handle& vb = fh->vertex(ccw(edge_index));

  bool 
    mark_at_right = fh->is_marked(),
    mark_at_left = fh->neighbor(edge_index)->is_marked();

  List_of_face_handles zone_of_p;

  // deconstrain the edge
  remove_constrained_edge(fh, edge_index);

  get_conflicts_and_boundary(p, 
			     std::back_inserter(zone_of_p), 
			     Emptyset_iterator(), fh);

  for(typename List_of_face_handles::iterator fh_it = zone_of_p.begin();
      fh_it != zone_of_p.end();
      ++fh_it)
    bad_faces.erase(*fh_it);

  extras().signal_before_inserted_vertex_in_edge(static_cast<const Tr&>(*this),
					       fh, edge_index, p);

  Vertex_handle vp = insert(p, fh);

  // re-insert the two constrained edges
  insert_constraint(va, vp);
  insert_constraint(vp, vb);

  extras().signal_after_inserted_vertex_in_edge(static_cast<const Tr&>(*this),
						vp);

  // now, let's update 'in-domain' markers
  int dummy;
  // if we put edge_index instead of dummy, Intel C++ does not find
  // a matching function for is_edge
  is_edge(va, vp, fh, dummy); 
  // set fh to the face at the right of [va,vp]

  Face_circulator fc = incident_faces(vp, fh), fcbegin(fc);
  // circulators are counter-clockwise, so we start at the right of
  // [va,vp]
  do {
    if( !is_infinite(fc) )
      fc->set_marked(mark_at_right);
    ++fc;
  } while ( fc->vertex(ccw(fc->index(vp))) != vb );
  // we are now at the left of [va,vb]
  do {
    if( !is_infinite(fc) )
      fc->set_marked(mark_at_left);
    ++fc;
  } while ( fc != fcbegin );

  // then let's update bad faces
  compute_new_bad_faces(vp);

  return vp;
}


template <class Tr, class Extras>
bool Delaunay_mesh_2<Tr, Extras>::
is_bad_faces_valid()
{
  typedef std::list<std::pair<Quality, Face_handle> > Bad_faces_list;
  
  bool result = true;

  Bad_faces_list bad_faces_list;

  while(!bad_faces.empty())
    {
      Quality d = bad_faces.front()->first;
      Face_handle fh = bad_faces.front()->second;
      bad_faces.pop_front();
      
      bad_faces_list.push_back(std::make_pair(d, fh));

      const Vertex_handle
	& va = fh->vertex(0),
	& vb = fh->vertex(1),
	& vc = fh->vertex(2);

      Face_handle fh2;
      Quality q;
      bool face = ( is_face(va, vb, vc, fh2) && fh == fh2 );
      bool marked = fh->is_marked();
      bool bad = is_bad(fh, q);
      if( ! ( face && marked && bad ) )
	{
	  result = false;
	  std::cerr << "Invalid bad face: (" << va->point() << ", "
		    << vb->point() << ", " << vc->point() << ")" << std::endl;
	  if( ! face )
	    std::cerr << "(not a face)" << std::endl;
	  if( ! marked )
	    std::cerr << "(not marked)" << std::endl;
	  if( ! bad )
	    std::cerr << "(not bad, quality=" << q << ")" << std::endl;
	}
    }

  for(typename Bad_faces_list::iterator it = bad_faces_list.begin();
      it != bad_faces_list.end();
      ++it)
    bad_faces.insert(it->second, it->first);
  
  return result;
}

// ->traits?
//the shortest edge that are in a triangle
// # used by: refine_face, squared_minimum_sine
template <class Tr, class Extras>
typename Delaunay_mesh_2<Tr, Extras>::FT
Delaunay_mesh_2<Tr, Extras>::
shortest_edge_squared_length(Face_handle f)
{
  Compute_squared_distance_2 squared_distance = 
    geom_traits().compute_squared_distance_2_object();
  const Point 
    & pa = (f->vertex(0))->point(),
    & pb = (f->vertex(1))->point(),
    & pc = (f->vertex(2))->point();
  FT a, b, c;
  a = squared_distance(pb, pc);
  b = squared_distance(pc, pa);
  c = squared_distance(pa, pb);
  return (min(a, min(b, c)));
}

// --- GLOBAL FUNCTIONS ---

// this a workaround, used just below, to fix a bug in Sun CC.
template <class Tr>
struct Refine_mesh_2_default_argument_helper : public Tr::Geom_traits {};

template <class Tr>
void
refine_Delaunay_mesh_2(Tr& t,
		       const typename Tr::Geom_traits& gt
		       = Refine_mesh_2_default_argument_helper<Tr>() )
{
  typedef Delaunay_mesh_2<Tr,
    Conforming_Delaunay_triangulation_2_default_extras<Tr> > Mesh;

  Mesh mesh;
  mesh.swap(t);
  mesh.set_geom_traits(gt);
  mesh.refine_mesh();
  t.swap(mesh);
}

template <class Tr, typename InputIterator>
void
refine_Delaunay_mesh_2(Tr& t,
		       InputIterator b, InputIterator e,
		       bool mark = false,
		       const typename Tr::Geom_traits& gt
		       = Refine_mesh_2_default_argument_helper<Tr>())
{
  typedef Delaunay_mesh_2<Tr,
    Conforming_Delaunay_triangulation_2_default_extras<Tr> > Mesh;

  Mesh mesh;
  mesh.swap(t);
  mesh.set_geom_traits(gt);
  mesh.set_seeds(b, e, mark);
  mesh.refine_mesh();
  t.swap(mesh);
}

}

#endif // CGAL_MESH_2_H
