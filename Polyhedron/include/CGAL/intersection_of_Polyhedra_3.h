// Copyright (c) 2011 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERSECTION_OF_POLYHEDRA_3_H
#define CGAL_INTERSECTION_OF_POLYHEDRA_3_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>

#include <boost/next_prior.hpp>
#include <set>
#include <vector>
#include <list>
#include <algorithm>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>

#include <CGAL/Modifier_base.h>
#include <CGAL/internal/corefinement/Polyhedron_constness_types.h>
#include <CGAL/internal/corefinement/intersection_triangle_segment_3.h>
#include <CGAL/internal/corefinement/intersection_coplanar_triangles_3.h>

#include <boost/utility.hpp>
#include <boost/type_traits/is_floating_point.hpp>

#ifdef CGAL_COREFINEMENT_DEBUG
#warning look at CGAL/Mesh_3/Robust_intersection_traits.h and the statically filtered decision tree
#endif

namespace CGAL{

//This functor computes the pairwise intersection of polyhedral surfaces.
//Intersection are given as a set of polylines
//The algorithm works as follow:
//From each polyhedral surface we can get it as a set of segments or as a set of triangles.
//We first use Box_intersection_d to filter intersection between all polyhedral
//surface segments and polyhedral triangles.
//From this filtered set, for each pair (segment,triangle), we look at the
//intersection type. If not empty, we can have three different cases
//  1)the segment intersect the interior of the triangle:
//        We compute the intersection point and for each triangle incident
//        to the segment, we write the fact that the point belong to the intersection
//        of these two triangles.
//  2)the segment intersect the triangle on an edge
//        We do the same thing as described above but
//        for all triangle incident to the edge intersected
//  3)the segment intersect the triangle at a vertex
//        for each edge incident to the vertex, we do
//        the same operations as in 2)
//
//In case the segment intersect the triangle at one of the segment endpoint,  
//we repeat the same procedure for each segment incident to this
//endpoint.
//
//Note that given a pair (segment,triangle)=(S,T), if S belongs 
//to the plane of T, we have nothing to do in the following cases:
//  -- no triangle T' contains S such that T and T' are coplanar
//  -- at least one triangle contains S
// Indeed, the intersection points of S and T will be found using segments
// of T or segments adjacent to S.


namespace internal_IOP {
  //an enum do decide which kind of intersection points are needed
  struct No_predicates_on_constructions{};
  struct Predicates_on_constructions{};
} // namespace internal_IOP

template<class Polyhedron>
struct Empty_node_visitor{
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_handle;
  void new_node_added(int,internal_IOP::Intersection_type,Halfedge_handle,Halfedge_handle,bool,bool){}
  template<class Iterator>
  void annotate_graph(Iterator,Iterator){}
  void update_terminal_nodes(std::vector<bool>&){}
  void set_number_of_intersection_points_from_coplanar_facets(int){};
  void add_filtered_intersection(Halfedge_handle,Halfedge_handle,const Polyhedron&,const Polyhedron&){}
  void start_new_polyline(int,int){}
  void add_node_to_polyline(int){}
  void new_input_polyhedron(const Polyhedron&){}
  template<class T>
  void finalize(T&){}
  typedef internal_IOP::No_predicates_on_constructions Node_storage_type;
  typedef Tag_true Is_polyhedron_const;
  static const bool do_need_vertex_graph = false;
};

namespace internal_IOP{
  
template <class Polyhedron,class Is_const>
struct Compare_handles{
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle   Halfedge_handle;
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge          Halfedge;
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Vertex            Vertex;
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Facet             Facet;
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Facet_handle      Facet_handle;
  typedef std::pair<const Vertex*,const Vertex*>                            Vertex_handle_pair;
  
  static inline Vertex_handle_pair 
  make_sorted_pair_of_vertices(Halfedge_handle h) {
    const Vertex* v1=&(* h->vertex() );
    const Vertex* v2=&(* h->opposite()->vertex() );
    if ( v1 < v2 )
      return Vertex_handle_pair(v1,v2);
    return Vertex_handle_pair(v2,v1);
  }
  
  bool operator()(Halfedge_handle h1,Halfedge_handle h2) const {
    Vertex_handle_pair p1=make_sorted_pair_of_vertices(h1);
    Vertex_handle_pair p2=make_sorted_pair_of_vertices(h2);
    return  p1 < p2; 
  }
  
  bool operator()(Facet_handle f1,Facet_handle f2) const {
    return &(*f1) < &(*f2);
  }
  
  bool operator()(const std::pair<Halfedge_handle,Polyhedron*>& p1, const std::pair<Halfedge_handle,Polyhedron*>& p2) const{
    Halfedge* h1= (std::min) ( &(*(p1.first)), &(*(p1.first->opposite())) );
    Halfedge* h2= (std::min) ( &(*(p2.first)), &(*(p2.first->opposite())) );
    return h1<h2;
  }
};

//could not put it in the above struct (gcc complains about an ambiguous call)
template <class Polyhedron,class Is_const>
struct Compare_handle_pairs{
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Facet             Facet;
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Facet_handle      Facet_handle;  
  typedef std::pair<Facet_handle,Facet_handle>                              Facet_pair;
  typedef std::pair<Facet_pair,int>                                         Facet_pair_and_int;
  
  bool operator()(const Facet_pair& p1, const Facet_pair& p2) const{
    Facet* f1=&(*p1.first);
    Facet* f2=&(*p2.first);
    if (f1==f2){
      f1=&(*p1.second);
      f2=&(*p2.second);
    }
    return f1<f2;
  }

  bool operator()(const Facet_pair_and_int& p1, const Facet_pair_and_int& p2) const{
    Facet* f1=&(*p1.first.first);
    Facet* f2=&(*p2.first.first);
    if (f1==f2){
      f1=&(*p1.first.second);
      f2=&(*p2.first.second);
    }
    if (f1<f2) return true;
    if (f1>f2) return false;
    return p1.second<p2.second;
  }
};

template<class Polyhedron,class Node_vector,class Is_const>
struct Order_along_a_halfedge{
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle Halfedge_handle;
  const Node_vector& nodes;
  Halfedge_handle hedge;
  
  Order_along_a_halfedge(Halfedge_handle hedge_,const Node_vector& nodes_):nodes(nodes_),hedge(hedge_){}
  bool operator()(int i,int j) const {
    //returns true, iff q lies strictly between p and r.
    try{
      typename Node_vector::Protector p;
      return CGAL::collinear_are_strictly_ordered_along_line(nodes.to_interval(hedge->vertex()->point()),
                                                             nodes.interval_node(j),
                                                             nodes.interval_node(i));
    }
    catch(CGAL::Uncertain_conversion_exception&){
      return CGAL::collinear_are_strictly_ordered_along_line(nodes.to_exact(hedge->vertex()->point()),
                                                             nodes.exact_node(j),
                                                             nodes.exact_node(i));      
    }
  }
};


template <class HDS>
class Split_halfedge_at_point : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Vertex_handle   Vertex_handle;
  typedef typename HDS::Vertex          Vertex;
  Halfedge_handle hedge;
  Vertex          vertex;
  
  typename HDS::Halfedge::Base*
  unlock_halfedge(Halfedge_handle h){
    return static_cast<typename HDS::Halfedge::Base*>(&(*h));
  }
  
public:
  
  template <class Point_3>
  Split_halfedge_at_point( Halfedge_handle h,const Point_3& point):hedge(h),vertex(point){}

  void operator()( HDS& hds) {
    
    Vertex_handle v=hds.vertices_push_back(vertex);
    Halfedge_handle opposite=hedge->opposite();
    
    Halfedge_handle new_hedge=hds.edges_push_back(*hedge);
    Halfedge_handle new_opposite=new_hedge->opposite();
    
    //update next relations
    unlock_halfedge(new_hedge)->set_next(hedge);
    unlock_halfedge(new_hedge->prev())->set_next(new_hedge);
    unlock_halfedge(hedge)->set_prev(new_hedge);

    unlock_halfedge(opposite)->set_next(new_opposite);
    unlock_halfedge(new_opposite)->set_prev(opposite);    
    unlock_halfedge(new_opposite->next())->set_prev(new_opposite);

    unlock_halfedge(opposite)->set_vertex(v);
    unlock_halfedge(new_hedge)->set_vertex(v);

    v->set_halfedge(new_hedge);
    new_opposite->vertex()->set_halfedge(new_opposite);
  }
};



} //namespace internal_IOP



//WARNING THIS IS DONE ONLY FOR POLYHEDRON
template<class Polyhedron,class Halfedge_predicate,
         class Set_vertex_corner, class Kernel=typename Polyhedron::Traits::Kernel>
class Node_visitor_for_polyline_split{
//typedefs  
  typedef typename Polyhedron::Halfedge_handle                         Halfedge_handle;
  typedef typename Polyhedron::Halfedge                                Halfedge;
  typedef typename Polyhedron::Vertex_handle                           Vertex_handle;
  //Info stores information about a particular intersection on an
  //edge or a vertex of a polyhedron. The two first elements in 
  //the template describe the intersected simplex of the considered 
  //polyhedron; the two last elements describe the element of the 
  //second polyhedron (can be either a vertex, an edge of a facet)
  //involved in the intersection
  typedef CGAL::cpp0x::tuple<internal_IOP::Intersection_type,
                             Halfedge_handle,
                             internal_IOP::Intersection_type,
                             Halfedge_handle>                          Info;
  typedef std::map<Halfedge*,Polyhedron*>                              Hedge_to_polyhedron_map;
  typedef std::vector<Info>                                            Infos;
  typedef std::map<int,Infos >                                         Node_to_infos_map;
//data members
  Node_to_infos_map   node_infos;
  Hedge_to_polyhedron_map hedge_to_polyhedron;
  Halfedge_predicate is_on_polyline;
  Set_vertex_corner set_as_corner;
//functions  
  void handle_principal_edge(int node_id,
                             internal_IOP::Intersection_type type,
                             Halfedge_handle principal_edge,
                             Halfedge_handle additional_edge,
                             bool is_vertex_coplanar,
                             bool is_vertex_opposite_coplanar)
  {
    bool coplanar_v=false;
    if (is_vertex_coplanar) coplanar_v=true;
    else if (is_vertex_opposite_coplanar){
      principal_edge=principal_edge->opposite();
      coplanar_v=true;
    }
    
    if (coplanar_v)
      handle_on_vertex(node_id,principal_edge,type,additional_edge);
    else{
      if ( is_on_polyline(principal_edge) ){
        typename Node_to_infos_map::iterator it_res=
          node_infos.insert(std::make_pair(node_id,Infos())).first; 
          it_res->second.push_back( Info(internal_IOP::EDGE,principal_edge,type,additional_edge) );
      }
    }
  }
  
  void handle_on_vertex(int node_id,Halfedge_handle edge,internal_IOP::Intersection_type type,Halfedge_handle additional_edge){
    Halfedge_handle current=edge;
    do{
      if (is_on_polyline(current)){
        typename Node_to_infos_map::iterator it_res=
          node_infos.insert(std::make_pair(node_id,Infos())).first; 
          it_res->second.push_back( Info(internal_IOP::VERTEX,current,type,additional_edge) );        
        break;
      }
      current=current->next()->opposite();
    }
    while(current!=edge);    
  }
  
  Halfedge* make_unique_key(Halfedge_handle h){
    if (&(*h) < &(*h->opposite()))
      return &(*h);
    else
      return &(*h->opposite());
  }
  
  //   new_hedge    hedge
  //  ----------->   ----------->
  //               v
  //  <-----------   <-----------
  //   new_opposite     opposite 
  //  
  void split_edge_and_retriangulate(Halfedge_handle hedge,const typename Kernel::Point_3& point,Polyhedron& P){
    internal_IOP::Split_halfedge_at_point<typename Polyhedron::HalfedgeDS> delegated(hedge,point);
    P.delegate( delegated );
    CGAL_assertion(P.is_valid());
    //triangulate the two adjacent facets
    if (!hedge->is_border())
      P.split_facet(hedge->prev(),hedge->next());
    if (!hedge->opposite()->is_border())
      P.split_facet(hedge->opposite(),hedge->opposite()->next()->next());
    CGAL_assertion(P.is_valid());
  }

  //sort node ids so that we can split the hedge
  //consecutively
  template <class Node_vector>
  void sort_vertices_along_hedge(std::vector<int>& node_ids,Halfedge_handle hedge,const Node_vector& nodes)
  {
    std::sort(node_ids.begin(),
              node_ids.end(),
              internal_IOP::Order_along_a_halfedge<Polyhedron,Node_vector,Is_polyhedron_const>(hedge,nodes)
    );
  }
  
public:
  static const bool do_need_vertex_graph = false;  
  typedef internal_IOP::Predicates_on_constructions  Node_storage_type;  
  typedef Tag_false Is_polyhedron_const;

  Node_visitor_for_polyline_split(){}
  Node_visitor_for_polyline_split(const Halfedge_predicate& getting,
                                  const Set_vertex_corner& setting)
    :is_on_polyline(getting),set_as_corner(setting){}

  void new_node_added(int node_id,
                      internal_IOP::Intersection_type type,
                      Halfedge_handle principal_edge,
                      Halfedge_handle additional_edge,
                      bool is_vertex_coplanar,
                      bool is_vertex_opposite_coplanar)
  {
    switch(type)
    {
      case internal_IOP::FACET: //Facet intersected by an edge
        handle_principal_edge(node_id,type,principal_edge,additional_edge,is_vertex_coplanar,is_vertex_opposite_coplanar);
        break;
      case internal_IOP::EDGE: //Edge intersected by an edge
        handle_principal_edge(node_id,type,principal_edge,additional_edge,is_vertex_coplanar,is_vertex_opposite_coplanar);
        if ( is_on_polyline(additional_edge) ){
          typename Node_to_infos_map::iterator it_res=
            node_infos.insert(std::make_pair(node_id,Infos())).first; 
            it_res->second.push_back( Info(type,additional_edge,
                                           ( is_vertex_coplanar||is_vertex_opposite_coplanar ) ? internal_IOP::VERTEX:internal_IOP::EDGE,
                                           is_vertex_opposite_coplanar?principal_edge->opposite():principal_edge) );
        }
        break;
      case internal_IOP::VERTEX://Vertex intersected by an edge
        handle_principal_edge(node_id,type,principal_edge,additional_edge,is_vertex_coplanar,is_vertex_opposite_coplanar);
        handle_on_vertex( node_id,additional_edge,
                          ( is_vertex_coplanar||is_vertex_opposite_coplanar ) ? internal_IOP::VERTEX:internal_IOP::EDGE,
                          is_vertex_opposite_coplanar?principal_edge->opposite():principal_edge);
        break;
      default:
        break;
    }
  }
  
  template<class Iterator>
  void annotate_graph(Iterator begin,Iterator end){
    for(Iterator it=begin;it!=end;++it){
      typename Node_to_infos_map::iterator it_res=node_infos.find(it->first);
      if (it_res!=node_infos.end())
        it->second.make_terminal();
    }
  }
  
  void update_terminal_nodes(std::vector<bool>& terminal_bools){
    for (typename Node_to_infos_map::iterator it=node_infos.begin();it!=node_infos.end();++it){
      terminal_bools[it->first]=true;
    }
  }
  
  void new_input_polyhedron(const Polyhedron&){}
  void start_new_polyline(int,int){}
  void add_node_to_polyline(int){}
  void set_number_of_intersection_points_from_coplanar_facets(int){}
  
  void add_filtered_intersection(Halfedge_handle eh,Halfedge_handle fh,Polyhedron& Pe,Polyhedron& Pf){
    hedge_to_polyhedron.insert(std::make_pair(make_unique_key(eh),&Pe));
    hedge_to_polyhedron.insert(std::make_pair(make_unique_key(fh),&Pf));
    hedge_to_polyhedron.insert(std::make_pair(make_unique_key(fh->next()),&Pf));
    hedge_to_polyhedron.insert(std::make_pair(make_unique_key(fh->next()->next()),&Pf));
  }
    
  //split_halfedges
  template <class Node_vector>
  void finalize(const Node_vector& nodes){
    typedef std::map<std::pair<Halfedge_handle,Polyhedron*>,
                     std::vector<int>,internal_IOP::Compare_handles<Polyhedron,Is_polyhedron_const> >  Halfedges_to_split;
    
    Halfedges_to_split halfedges_to_split;
    
    for (typename Node_to_infos_map::iterator it=node_infos.begin();it!=node_infos.end();++it){
      int node_id=it->first;
      const Infos& infos=it->second;
      std::map<Polyhedron*,std::vector<Halfedge_handle> > hedges_to_split;
      
      //collect information about halfedge to split
      typename Infos::const_iterator it_info=infos.begin();
      for (;it_info!=infos.end();++it_info)
      {
        typename Hedge_to_polyhedron_map::iterator  it_poly=
          hedge_to_polyhedron.find(make_unique_key(CGAL::cpp0x::get<1>(*it_info)));
        CGAL_assertion(it_poly!=hedge_to_polyhedron.end());
        //associate information to an intersection point:
        //we give which simplex of the other polyhedron intersect the simplex considered
        set_as_corner.add_info_to_node(node_id,it_poly->second,*it_info);
        switch(CGAL::cpp0x::get<0>(*it_info))
        {
          case internal_IOP::EDGE:
          {
            halfedges_to_split.insert(
              std::make_pair( std::make_pair(CGAL::cpp0x::get<1>(*it_info),&(*(it_poly->second))),std::vector<int>() )
            ).first->second.push_back(node_id);
          break;
          }
          case internal_IOP::VERTEX:
            set_as_corner(CGAL::cpp0x::get<1>(*it_info)->vertex(),node_id,it_poly->second);
          break;
          default:
            CGAL_assertion(false);
            //should never be here
        }
      }
    }
    
    
    //do the split
    for(typename Halfedges_to_split::iterator it=halfedges_to_split.begin();it!=halfedges_to_split.end();++it){
      Halfedge_handle hedge=it->first.first;
      Polyhedron* P=it->first.second;
      std::vector<int>& node_ids=it->second;
      
      sort_vertices_along_hedge(node_ids,hedge,nodes);
      for (std::vector<int>::iterator it_id=node_ids.begin();it_id!=node_ids.end();++it_id){
        split_edge_and_retriangulate(hedge,nodes[*it_id],*P);
        set_as_corner(hedge->opposite()->vertex(),*it_id,(Polyhedron*)(0));
      }
    }
  }
};


namespace internal_IOP{
  
  template <class Polyhedron,class In_kernel,class Exact_kernel>
  typename Exact_kernel::Point_3
  compute_triangle_segment_intersection_point(
    typename Polyhedron::Vertex_const_handle vh1,typename Polyhedron::Vertex_const_handle vh2,
    typename Polyhedron::Vertex_const_handle vf1,typename Polyhedron::Vertex_const_handle vf2,typename Polyhedron::Vertex_const_handle vf3,
    const Exact_kernel& ek)       
  {
    CGAL::Cartesian_converter<In_kernel,Exact_kernel> to_exact;
    typename Exact_kernel::Triangle_3 t(to_exact( vf1->point() ),
                                        to_exact( vf2->point() ),
                                        to_exact( vf3->point() )
    );
    
    typename Exact_kernel::Segment_3 s (to_exact( vh1->point() ),
                                        to_exact( vh2->point() )
    );
    
    typename Exact_kernel::Intersect_3 exact_intersect=ek.intersect_3_object();
    CGAL::Object inter=exact_intersect(t,s);
    CGAL_assertion(CGAL::do_intersect(t,s));
    const typename Exact_kernel::Point_3* e_pt=CGAL::object_cast<typename Exact_kernel::Point_3>(&inter);
    CGAL_assertion(e_pt!=NULL);
    return *e_pt;
  }
  
  template <class Polyhedron,class In_kernel,class Exact_kernel>
  typename Exact_kernel::Point_3
  compute_triangle_segment_intersection_point(
    typename Polyhedron::Halfedge_const_handle edge,
    typename Polyhedron::Facet_const_handle facet,
    const Exact_kernel& ek) 
  {
    return compute_triangle_segment_intersection_point<Polyhedron,In_kernel,Exact_kernel>(
            edge->vertex(),edge->opposite()->vertex(),
            facet->halfedge()->vertex(),facet->halfedge()->next()->vertex(),facet->halfedge()->opposite()->vertex(),
            ek);
            
  }    
    
  
  //A class containing a vector of the intersection points.
  //The third template parameter indicates whether an
  //exact representation is required
  template <class Polyhedron,class Kernel,class Node_storage,bool Has_exact_constructions=!boost::is_floating_point<typename Kernel::FT>::value>
  class Triangle_segment_intersection_point;
  
  
  //Store only the double version of the intersection points.
  template <class Polyhedron,class Kernel>
  class Triangle_segment_intersection_point<Polyhedron,Kernel,No_predicates_on_constructions,false>
  {
  //typedefs
    typedef std::vector <typename Kernel::Point_3>             Node_vector;
    typedef typename Polyhedron::Halfedge_const_handle         Halfedge_handle;
    typedef typename Polyhedron::Facet_const_handle            Facet_handle;
    typedef CGAL::Exact_predicates_exact_constructions_kernel  Exact_kernel;
    typedef CGAL::Cartesian_converter<Exact_kernel,Kernel>     Exact_to_double;    
  //members
    Node_vector nodes;
    Exact_kernel ek;
    Exact_to_double exact_to_double;
  public:
    typedef CGAL::Interval_nt<true>::Protector                 Protector;
  
    const typename Kernel::Point_3&
    operator[](int i) const {
      return nodes[i];
    }
    
    const typename Kernel::Point_3& exact_node(int i) const {return nodes[i];}
    const typename Kernel::Point_3& interval_node(int i) const {return nodes[i];}
    const typename Kernel::Point_3& to_exact(const typename Kernel::Point_3& p) const {return p;}
    const typename Kernel::Point_3& to_interval(const typename Kernel::Point_3& p) const {return p;}
    
    size_t size() const {return nodes.size();}
    
    //add a new node in the final graph.
    //it is the intersection of the triangle with the segment
    void add_new_node(Halfedge_handle edge,Facet_handle facet)
    {
      nodes.push_back (  exact_to_double(
        compute_triangle_segment_intersection_point<Polyhedron,Kernel>(edge,facet,ek)
      ));
    }
    
    void add_new_node(const typename Exact_kernel::Point_3& p)
    {
      nodes.push_back(  exact_to_double(p) );
    }    

    void add_new_node(const typename Kernel::Point_3& p)
    {
      nodes.push_back(p);
    }    
  };
  
  
  //second specializations: store an exact copy of the points so that we can answer exactly predicates
  //FYI, it used to have two specializations (one in the case the polyhedron
  //can be edited and on if it cannot) building exact representation on demand.
  //In the former case, we were using facet and halfedge while in the latter
  //triple of vertex_handle and pair of vertex_handle
  template <class Polyhedron,class Kernel>
  class Triangle_segment_intersection_point<Polyhedron,Kernel,Predicates_on_constructions,false>
  {
  //typedefs
  public: 
    typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >  Ikernel;
    typedef CGAL::Exact_predicates_exact_constructions_kernel  Exact_kernel;
  private:
    typedef CGAL::Cartesian_converter<Ikernel,Kernel>          Interval_to_double;
    typedef CGAL::Cartesian_converter<Kernel,Ikernel>          Double_to_interval;
    typedef CGAL::Cartesian_converter<Exact_kernel,Ikernel>    Exact_to_interval;
    typedef CGAL::Cartesian_converter<Kernel,Exact_kernel>     Double_to_exact;
  
    typedef typename Polyhedron::Vertex_const_handle           Vertex_handle;
    typedef typename Polyhedron::Halfedge_const_handle         Halfedge_handle;
    typedef typename Polyhedron::Facet_const_handle            Facet_handle; 
    
    typedef std::vector <Ikernel::Point_3>      Interval_nodes;
    typedef std::vector <Exact_kernel::Point_3> Exact_nodes;
 
        
  //members
    Interval_nodes  inodes;
    Exact_nodes enodes;
  
    Interval_to_double  interval_to_double;
    Exact_to_interval   exact_to_interval;
    Double_to_interval  double_to_interval;
    Double_to_exact double_to_exact;
    Exact_kernel        ek;
    
  public:
    typedef CGAL::Interval_nt<false>::Protector                 Protector;  
  
    typename Kernel::Point_3
    operator[](int i) const {
      return interval_to_double(inodes[i]);
    }
    
    const typename Ikernel::Point_3&
    interval_node(int i) const {
      return inodes[i];
    }
    
    typename Ikernel::Point_3
    to_interval(const typename Kernel::Point_3& p) const {
      return double_to_interval(p);
    }
    
    const Exact_kernel::Point_3
    exact_node(int i) const {
      return enodes[i];
    }

    typename Exact_kernel::Point_3
    to_exact(const typename Kernel::Point_3& p) const {
      return double_to_exact(p);
    }    
    
    
    size_t size() const {return enodes.size();}

    void add_new_node(Halfedge_handle edge,Facet_handle facet)
    {
      enodes.push_back(compute_triangle_segment_intersection_point<Polyhedron,Kernel>(edge,facet,ek) );
      inodes.push_back( exact_to_interval(enodes.back()) );
    }

    void add_new_node(const Exact_kernel::Point_3& p){
      enodes.push_back(p);
      inodes.push_back( exact_to_interval(p) );
    }

    //the point is an input
    void add_new_node(const typename Kernel::Point_3& p){
      enodes.push_back(to_exact(p));
      inodes.push_back( double_to_interval(p) );
    }
  };  
  
  //Third specialization: The kernel already has exact constructions.
  template <class Polyhedron,class Kernel,class Node_storage>
  class Triangle_segment_intersection_point<Polyhedron,Kernel,Node_storage,true>
  {
  //typedefs
    typedef std::vector <typename Kernel::Point_3>             Node_vector;
    typedef typename Polyhedron::Halfedge_const_handle         Halfedge_handle;
    typedef typename Polyhedron::Facet_const_handle            Facet_handle;
  //members
    Node_vector nodes;
    Kernel k;
  public:
    typedef Kernel Ikernel;
    typedef Kernel Exact_kernel;
    typedef void* Protector;
    const typename Kernel::Point_3&
    operator[](int i) const {
      return nodes[i];
    }
   
    size_t size() const {return nodes.size();}
    const typename Kernel::Point_3& exact_node(int i) const {return nodes[i];}
    const typename Kernel::Point_3& interval_node(int i) const {return nodes[i];}
    
    //add a new node in the final graph.
    //it is the intersection of the triangle with the segment
    void add_new_node(Halfedge_handle edge,Facet_handle facet)
    {
      nodes.push_back (  
        compute_triangle_segment_intersection_point<Polyhedron,Kernel>(edge,facet,k)
      );
    }

    void add_new_node(const typename Kernel::Point_3& p)
    {
      nodes.push_back(p);
    }
    
    const typename Kernel::Point_3& to_interval(const typename Kernel::Point_3& p) const { return p; }
    const typename Kernel::Point_3& to_exact(const typename Kernel::Point_3& p) const { return p; }

  };  
  
}

//TODO an important requirement is that the Polyhedron should be based on a list-based
//HDS. We use a lot of maps that use the address of Facet,Halfedge and a reallocation would
//be dramatic.

template< class Polyhedron,
          class Kernel=typename Polyhedron::Traits::Kernel,
          class Node_visitor=Empty_node_visitor<Polyhedron>,
          class Node_storage_type=typename Node_visitor::Node_storage_type,
          class Use_const_polyhedron=typename Node_visitor::Is_polyhedron_const
         >
class Intersection_of_Polyhedra_3{

//typedefs  
  typedef typename Kernel::Triangle_3                        Triangle;
  typedef typename Kernel::Segment_3                         Segment;
  typedef internal_IOP::
    Polyhedron_types<Polyhedron,Use_const_polyhedron>        Polyhedron_types;
  
  typedef typename Polyhedron_types::Polyhedron_ref          Polyhedron_ref;
  typedef typename Polyhedron_types::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron_types::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron_types::Facet_iterator          Facet_iterator;
  typedef typename Polyhedron_types::Facet_handle            Facet_handle;
  typedef typename Polyhedron_types::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron_types::Vertex                  Vertex;
  typedef typename Polyhedron_types::Facet                   Facet;
  typedef CGAL::Box_intersection_d::Box_with_handle_d<
            double, 3, Halfedge_handle>                      Box;
  typedef std::pair<Facet_handle,Facet_handle>               Facet_pair;
  typedef std::pair<Facet_pair,int>                          Facet_pair_and_int;
  
  typedef internal_IOP::
            Compare_handles<Polyhedron,Use_const_polyhedron> Compare_handles;

  typedef internal_IOP::
       Compare_handle_pairs<Polyhedron,Use_const_polyhedron> Compare_handle_pairs;
  

  typedef std::map<Facet_pair_and_int,                                           //we use Facet_pair_and_int and not Facet_pair to handle coplanar case.
                   std::set<int>,Compare_handle_pairs>       Facets_to_nodes_map;//Indeed the boundary of the intersection of two coplanar triangles may contain several segments.
  typedef std::set<Facet_pair,Compare_handle_pairs>          Coplanar_facets_set;//any insertion should be done with make_sorted_pair_of_facets
  typedef typename Kernel::Point_3                           Node;
  typedef internal_IOP::Triangle_segment_intersection_point
            <Polyhedron,Kernel,Node_storage_type>            Node_vector;

  typedef typename internal_IOP::
    Intersection_types<Polyhedron,Use_const_polyhedron>
      ::Intersection_result                                  Intersection_result;

  typedef std::set<Facet_handle,Compare_handles>             Facet_set;

  typedef std::map
            <Halfedge_handle,Facet_set,Compare_handles>      Edge_to_intersected_facets;
  #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
  typedef std::set<Facet_pair,Compare_handle_pairs>          Coplanar_duplicated_intersection_set;
  #endif
//helper functions
  static inline Facet_pair 
  make_sorted_pair_of_facets(Facet_handle fh1,Facet_handle fh2) {
    const Facet* f1=&(*fh1);
    const Facet* f2=&(*fh2);
    if (f1 < f2)
      return Facet_pair(fh1,fh2);
    return Facet_pair(fh2,fh1);
  }

  static inline Facet_pair_and_int
  make_sorted_pair_of_facets_with_int(Facet_handle fh1,Facet_handle fh2,int i) {
    return std::make_pair(make_sorted_pair_of_facets(fh1,fh2),i);
  }

  static inline std::pair<void*,void*> make_sorted_void_pair(void* v1,void* v2){
    if (v1<v2) return std::make_pair(v1,v2);
    return std::make_pair(v2,v1);
  }

  static inline Halfedge_handle smaller_handle(Halfedge_handle h){
    if ( &(*h)<&(*h->opposite()) ) return h;
    return h->opposite();
  }
  
//member variables
  Edge_to_intersected_facets  edge_to_sfacet; //Associate a segment to a filtered set of facets that may be intersected
  Facets_to_nodes_map         f_to_node;      //Associate a pair of triangle to their intersection points
  Coplanar_facets_set         coplanar_facets;//Contains all pairs of triangular facets intersecting that are coplanar
  Node_vector                 nodes;          //Contains intersection points of polyhedra
  Node_visitor*               visitor;
  bool                        is_default_visitor; //indicates whether the visitor need to be deleted
  #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
  //this does not occur only when one extremity of an edge is inside a face.
  // The problem occur every time an edge or a part of an edge with two incident triangles
  // is on the intersection polyline, I choose to direcly filter the output by removing duplicated edges
  Coplanar_duplicated_intersection_set coplanar_duplicated_intersection;//Set containing edges that are duplicated because of edges (partially) included in a triangle
  #endif
  
//functions that should come from a traits class
  bool has_at_least_two_incident_faces(Halfedge_handle edge)
  {
    return !edge->is_border_edge();
  }
  
  template <class Output_iterator>
  void get_incident_facets(Halfedge_handle edge,Output_iterator out){
    if (!edge->is_border()) *out++=edge->facet();
    if (!edge->opposite()->is_border()) *out++=edge->opposite()->facet();
  }

  template <class Output_iterator>
  void get_incident_edges_to_vertex(Halfedge_handle edge,Output_iterator out){
    Halfedge_handle current=edge;
    do{
      *out++=current;
      current=current->next()->opposite();
    }
    while(current!=edge);
  }

//internal functions
  
  class Map_edge_facet_bbox_intersection {
    Edge_to_intersected_facets& edge_to_sfacet;
    Polyhedron_ref polyhedron_triangle;
    Polyhedron_ref polyhedron_edge;
    Node_visitor& visitor;
  public:
    Map_edge_facet_bbox_intersection(Edge_to_intersected_facets& map_,
                                     Polyhedron_ref  P,
                                     Polyhedron_ref Q,
                                     Node_visitor& visitor_)
      :edge_to_sfacet(map_),polyhedron_triangle(P),polyhedron_edge(Q),visitor(visitor_){}

    void operator()( const Box* fb, const Box* eb) const {
      Halfedge_handle fh = fb->handle();
      Halfedge_handle eh = eb->handle();
      
      typename Edge_to_intersected_facets::iterator res=
        edge_to_sfacet.insert(std::make_pair(eh,Facet_set())).first;
      res->second.insert(fh->facet());
      visitor.add_filtered_intersection(eh,fh,polyhedron_edge,polyhedron_triangle);
    }
  };

  class Map_edge_facet_bbox_intersection_extract_coplanar {
    Edge_to_intersected_facets& edge_to_sfacet;
    Coplanar_facets_set& coplanar_facets;
    Polyhedron_ref polyhedron_triangle;
    Polyhedron_ref polyhedron_edge;
    Node_visitor& visitor;
  public:
    Map_edge_facet_bbox_intersection_extract_coplanar(
      Edge_to_intersected_facets& map_,
      Coplanar_facets_set& coplanar_facets_,
      Polyhedron_ref  P,
      Polyhedron_ref Q,
      Node_visitor& visitor_)
      :edge_to_sfacet(map_),coplanar_facets(coplanar_facets_),polyhedron_triangle(P),polyhedron_edge(Q),visitor(visitor_)
    {}

    void operator()( const Box* fb, const Box* eb) const {
      Halfedge_handle fh = fb->handle();
      Halfedge_handle eh = eb->handle();
      
      
      //check if the segment intersects the plane of the facet or if it is included in the plane
      const typename Kernel::Point_3 & a = fh->vertex()->point();
      const typename Kernel::Point_3 & b = fh->next()->vertex()->point();
      const typename Kernel::Point_3 & c = fh->next()->next()->vertex()->point();
      const Orientation abcp = orientation(a,b,c,eh->vertex()->point());
      const Orientation abcq = orientation(a,b,c,eh->opposite()->vertex()->point());
      if (abcp==abcq){
        if (abcp!=COPLANAR){
//          std::cout << "rejected " << &(*fh->facet()) << "{" << &(*eh->facet()) << " " <<&(*eh->opposite()->facet()) << " "<< eh->vertex()->point() << " " << eh->opposite()->vertex()->point() << "}" <<std::endl;
          return; //no intersection
        }
        //WARNING THIS IS DONE ONLY FOR POLYHEDRON (MAX TWO INCIDENT FACETS TO EDGE)
        if (!eh->is_border() && orientation(a,b,c,eh->next()->vertex()->point())==COPLANAR){
          coplanar_facets.insert(make_sorted_pair_of_facets(eh->facet(),fh->facet()));
        }
        if (!eh->opposite()->is_border() && orientation(a,b,c,eh->opposite()->next()->vertex()->point())==COPLANAR){
          coplanar_facets.insert(make_sorted_pair_of_facets(eh->opposite()->facet(),fh->facet()));
        }
        visitor.add_filtered_intersection(eh,fh,polyhedron_edge,polyhedron_triangle);
        //in case only the edge is coplanar, the intersection points will be detected using an incident facet 
        //(see remark at the beginning of the file)
        return; 
      }
      
      typename Edge_to_intersected_facets::iterator res=
        edge_to_sfacet.insert(std::make_pair(eh,Facet_set())).first;
      res->second.insert(fh->facet());
      visitor.add_filtered_intersection(eh,fh,polyhedron_edge,polyhedron_triangle);
    }
  };
  
  // This function tests the intersection of the faces of P with the edges of Q
  void filter_intersections( Polyhedron_ref P, Polyhedron_ref Q) {
    std::vector<Box> facet_boxes, edge_boxes;
    facet_boxes.reserve( P.size_of_facets());
    for ( Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i){
        facet_boxes.push_back(
            Box( i->halfedge()->vertex()->point().bbox()
               + i->halfedge()->next()->vertex()->point().bbox()
               + i->halfedge()->next()->next()->vertex()->point().bbox(),
                 i->halfedge()));
    }
    std::vector<const Box*> facet_box_ptr;
    facet_box_ptr.reserve( P.size_of_facets());
    for ( typename std::vector<Box>::iterator j = facet_boxes.begin(); j != facet_boxes.end(); ++j){
        facet_box_ptr.push_back( &*j);
    }
      
    for ( Halfedge_iterator i = Q.halfedges_begin(); i != Q.halfedges_end(); ++i){
      if(&*i < &*(i->opposite())){
        edge_boxes.push_back(
            Box( i->vertex()->point().bbox()
                 + i->opposite()->vertex()->point().bbox(),
                 i));
      }
   }
   
    std::vector<const Box*> edge_box_ptr;
    edge_box_ptr.reserve( Q.size_of_halfedges()/2);
    for ( typename std::vector<Box>::iterator j = edge_boxes.begin(); j != edge_boxes.end(); ++j){
        edge_box_ptr.push_back( &*j);
    }

    CGAL::box_intersection_d( facet_box_ptr.begin(), facet_box_ptr.end(),
                              edge_box_ptr.begin(), edge_box_ptr.end(),
    #ifdef DO_NOT_HANDLE_COPLANAR_FACETS
                              Map_edge_facet_bbox_intersection(edge_to_sfacet,P,Q,*visitor),
    #else
                              Map_edge_facet_bbox_intersection_extract_coplanar(edge_to_sfacet,coplanar_facets,P,Q,*visitor),
    #endif
                              std::ptrdiff_t(2000)
    );
  }


  void add_intersection_point_to_facet_and_all_edge_incident_facets(Facet_handle facet,
                                                                    Halfedge_handle edge,
                                                                    int node_id)
  {
    std::vector<Facet_handle> incident_facets;
    get_incident_facets(edge,std::back_inserter(incident_facets));
    for (typename std::vector<Facet_handle>::iterator it=incident_facets.begin();
         it!=incident_facets.end();++it)
    {
      CGAL_assertion(cgal_do_intersect_debug(facet,*it));
      
      Facet_pair facet_pair = make_sorted_pair_of_facets(facet,*it);
      if ( !coplanar_facets.empty() && coplanar_facets.find(facet_pair)!=coplanar_facets.end() ) continue;
      typename Facets_to_nodes_map::iterator it_list=
        f_to_node.insert( std::make_pair( Facet_pair_and_int(facet_pair,0),std::set<int>()) ).first;
      it_list->second.insert(node_id);
    }
  }

  void cip_handle_case_edge(int node_id,
                            Facet_set* fset,
                            Halfedge_handle edge,
                            Halfedge_handle edge_intersected)
  {
    //associate the intersection point to all facets incident to the intersected edge using edge
    std::vector<Facet_handle> incident_facets;
    get_incident_facets(edge_intersected,std::back_inserter(incident_facets));
    for (typename std::vector<Facet_handle>::iterator it=incident_facets.begin();
         it!=incident_facets.end();++it)
    {
      add_intersection_point_to_facet_and_all_edge_incident_facets(*it,edge,node_id);
      if (fset!=NULL) fset->erase(*it);
    }
    incident_facets.clear();

    //associate the intersection point to all facets incident to edge using the intersected edge
    //at least one pair of facets is already handle above
    
    typename Edge_to_intersected_facets::iterator it_fset=edge_to_sfacet.find(edge_intersected);
    if (it_fset==edge_to_sfacet.end()) return;
    Facet_set& fset_bis=it_fset->second;
    get_incident_facets(edge,std::back_inserter(incident_facets));
    for (typename std::vector<Facet_handle>::iterator it=incident_facets.begin();
         it!=incident_facets.end();++it)
    {
//      add_intersection_point_to_facet_and_all_edge_incident_facets(*it,edge_intersected,node_id); //this call is not needed, already done in the first loop
      fset_bis.erase(*it);
    }
  }

  void cip_handle_case_vertex(int node_id,
                              Facet_set* fset,
                              Halfedge_handle edge,
                              Halfedge_handle vertex_intersected)
  {
    std::vector<Halfedge_handle> incident_halfedges;
    get_incident_edges_to_vertex(vertex_intersected,std::back_inserter(incident_halfedges));
    for (typename std::vector<Halfedge_handle>::iterator 
      it=incident_halfedges.begin();it!=incident_halfedges.end();++it)
    {
      cip_handle_case_edge(node_id,fset,edge,*it);
    }
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(Halfedge_handle edge,
                    Facet_handle facet,
                    const Intersection_result& inter_res,
                    Node_vector& nodes)
  {
    bool is_vertex_coplanar = CGAL::cpp0x::get<2>(inter_res);
    if (is_vertex_coplanar)
      nodes.add_new_node(edge->vertex()->point());
    else{
      bool is_opposite_vertex_coplanar = CGAL::cpp0x::get<3>(inter_res);
      if (is_opposite_vertex_coplanar)
        nodes.add_new_node(edge->opposite()->vertex()->point());
      else
        nodes.add_new_node(edge,facet);
    }
  }

  //either the exact or input point can be used to create a node
  //with this function
  template<class Point>
  void add_new_node(const Point& pt)
  {
    nodes.add_new_node(pt);
  }

  
  #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
  void check_coplanar_edge(Halfedge_handle hedge,Facet_handle facet)
  {
    const typename Kernel::Point_3& p0=facet->halfedge()->vertex()->point();
    const typename Kernel::Point_3& p1=facet->halfedge()->next()->vertex()->point();
    const typename Kernel::Point_3& p2=facet->halfedge()->opposite()->vertex()->point();
    CGAL_precondition( orientation( p0,p1,p2,hedge->vertex()->point() ) == COPLANAR );

    if ( has_at_least_two_incident_faces(hedge) &&  orientation( p0,p1,p2,hedge->opposite()->vertex()->point() ) == COPLANAR )
    {
      //In case two facets are incident along such this edge, the intersection
      //will be reported twice. We keep track of this so that at the end, we can remove one intersecting edge out of the two
      //choose the smaller of the two faces (only one need to de deleted)
      Facet_handle smaller=make_sorted_pair_of_facets(hedge->face(),hedge->opposite()->face()).first;
      coplanar_duplicated_intersection.insert(make_sorted_pair_of_facets(smaller,facet));
    }    
  }
  
  bool are_incident_facets_coplanar(Halfedge_handle hedge){
    const typename Kernel::Point_3& p0=hedge->vertex()->point();
    const typename Kernel::Point_3& p1=hedge->next()->vertex()->point();
    const typename Kernel::Point_3& p2=hedge->opposite()->vertex()->point();
    const typename Kernel::Point_3& p3=hedge->opposite()->next()->vertex()->point();
    return orientation( p0,p1,p2,p3 ) == COPLANAR;
  }

  void check_coplanar_edge(Halfedge_handle hedge,Halfedge_handle additional_edge,internal_IOP::Intersection_type type)
  {
    switch(type){
      case internal_IOP::FACET:
        check_coplanar_edge(hedge,additional_edge->face());
      break;
      
      case internal_IOP::EDGE:
        if ( !additional_edge->is_border() ){
          check_coplanar_edge(hedge,additional_edge->face());
        }
        if (!additional_edge->opposite()->is_border())
          check_coplanar_edge(hedge,additional_edge->opposite()->face());
      break;        
      case internal_IOP::VERTEX:
      {
        //consider  all incident faces
        Halfedge_handle current=additional_edge;
        do{
          if( !current->is_border() )
          check_coplanar_edge(hedge,current->face());
          current=current->next()->opposite();
        }
        while(current!=additional_edge);
      }
      break;
      
      default:
        CGAL_assertion(type==internal_IOP::COPLNR);
      break;
    }    
  }
  
  void check_coplanar_edge_old(Halfedge_handle hedge,Halfedge_handle additional_edge,internal_IOP::Intersection_type type)
  {
    switch(type){
      case internal_IOP::FACET:
        check_coplanar_edge(hedge,additional_edge->face());
      break;
      
      case internal_IOP::EDGE:
      {
        if ( !additional_edge->is_border() ){
          if (!additional_edge->opposite()->is_border()){
            if ( are_incident_facets_coplanar(additional_edge) )
            {
              Facet_handle facet=additional_edge->face();
              const typename Kernel::Point_3& p0=facet->halfedge()->vertex()->point();
              const typename Kernel::Point_3& p1=facet->halfedge()->next()->vertex()->point();
              const typename Kernel::Point_3& p2=facet->halfedge()->opposite()->vertex()->point();
              CGAL_precondition( orientation( p0,p1,p2,hedge->vertex()->point() ) == COPLANAR );

              if ( has_at_least_two_incident_faces(hedge) &&  orientation( p0,p1,p2,hedge->opposite()->vertex()->point() ) == COPLANAR )
              {
                //In case two facets are incident along a common edge of two coplanar triangles.
                //We need to remove three out of the four reported pair
                Facet_handle smaller=make_sorted_pair_of_facets(hedge->face(),hedge->opposite()->face()).first;
                coplanar_duplicated_intersection.insert(make_sorted_pair_of_facets(hedge->face(),facet));
                coplanar_duplicated_intersection.insert(make_sorted_pair_of_facets(hedge->opposite()->face(),facet));
                coplanar_duplicated_intersection.insert(make_sorted_pair_of_facets(hedge->opposite()->face(),additional_edge->opposite()->face()));
              }              
            }
            else
            {
              check_coplanar_edge(hedge,additional_edge->face());
              check_coplanar_edge(hedge,additional_edge->opposite()->face());
            }
          }
          else
            check_coplanar_edge(hedge,additional_edge->face());
        }
        else{
          CGAL_assertion(!additional_edge->opposite()->is_border());
          check_coplanar_edge(hedge,additional_edge->opposite()->face());
        }
      }
      break;        
      case internal_IOP::VERTEX:
        
      break;
      
      default:
        CGAL_assertion(type==internal_IOP::COPLNR);
      break;
    }    
  }
  
  template <class Hedge_iterator>
  void check_coplanar_edges(Hedge_iterator begin,Hedge_iterator end,Halfedge_handle additional_edge,internal_IOP::Intersection_type type)
  {
    for (Hedge_iterator it=begin;it!=end;++it)
      check_coplanar_edge(*it,additional_edge,type);
  }
  #endif
  
  void print_type_debug(internal_IOP::Intersection_type type,bool cpl,bool opp_cpl)
  {
    switch(type){
      case internal_IOP::COPLNR:
        std::cout << "COPLNR " << cpl << " " <<  opp_cpl << std::endl;
      break;
      case internal_IOP::EMPTY:
        std::cout << "EMPTY " << cpl << " " <<  opp_cpl << std::endl;
      break;
      
      case internal_IOP::FACET:
        std::cout << "FACET " << cpl << " " <<  opp_cpl << std::endl;
      break;
      
      case internal_IOP::EDGE:
        std::cout << "EDGE " << cpl << " " <<  opp_cpl << std::endl;
      break;        
      case internal_IOP::VERTEX:
        std::cout << "VERTEX " << cpl << " " <<  opp_cpl << std::endl;
      break;
    }
  }
  

  void handle_coplanar_case_VERTEX_FACET(Halfedge_handle vertex,Halfedge_handle facet,int node_id){
    visitor->new_node_added(node_id,internal_IOP::FACET,vertex,facet,true,false);
    std::vector<Halfedge_handle> all_edges;
    get_incident_edges_to_vertex(vertex,std::back_inserter(all_edges));
    typename std::vector<Halfedge_handle>::iterator it_edge=all_edges.begin();
    for (;it_edge!=all_edges.end();++it_edge){
      add_intersection_point_to_facet_and_all_edge_incident_facets(facet->facet(),*it_edge,node_id);
      typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
      if (it_ets!=edge_to_sfacet.end()) it_ets->second.erase(facet->facet());
    }
  }
  
  void handle_coplanar_case_VERTEX_EDGE(Halfedge_handle vertex,Halfedge_handle edge,int node_id){
    visitor->new_node_added(node_id,internal_IOP::VERTEX,edge,vertex,false,false);
    std::vector<Halfedge_handle> all_edges;
    get_incident_edges_to_vertex(vertex,std::back_inserter(all_edges));
    typename std::vector<Halfedge_handle>::iterator it_edge=all_edges.begin();
    for (;it_edge!=all_edges.end();++it_edge){
      typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
      Facet_set* fset = (it_ets!=edge_to_sfacet.end())?&(it_ets->second):NULL;
      cip_handle_case_edge(node_id,fset,*it_edge,edge);
    }
  }

  void handle_coplanar_case_VERTEX_VERTEX(Halfedge_handle vertex1,Halfedge_handle vertex2,int node_id){
    visitor->new_node_added(node_id,internal_IOP::VERTEX,vertex2,vertex1,true,false);
    std::vector<Halfedge_handle> all_edges;
    get_incident_edges_to_vertex(vertex1,std::back_inserter(all_edges));
    typename std::vector<Halfedge_handle>::iterator it_edge=all_edges.begin();          
    for (;it_edge!=all_edges.end();++it_edge){
      typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
      Facet_set* fset = (it_ets!=edge_to_sfacet.end())?&(it_ets->second):NULL;
      cip_handle_case_vertex(node_id,fset,*it_edge,vertex2);
    }
  }
  
  
  template<class Cpl_inter_pt,class Coplanar_node_map>
  int get_or_create_node(Cpl_inter_pt& ipt,int& current_node,Coplanar_node_map& coplanar_node_map){
    void *v1, *v2;
    switch(ipt.type_1){
      case internal_IOP::VERTEX: v1=(void*)( &(*(ipt.info_1->vertex()))     );  break;
      case internal_IOP::EDGE  : v1=(void*)( &(*smaller_handle(ipt.info_1)) );  break;
      case internal_IOP::FACET : v1=(void*)( &(*(ipt.info_1->facet()))      );  break;
      default: CGAL_assertion(!"Should not get there!");
    }

    switch(ipt.type_2){
      case internal_IOP::VERTEX: v2=(void*)( &(*(ipt.info_2->vertex()))     );  break;
      case internal_IOP::EDGE  : v2=(void*)( &(*smaller_handle(ipt.info_2)) );  break;
      case internal_IOP::FACET : v2=(void*)( &(*(ipt.info_2->facet()))      );  break;
      default: CGAL_assertion(!"Should not get there!");
    }

    std::pair<void*,void*> key=make_sorted_void_pair(v1,v2);

    std::pair<typename Coplanar_node_map::iterator,bool> res=coplanar_node_map.insert(std::make_pair(key,current_node+1));
    if (res.second){ //insert a new node
      
      if (ipt.type_1==internal_IOP::VERTEX)
        add_new_node(ipt.info_1->vertex()->point());
      else{
        if(ipt.type_2==internal_IOP::VERTEX)
          add_new_node(ipt.info_2->vertex()->point());
        else
          add_new_node(ipt.point);
      }
      return ++current_node;
    }
    return res.first->second;
  }
  
  void compute_intersection_of_coplanar_facets(int& current_node){
    typedef std::map<std::pair<void*,void*>,int> Coplanar_node_map;
    Coplanar_node_map coplanar_node_map;
    
    for (typename Coplanar_facets_set::iterator it=coplanar_facets.begin();it!=coplanar_facets.end();++it){
      Facet_handle f1=it->first;
      Facet_handle f2=it->second;
      typedef internal_IOP::Intersection_point_with_info<Kernel,Halfedge_handle> Cpl_inter_pt;
      std::list<Cpl_inter_pt> inter_pts;
      internal_IOP::intersection_coplanar_facets<Kernel>(f1->halfedge(),f2->halfedge(),inter_pts);
//      std::cout << "found " << inter_pts.size() << " inter pts: "; 
      std::size_t nb_pts=inter_pts.size();
      std::vector<int> cpln_nodes; cpln_nodes.reserve(nb_pts);
      for (typename std::list<Cpl_inter_pt>::iterator iti=inter_pts.begin();iti!=inter_pts.end();++iti){
        #ifdef CGAL_COREFINEMENT_DEBUG
        //iti->print_debug();
        #endif
        CGAL_assertion(iti->is_valid());
        int node_id=get_or_create_node(*iti,current_node,coplanar_node_map);
        cpln_nodes.push_back(node_id);

        switch(iti->type_1){
          case internal_IOP::VERTEX:
          {
            switch(iti->type_2){
              case internal_IOP::VERTEX: handle_coplanar_case_VERTEX_VERTEX(iti->info_1,iti->info_2,node_id); break;
              case internal_IOP::EDGE  : handle_coplanar_case_VERTEX_EDGE(iti->info_1,iti->info_2,node_id);   break;
              case internal_IOP::FACET : handle_coplanar_case_VERTEX_FACET(iti->info_1,iti->info_2,node_id);  break;
              default: CGAL_assertion(!"Should not get there!");
            }
          }
          break;
          case internal_IOP::EDGE:{
            switch(iti->type_2){
              case internal_IOP::VERTEX:handle_coplanar_case_VERTEX_EDGE(iti->info_2,iti->info_1,node_id);break;
              case internal_IOP::EDGE:
              {
                visitor->new_node_added(node_id,internal_IOP::EDGE,iti->info_1,iti->info_2,false,false);
                typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(iti->info_1);
                Facet_set* fset = (it_ets!=edge_to_sfacet.end())?&(it_ets->second):NULL;
                cip_handle_case_edge(node_id,fset,iti->info_1,iti->info_2);                
              }
              break;
              default: CGAL_assertion(!"Should not get there!");
            }
          }            
          break;
          
          case internal_IOP::FACET:
          {
            CGAL_assertion(iti->type_2==internal_IOP::VERTEX);
            handle_coplanar_case_VERTEX_FACET(iti->info_2,iti->info_1,node_id);
          }
          break;
          
          default: CGAL_assertion(!"Should not get there!");
        }
      }
      switch (nb_pts){
        case 0: break;
        case 1:
        {
          typename Facets_to_nodes_map::iterator it_list=  
            f_to_node.insert( std::make_pair( Facet_pair_and_int(*it,1),std::set<int>()) ).first;
          it_list->second.insert(cpln_nodes[0]);
        }
        break;
        default:
        {
          int i=0;
          std::size_t stop=nb_pts + (nb_pts<3?-1:0);
          for (std::size_t k=0;k<stop;++k){
            typename Facets_to_nodes_map::iterator it_list=
              f_to_node.insert( std::make_pair( Facet_pair_and_int(*it,++i),std::set<int>()) ).first;
            it_list->second.insert( cpln_nodes[k] );
            it_list->second.insert( cpln_nodes[(k+1)%nb_pts] );
          }
        }
      }
//      std::cout << std::endl;
    }
  }
  
  void compute_intersection_points(int& current_node){
    for(typename Edge_to_intersected_facets::iterator it=edge_to_sfacet.begin();it!=edge_to_sfacet.end();++it){
      Halfedge_handle edge=it->first;
      Facet_set& fset=it->second;
      while (!fset.empty()){
        Facet_handle facet=*fset.begin();
        
        Intersection_result res=internal_IOP::do_intersect<Polyhedron,Kernel,Use_const_polyhedron>(edge,facet);
        internal_IOP::Intersection_type type=CGAL::cpp0x::get<0>(res);
        
        //handle degenerate case: one extremity of edge below to facet
        std::vector<Halfedge_handle> all_edges;
        if ( CGAL::cpp0x::get<2>(res) )
          get_incident_edges_to_vertex(edge,std::back_inserter(all_edges));
        else{
          if ( CGAL::cpp0x::get<3>(res) )
            get_incident_edges_to_vertex(edge->opposite(),std::back_inserter(all_edges));
          else
            all_edges.push_back(edge);
        }
        
        CGAL_precondition(*all_edges.begin()==edge || *all_edges.begin()==edge->opposite());
//        print_type_debug(type,CGAL::cpp0x::get<2>(res),CGAL::cpp0x::get<3>(res));
       
        #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
        check_coplanar_edges(boost::next(all_edges.begin()),all_edges.end(),CGAL::cpp0x::get<1>(res),type);
        #endif
        
        typename std::vector<Halfedge_handle>::iterator it_edge=all_edges.begin();
        switch(type){
          case internal_IOP::COPLNR:
            #ifndef DO_NOT_HANDLE_COPLANAR_FACETS
            assert(!"COPLNR : this point should never be reached!");
            #else
            //nothing need to be done, cf. comments at the beginning of the file
            #endif
          break;
          case internal_IOP::EMPTY:
            fset.erase(fset.begin());
            CGAL_assertion(!cgal_do_intersect_debug(edge,facet));
          break;
          
          case internal_IOP::FACET:
          {
            CGAL_assertion(cgal_do_intersect_debug(edge,facet));
            CGAL_assertion(facet==CGAL::cpp0x::get<1>(res)->face());
            
            int node_id=++current_node;
            add_new_node(edge,facet,res,nodes);
            visitor->new_node_added(node_id,internal_IOP::FACET,edge,facet->halfedge(),CGAL::cpp0x::get<2>(res),CGAL::cpp0x::get<3>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              add_intersection_point_to_facet_and_all_edge_incident_facets(facet,*it_edge,node_id);
              //erase facet from the list to test intersection with it_edge
              if ( it_edge==all_edges.begin() )
                fset.erase(fset.begin());
              else
              {
                typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
                if(it_ets!=edge_to_sfacet.end()) it_ets->second.erase(facet);
              }
            }
          }
          break;
          
          case internal_IOP::EDGE:
          {
            CGAL_assertion(cgal_do_intersect_debug(edge,facet));
            int node_id=++current_node;
            add_new_node(edge,facet,res,nodes);
            Halfedge_handle edge_intersected=CGAL::cpp0x::get<1>(res);
            visitor->new_node_added(node_id,internal_IOP::EDGE,edge,edge_intersected,CGAL::cpp0x::get<2>(res),CGAL::cpp0x::get<3>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
                Facet_set* fset_bis = (it_ets!=edge_to_sfacet.end())?&(it_ets->second):NULL;
                cip_handle_case_edge(node_id,fset_bis,*it_edge,edge_intersected);
              }
              else
                cip_handle_case_edge(node_id,&fset,*it_edge,edge_intersected);
            }
          }
          break;        
          
          case internal_IOP::VERTEX:
          {
            CGAL_assertion(cgal_do_intersect_debug(edge,facet));
            int node_id=++current_node;
            Halfedge_handle vertex_intersected=CGAL::cpp0x::get<1>(res);
            add_new_node(vertex_intersected->vertex()->point()); //we use the original vertex to create the node
            //before it was internal_IOP::FACET but do not remember why, probably a bug...
            visitor->new_node_added(node_id,internal_IOP::VERTEX,edge,vertex_intersected,CGAL::cpp0x::get<2>(res),CGAL::cpp0x::get<3>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_intersected_facets::iterator it_ets=edge_to_sfacet.find(*it_edge);
                Facet_set* fset_bis = (it_ets!=edge_to_sfacet.end())?&(it_ets->second):NULL;                
                cip_handle_case_vertex(node_id,fset_bis,*it_edge,vertex_intersected);
              }
              else
                cip_handle_case_vertex(node_id,&fset,*it_edge,vertex_intersected);
            }
          }
          break;          

        }
      }
    }
    CGAL_assertion(nodes.size()==unsigned(current_node+1));
  }
  
  #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
  void remove_duplicated_intersecting_edges()
  {
    for (typename Coplanar_duplicated_intersection_set::iterator 
        it=coplanar_duplicated_intersection.begin();
        it!=coplanar_duplicated_intersection.end();
        ++it)
    {
      typename Facets_to_nodes_map::iterator res=f_to_node.find(*it);
      //CGAL_assertion(res!=f_to_node.end());
      //we cannot use a precondition here: in some cases the coplanarity test is not sufficient.
      //if on surface1 we have to coplanar triangle incident along edge e. Then in surface2, take a triangle
      //with one edge inside one triangle of surface1 such that one vertex in on e. Then the collinearity test
      //return true for both faces but only one implies a duplicated edge
      if(res!=f_to_node.end())
        f_to_node.erase(res);
    }
  }
  #else
  void remove_duplicated_intersecting_edges()
  {
    std::set< std::pair<int,int> > already_seen;
    std::list<typename Facets_to_nodes_map::iterator> to_erase;
    for (typename Facets_to_nodes_map::iterator 
          it=f_to_node.begin();
          it!=f_to_node.end();
          ++it
    )
    {
      if (it->second.size()==1) continue;
      CGAL_precondition(it->second.size()==2);
      //it->second is a set so int are already sorted
      bool is_new=already_seen.insert(  std::make_pair(
                                       *(it->second.begin()),
                                       *boost::next(it->second.begin()) )
      ).second;
      
      if (!is_new)
        to_erase.push_back(it);
    }
    
    for ( typename std::list<typename Facets_to_nodes_map::iterator>::iterator
          it=to_erase.begin();it!=to_erase.end();++it
    )
      f_to_node.erase(*it);
  }
  #endif


  struct Graph_node{
    std::set<int> neighbors;
    unsigned size;
    
    Graph_node():size(0){}
    
    void insert(int i){
      ++size;
      CGAL_assertion(neighbors.find(i)==neighbors.end());
      neighbors.insert(i);
    }
    
    void erase(int i){
      CGAL_assertion(neighbors.find(i)!=neighbors.end());
      neighbors.erase(i);
    }
    void make_terminal() {size=45;}
    bool is_terminal()const {return size!=2;}
    bool empty() const {return neighbors.empty();}
    int top() const {return *neighbors.begin();}
    void pop() {
      CGAL_assertion(!neighbors.empty());
      neighbors.erase(neighbors.begin());
    }
  };


  template <class Output_iterator>
  void construct_polylines(Node_vector& nodes,Output_iterator out){
    typedef std::map<int,Graph_node> Graph;
    Graph graph;
    
    //counts the number of time each node has been seen
    std::size_t nb_nodes=nodes.size();
    std::vector<int> node_mult(nb_nodes,0);
    bool isolated_point_seen=false;
    for (typename Facets_to_nodes_map::iterator it=f_to_node.begin();it!=f_to_node.end();++it){
      const std::set<int>& segment=it->second;
      CGAL_assertion(segment.size()==2 || segment.size()==1);
      if (segment.size()==2){
        int i=*segment.begin();
        int j=*boost::next(segment.begin());
        typename Graph::iterator ins_res=graph.insert(std::make_pair(i,Graph_node())).first;
        ins_res->second.insert(j);
        ins_res=graph.insert(std::make_pair(j,Graph_node())).first;
        ins_res->second.insert(i);
        ++(node_mult[i]);
        ++(node_mult[j]);
      }
      else{
        CGAL_assertion(segment.size()==1);
        isolated_point_seen=true;
      }
    }
    
    //add isolated points
    if (isolated_point_seen){
      for (unsigned index=0;index<nb_nodes;++index)
        if (node_mult[index]==0){
          *out++=std::vector<typename Kernel::Point_3>(1,nodes[index]);
          visitor->start_new_polyline(index,index);
        }
    }
    
    //visitor call
    visitor->annotate_graph(graph.begin(),graph.end());
    
    bool only_cycle=false;
    while (!graph.empty()){
      typename Graph::iterator it=graph.begin();
      for (;!only_cycle && it!=graph.end();++it){
        if (it->second.is_terminal()) break;
      }
      
      std::vector<typename Kernel::Point_3> polyline;
      
      if(!only_cycle && it!=graph.end()){
        //this is a polyline
        int i=it->first;
        int j=it->second.top();
        visitor->start_new_polyline(i,j);
        CGAL_assertion(i!=j);
        it->second.pop();
        if (it->second.empty())
          graph.erase(it);
        polyline.push_back(nodes[i]);
        visitor->add_node_to_polyline(i);
        while(true){
          it=graph.find(j);
          CGAL_assertion(it!=graph.end());
          it->second.erase(i);
          i=j;
          polyline.push_back(nodes[i]);
          visitor->add_node_to_polyline(i);
          if (it->second.empty()){
            graph.erase(it);
            break;
          }
          if (it->second.is_terminal()) break;
          j=it->second.top();
          it->second.pop();
          if (it->second.empty())
            graph.erase(it);
        }
        *out++=polyline;
      }
      else{
        //it remains only cycles
        only_cycle=true;
        it=graph.begin();
        int i=it->first;
        int j=it->second.top();
        visitor->start_new_polyline(i,j);
        graph.erase(it);
        polyline.push_back(nodes[i]);
        visitor->add_node_to_polyline(i);
        int first=i;
        do{
          it=graph.find(j);
          it->second.erase(i);
          i=j;
          polyline.push_back(nodes[i]);
          visitor->add_node_to_polyline(i);
          j=it->second.top();
          graph.erase(it);
        }while(j!=first);
        polyline.push_back(nodes[j]);// we duplicate first point for cycles
        visitor->add_node_to_polyline(j);
        *out++=polyline;      
      }
    }
  }
  
  int get_other_int(const std::set<int>& s,int i) const {
    if (*s.begin()!=i) return *s.begin();
    return *boost::next(s.begin());
  }
  
  template <class T,class Dispatch_out_it>
  bool is_grabbed(const Dispatch_out_it&) const{
    return CGAL::Is_in_tuple<T,typename Dispatch_out_it::Value_type_tuple>::value;
  }
  
  
  template <class OutputIterator>
  struct Is_dispatch_based_ouput_iterator{
    typedef  boost::false_type type;
  };

  template <template<class V_,class O_> class Dispatch_based_output_it,class V,class O>
  struct Is_dispatch_based_ouput_iterator<Dispatch_based_output_it<V,O> >{
    typedef typename boost::is_base_of< Dispatch_output_iterator<V,O>,
                                        Dispatch_based_output_it<V,O>  >::type type;
  };
  
  template <class Output_iterator>
  inline void construct_polylines_with_info(Node_vector& nodes,Output_iterator out){
    return construct_polylines_with_info(nodes,out,typename Is_dispatch_based_ouput_iterator<Output_iterator>::type());
  }
    
  template <class Output_iterator>
  void construct_polylines_with_info(Node_vector& nodes,Output_iterator out,boost::false_type){
    construct_polylines_with_info(nodes,
                                  dispatch_or_drop_output<std::vector<typename Kernel::Point_3> >(out),boost::true_type());
  }
  
  template <template<class V_,class O_> class Dispatch_based_output_it,class V,class O>
  void construct_polylines_with_info(Node_vector& nodes,
                                     Dispatch_based_output_it<V,O> out,boost::true_type)
  {
    typedef typename Facets_to_nodes_map::value_type Edge;
    typedef std::list<Facet_pair> Polyline_info;
    
    
    std::size_t nb_nodes=nodes.size();
    std::vector<int> node_mult(nb_nodes,0);
    std::vector<bool> terminal_bools(nb_nodes,false);
    std::vector< std::set<Edge*> > connections(nb_nodes);
    // --counts the number of time each node has been seen
    // --associate to each node its incident edges. An edge = a pair of Facet_handle+2 indices of intersection points
    for (typename Facets_to_nodes_map::iterator it=f_to_node.begin();it!=f_to_node.end();++it){
      const std::set<int>& segment=it->second;
      CGAL_assertion(segment.size()==2 || segment.size()==1);
      if (segment.size()==2){
        int i=*segment.begin();
        int j=*boost::next(segment.begin());
        connections[i].insert(&(*it));
        connections[j].insert(&(*it));
        ++(node_mult[i]);
        ++(node_mult[j]);
      }
    }

    //detect terminal nodes and isolated nodes
    for (unsigned k=0;k<nb_nodes;++k){
      if (node_mult[k]!=2) terminal_bools[k]=true;
      if (node_mult[k]==0){
        *out++=std::vector<typename Kernel::Point_3>(1,nodes[k]);
        if ( is_grabbed<std::vector<Facet_pair> >(out))
          *out++=std::vector<Facet_pair>();
      }
    }
   
    //visitor call
    visitor->update_terminal_nodes(terminal_bools);
    
    //We start from a node N and recursively walk one edge to find other
    // node. If this is a cycle we stop when we are back at the node N;
    //otherwise we stop at a terminal node and restart a walk from N
    //following another edge (if N is not terminal) until finding a terminal
    //node. With this we can associate to each edge the pair of facet_handle
    //intersecting that define this edge.
    unsigned current_node=0;
    while (current_node!=nb_nodes){
      if (connections[current_node].empty()){
        ++current_node;
        continue;
      }
      
      Edge* edge=*connections[current_node].begin();
      connections[current_node].erase(connections[current_node].begin());
      
      Polyline_info polyline_info;
      std::list<typename Kernel::Point_3> polyline_embedding;
      
      if ( is_grabbed<std::vector<Facet_pair> >(out))
        polyline_info.push_back(edge->first.first);
      polyline_embedding.push_back(nodes[current_node]);
      
      unsigned i=current_node;
      unsigned start_node=current_node;
      bool reverse=false;
      while (true){
        i=get_other_int(edge->second,i);
        connections[i].erase(edge);
        
        if (reverse) polyline_embedding.push_front(nodes[i]);
        else         polyline_embedding.push_back(nodes[i]);
        
        if (i==start_node) break;
        if (terminal_bools[i]){
          if (reverse || terminal_bools[current_node]) break;
          reverse=true;
          i=current_node;
          if ( connections[i].empty() ) break;
        }
        
        edge=*connections[i].begin();
        connections[i].erase(connections[i].begin());
        
        if ( is_grabbed<std::vector<Facet_pair> >(out)){
          if (reverse) polyline_info.push_front(edge->first.first);
          else         polyline_info.push_back(edge->first.first);
        }
          
      }
      
      *out++=std::vector<typename Kernel::Point_3>(polyline_embedding.begin(),polyline_embedding.end());
      if ( is_grabbed<std::vector<Facet_pair> >(out)){
        CGAL_assertion(polyline_embedding.size()==polyline_info.size()+1);
        *out++=std::vector<Facet_pair>(polyline_info.begin(),polyline_info.end());
      }
    }
  }
  
//debug functions
  
  bool cgal_do_intersect_debug(Halfedge_handle eh,Facet_handle fh){
    Triangle t( fh->halfedge()->vertex()->point(),
                fh->halfedge()->next()->vertex()->point(),
                fh->halfedge()->next()->next()->vertex()->point());

    Segment s( eh->vertex()->point(),
               eh->opposite()->vertex()->point());

    return CGAL::do_intersect( s, t);
  }

  bool cgal_do_intersect_debug(Facet_handle fh1,Facet_handle fh2){
    Triangle t1( fh1->halfedge()->vertex()->point(),
                 fh1->halfedge()->next()->vertex()->point(),
                 fh1->halfedge()->next()->next()->vertex()->point());
    Triangle t2( fh2->halfedge()->vertex()->point(),
                 fh2->halfedge()->next()->vertex()->point(),
                 fh2->halfedge()->next()->next()->vertex()->point());


    return CGAL::do_intersect( t1, t2);
  }
  
  void print_f_to_node_debug(){
    std::cout << "print_f_to_node_debug " << &f_to_node << std::endl;
    for (typename Facets_to_nodes_map::iterator it=f_to_node.begin();it!=f_to_node.end();++it){
      std::cout << &(*(it->first.first.first)) << " " << &(*(it->first.first.second)) << " " << it->first.second << " -> {";
      std::copy(it->second.begin(),it->second.end(),std::ostream_iterator<int>(std::cout,","));
      std::cout << "}" <<std::endl;
    }
  }
  
  void print_graph_debug(const std::map<int,Graph_node>& graph){
    for (typename std::map<int,Graph_node>::const_iterator it=graph.begin();it!=graph.end();++it){
      std::cout << it->first << " -> {";
      std::copy(it->second.neighbors.begin(),it->second.neighbors.end(),std::ostream_iterator<int>(std::cout,","));
      std::cout << "}" <<std::endl;      
    }
  }
  
  void print_edge_to_sfacet_debug(){
    std::cout << "Potential intersection "<< edge_to_sfacet.size() << std::endl;
    for(typename Edge_to_intersected_facets::iterator it=edge_to_sfacet.begin();it!=edge_to_sfacet.end();++it){
      Facet_set& fset=it->second;
      std::cout << &fset << " fset size " << fset.size() << std::endl;
    }
  }
  
  
  template <class OutputIterator>
  OutputIterator main_run(OutputIterator out,bool build_polylines=true){
    // std::cout << edge_to_sfacet.size() << std::endl;
    int current_node=-1;
    
    //print_edge_to_sfacet_debug();
    #ifndef DO_NOT_HANDLE_COPLANAR_FACETS
    //first handle coplanar triangles
    compute_intersection_of_coplanar_facets(current_node);
    visitor->set_number_of_intersection_points_from_coplanar_facets(current_node+1);
    #endif
    //print_edge_to_sfacet_debug();
    //compute intersection points of segments and triangles.
    //build node of the graph
    //build connectivity info
    compute_intersection_points(current_node);
    
    if (!build_polylines){
      visitor->finalize(nodes);
      return out;
    }
    //remove duplicated intersecting edges:      
    //  In case two facets are incident along such an edge coplanar in a facet of another polyhedron (and one extremity inside the facet), the intersection
    //  will be reported twice. We kept track (check_coplanar_edge(s)) of this so that, we can remove one intersecting edge out of the two
    //print_f_to_node_debug();
    remove_duplicated_intersecting_edges();
    
    //std::cout << "f_to_node "<< f_to_node.size() << std::endl;
    //print_f_to_node_debug();
    //collect connectivity infos and create polylines
    if ( Node_visitor::do_need_vertex_graph )
      construct_polylines(nodes,out); //using the graph approach (at some point we know all connections between intersection points)
    else
      construct_polylines_with_info(nodes,out); //direct construction by propagation
    
    visitor->finalize(nodes);    
    
    return out;
  }
  
public:
  Intersection_of_Polyhedra_3():visitor(new Node_visitor()),is_default_visitor(true){}
  Intersection_of_Polyhedra_3(Node_visitor& v):visitor(&v),is_default_visitor(false){}
  ~Intersection_of_Polyhedra_3(){if (is_default_visitor) delete visitor;}
  //pairwise intersection between all elements in the range
  template <class InputIterator, class OutputIterator>
  OutputIterator
  operator()(InputIterator begin, InputIterator end, OutputIterator out) {
    for(InputIterator it1=begin;it1!=end;++it1){
      CGAL_precondition( it1->is_pure_triangle() );
      Polyhedron_ref P=*it1;
      visitor->new_input_polyhedron(P);
      for(InputIterator it2=boost::next(it1);it2!=end;++it2){  
        CGAL_precondition( it2->is_pure_triangle() );
        Polyhedron_ref Q=*it2;
        filter_intersections(P, Q);
        filter_intersections(Q, P);
      }
    }
    
    return main_run(out);
  }
  
  //pairwise intersection between all elements in the range
  //(pointers version)
  template <class InputIterator, class OutputIterator>
  OutputIterator
  operator()(InputIterator begin, InputIterator end, OutputIterator out, int) {
    for(InputIterator it1=begin;it1!=end;++it1){
      CGAL_precondition( (*it1)->is_pure_triangle() );
      Polyhedron_ref P=**it1;
      visitor->new_input_polyhedron(P);
      for(InputIterator it2=boost::next(it1);it2!=end;++it2){  
        CGAL_precondition( (*it2)->is_pure_triangle() );
        Polyhedron_ref Q=**it2;
        filter_intersections(P, Q);
        filter_intersections(Q, P);
      }
    }
    
    return main_run(out);
  }
  
  //intersection between P and each element in the range
  template <class InputIterator, class OutputIterator>
  OutputIterator
  operator()( Polyhedron_ref P, InputIterator begin, InputIterator end, OutputIterator out) {
    CGAL_precondition( P.is_pure_triangle() );
    visitor->new_input_polyhedron(P);
    for(InputIterator it=begin;it!=end;++it){
      CGAL_precondition( it->is_pure_triangle() );
      Polyhedron_ref Q=*it;
      visitor->new_input_polyhedron(Q);
      filter_intersections(P, Q);
      filter_intersections(Q, P);
    }
    return main_run(out);
  }
  
  //intersection between P and each element in the range
  //(pointers version)
  template <class InputIterator, class OutputIterator>
  OutputIterator
  operator()(Polyhedron_ref P, InputIterator begin, InputIterator end, OutputIterator out, int) {
    CGAL_precondition( P.is_pure_triangle() );
    visitor->new_input_polyhedron(P);
    for(InputIterator it=begin;it!=end;++it){
      CGAL_precondition( (*it)->is_pure_triangle() );
      Polyhedron_ref Q=**it;
      visitor->new_input_polyhedron(Q);
      filter_intersections(P, Q);
      filter_intersections(Q, P);
    }
    return main_run(out);
  }
  
  //intersection between P and Q
  template <class OutputIterator>
  OutputIterator
  operator()(Polyhedron_ref P, Polyhedron_ref Q, OutputIterator out) {
    visitor->new_input_polyhedron(P);
    visitor->new_input_polyhedron(Q);
    filter_intersections(P, Q);
    filter_intersections(Q, P);
    return main_run(out);
  }

  //intersection between P and Q, only visitor called not polyline is constructed
  void operator()(Polyhedron_ref P, Polyhedron_ref Q) {
    visitor->new_input_polyhedron(P);
    visitor->new_input_polyhedron(Q);
    filter_intersections(P, Q);
    filter_intersections(Q, P);
    main_run(Emptyset_iterator(),false);
  }
};

template <typename Polyhedron, typename OutputIterator>
OutputIterator
intersection_Polyhedron_3_Polyhedron_3(const Polyhedron& P, const Polyhedron& Q, OutputIterator out)
{
  return Intersection_of_Polyhedra_3<Polyhedron>()(P,Q,out);
}

}// namespace CGAL

#endif //CGAL_INTERSECTION_OF_POLYHEDRA_3_H
