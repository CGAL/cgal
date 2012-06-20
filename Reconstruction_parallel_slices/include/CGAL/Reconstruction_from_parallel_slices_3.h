// Copyright (c) 2011  GeometryFactory (France).
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
// $URL$
// $Id$
//
// Author(s)     : Andreas Fabri
//                 Sebastien Loriot
// 

#ifndef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3
#define CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3


//Notes:
// *if DO_NOT_FILTER_NOTCHES is not defined, a graph of Voronoi vertices to be inserted is constructed (graph_G)
//  and some parts are removed to avoid so called notches. The component to be removed are connected to the infinite cell.
//  The filter consider a quotient of the max length of the bbox of a component of the graph vs min length of the cdt where
//  points should be inserted
//
// *if DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS is not defined, we consider the Voronoi vertices to be inserted and we check those
//  that have a dual facet which dual center is outside the cdt (where dual points should be inserted). We can define a segment (or a ray
//  if the facet is infinite)  which is part of the medial axis. We intersect this segment/ray with the cdt while we do not reach the boundary.
//  the edge intersection points are inserted into the cdt. If DO_NOT_FILTER_NOTCHES is not defined, only points on a constrained edges are inserted 



//#define DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
#define CGAL_NO_EDGE_EDGE_EXTRA_REMOVAL
//#define DO_NOT_FILTER_NOTCHES
//#define DO_NOT_HANDLE_NON_MANIFOLD_POINT
#define CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS

#if defined(CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS) && defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT)
    #error this is an invalid configuration
#endif

#if defined(DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS) && !defined(DO_NOT_FILTER_NOTCHES)
    #error filtering notches needs contour intersection with media axis
#endif


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
//#include <CGAL/nearest_vertex.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <bitset>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
#include <boost/function_output_iterator.hpp>
#endif

#ifndef CGAL_NO_EDGE_EDGE_EXTRA_REMOVAL
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#endif

namespace CGAL{

template <class t_Point_3>
class Slice_writer_into_file{
  int m_last_point_index;
  std::list<t_Point_3> m_surface_points;
  std::list<cpp0x::tuple<int,int,int> > m_surface_indices;
  std::string m_fname;
public:
  typedef t_Point_3 Point_3;

  Slice_writer_into_file(){}
  Slice_writer_into_file(const std::string& fname): m_last_point_index(-1),m_fname(fname){}
  void surface_point_push_back(const Point_3& p){
    m_surface_points.push_back(p);
    ++m_last_point_index;
  }
  
  void surface_indices_push_back(const cpp0x::tuple<int,int,int>& indices){
    m_surface_indices.push_back(indices);
  }

  int last_point_index() {return m_last_point_index;};
  
  void finalize(){
    std::ofstream output(m_fname.c_str() );
    output << "OFF " << m_last_point_index+1 << " " << m_surface_indices.size() << " 0\n";
    std::copy(m_surface_points.begin(),m_surface_points.end(),std::ostream_iterator<Point_3>(output,"\n"));
    for( typename std::list<cpp0x::tuple<int,int,int> >::const_iterator it=m_surface_indices.begin();
          it!=m_surface_indices.end();++it)
      output << "3 " << cpp0x::get<0>(*it) << " " << cpp0x::get<1>(*it) << " " << cpp0x::get<2>(*it) << "\n";
  }
  
  void finalize_layer(){}
};

template <class Polyhedron,class Kernel>
struct Incremental_slice_writer_into_polyhedron{
  typedef typename Kernel::Point_3 Point_3;
private:
  template <class HDS>
  class Import_modifier : public CGAL::Modifier_base<HDS> {
    typedef CGAL::Polyhedron_incremental_builder_3<HDS> IBuilder;
    

    unsigned int m_npts, m_nfcs;
    const std::list<Point_3>& m_surface_points;
    const std::list<cpp0x::tuple<int,int,int> >& m_surface_indices;    
    
  public:
    Import_modifier(unsigned int npts, unsigned int nfcs,const std::list<Point_3>& surface_points,const std::list<cpp0x::tuple<int,int,int> >& surface_indices)
      :m_npts(npts), m_nfcs(nfcs), m_surface_points(surface_points), m_surface_indices(surface_indices){}
    void operator()( HDS& hds) {
      IBuilder B( hds,true );
      B.begin_surface( m_npts, m_nfcs,0,IBuilder::ABSOLUTE_INDEXING);

      for (typename std::list<Point_3>::const_iterator it=m_surface_points.begin(),
                                                       it_end=m_surface_points.end();it!=it_end;++it)
        B.add_vertex(*it);
      
      for(typename std::list<cpp0x::tuple<int,int,int> >::const_iterator it=m_surface_indices.begin(),
                                                                         it_end=m_surface_indices.end();it!=it_end; ++it)
      {
        B.begin_facet();
        B.add_vertex_to_facet( cpp0x::get<0>(*it) );
        B.add_vertex_to_facet( cpp0x::get<1>(*it) );
        B.add_vertex_to_facet( cpp0x::get<2>(*it) );
        B.end_facet();
      }
      B.end_surface();
    }
  };
  
  
  int m_last_point_index;
  unsigned int m_npts, m_nfcs;
  std::list<Point_3> m_surface_points;
  std::list<cpp0x::tuple<int,int,int> > m_surface_indices;
  Polyhedron* m_poly_ptr;
public:
  Incremental_slice_writer_into_polyhedron(){}
  Incremental_slice_writer_into_polyhedron(Polyhedron& poly): m_last_point_index(-1),m_npts(0), m_nfcs(0), m_poly_ptr(&poly){}
    
  void surface_point_push_back(const Point_3& p){
    ++m_npts;
    m_surface_points.push_back(p);
    ++m_last_point_index;
  }
  
  void surface_indices_push_back(const cpp0x::tuple<int,int,int>& indices){
    ++m_nfcs;
    m_surface_indices.push_back(indices);
  }

  int last_point_index() {return m_last_point_index;};
  
  void finalize_layer() {
    Import_modifier<typename Polyhedron::HalfedgeDS> modif(m_npts,m_nfcs,m_surface_points,m_surface_indices);
    m_poly_ptr->delegate(modif);
    m_surface_indices.clear();
    m_surface_points.clear();
  }
  
  void finalize(){}
};


template<class T1,class T2>
std::pair<T1,T2>
make_sorted_pair(T1 va,T2 vb) {
  if (va<vb) return std::make_pair(va,vb);
  return std::make_pair(vb,va);
}

template <class Slice_writer>
class Reconstruction_from_parallel_slices_3{
  #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  //state whether a non-constrained edge made of two contour points that is non-manifold should be split to become manifold
  static const bool split_non_manifold_incontour_edges=true;
  #endif
  //the minimum distance to insert a point on a constraint or inside a face
  static const double m_min_point_squared_distance=1;
  
  #ifndef DO_NOT_FILTER_NOTCHES
  //filtering of small medial axis part that can induce bad artifact on the reconstructed surface
  static const double m_bbox_ratio_ma_filtering=10;
  #endif
  
  #ifndef CGAL_NO_EDGE_EDGE_EXTRA_REMOVAL  
  //these angle bounds are used to filter out T22 cells attached by two facets to the volume to
  //remove sharp features on the surface when we expect it to be smooth
  static const double max_sharp_dihedral_angle=90;
  static const double min_smooth_dihedral_angle=90;
  #endif
  
  struct VertexInfo2
  {
    VertexInfo2():on_contour(true){} //set to true by default for make_conforming_Gabriel_2
    int index;
    double z;
    int poly_index;
    bool on_contour;
  };

#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
  template <class Face_handle>
  struct T_FaceInfo2
  {
    #warning see whether storing the T31 associated makes sense (all planar facet maps and setting non-manifoldness can be direct)
    //PRO:
    //mark already printed facets (bitset) -> remove previous_layer_printed_planar_facets
    //mark printed nm_edge (bitset)
    //previous_bottom_incontour_nm_vertices (info in the vertex!)
    T_FaceInfo2():other_face(NULL),other_index(-1){}
    int nesting_level;

    bool in_domain() 
    { 
      return nesting_level%2 == 1;
    }
    //face of the other triangulation containing the projected circumcenter of
    //the face
    Face_handle other_face;
    //the index is != -1 if the point is not on an edge
    int other_index;
    bool dual_point_on_constraint;
  };
  
  template < typename GT,
             typename Fb = CGAL::Triangulation_face_base_2<GT> >
  class Triangulation_face_base_with_FaceInfo2
    : public Fb
  {
  public:
    typedef typename Fb::Vertex_handle                   Vertex_handle;
    typedef typename Fb::Face_handle                     Face_handle;
    typedef T_FaceInfo2<Face_handle>                       Info;

    template < typename TDS2 >
    struct Rebind_TDS {
      typedef typename Fb::template Rebind_TDS<TDS2>::Other       Fb2;
      typedef Triangulation_face_base_with_FaceInfo2<GT, Fb2>  Other;
    };

    Triangulation_face_base_with_FaceInfo2()
      : Fb() {}

    Triangulation_face_base_with_FaceInfo2(Vertex_handle v0, 
                                           Vertex_handle v1,
                                           Vertex_handle v2): Fb(v0, v1, v2) {}

    Triangulation_face_base_with_FaceInfo2(Vertex_handle v0, 
                                           Vertex_handle v1,
                                           Vertex_handle v2, 
                                           Face_handle   n0, 
                                           Face_handle   n1,
                                           Face_handle   n2 ) : Fb(v0, v1, v2, n0, n1, n2) {}

    const Info& info() const { return _info; }
    Info&       info()       { return _info; }
  private:
    Info _info;
  };
#else
  struct FaceInfo2
  {
    int nesting_level;

    bool in_domain() 
    { 
      return nesting_level%2 == 1;
    }
  };
#endif  

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, K> Vb;
  #ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
  typedef Triangulation_face_base_with_FaceInfo2<K> Fbb;
  #else
  typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K> Fbb;
  #endif
  typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>   Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>      TDS;

  typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS> CDT2;
  #ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
  typedef T_FaceInfo2<typename CDT2::Face_handle> FaceInfo2;
  #endif
  typedef typename K::Point_2 Point_2;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Line_2 Line_2;


  typedef typename CDT2::Locate_type Locate_type;
  typedef typename CDT2::Vertex_handle Vertex_handle_2;
  typedef typename CDT2::Face_handle Face_handle_2;
  typedef typename CDT2::All_faces_iterator All_faces_iterator;
  typedef typename CDT2::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDT2::Finite_edges_iterator Finite_edges_iterator;
  typedef typename CDT2::Finite_faces_iterator Finite_faces_iterator;
  typedef typename CDT2::Edge Edge_2;
  typedef cpp0x::tuple<Vertex_handle_2,Vertex_handle_2,Vertex_handle_2> Vertex_handle_triple;

  typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_kernel;
  typedef CGAL::Cartesian_converter<K,Exact_kernel> To_exact;
  typedef CGAL::Cartesian_converter<Exact_kernel,K> To_input;

  // For triangle-vertex cells the type is TOP_TRIANGLE or BOTTOM_TRIANGLE
  // and the cell is defined by f0 and vh 
  // For edge-edge cells the cell is defined by Edge(f0,i0) and Edge(f1,i1)
  struct CellInfo3 {
    enum Type { TOP_TRIANGLE, BOTTOM_TRIANGLE, EDGE_EDGE };
    Type type;
    
    int index;
    Face_handle_2 f0, f1;
    int i0, i1;
    Vertex_handle_2 vh;
    bool cc;
    bool volume;
    #if !defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT) && defined(CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS)
    bool non_manifold_features;
    bool has_a_nm_edge_with_manifold_vertices;
    bool isolated_triangle;
    CellInfo3():non_manifold_features(false),has_a_nm_edge_with_manifold_vertices(false),isolated_triangle(false),m_nmfacet(false){}
  private:
    std::bitset<6> nm_edges;
    std::bitset<4> m_nmvertices;
    static const int nm_edge_index[4][4];//indices mapping in nm_edges
    bool m_nmfacet; //a simple bool as only one facet can be non-manifold (the one opposite to io)
  public:
    void set_nm_vertices(int i){
      non_manifold_features=true;
      m_nmvertices.set(i);
    }
  
    bool nm_vertices(int i) const{
      return m_nmvertices[i];
    }
    
    template <class Cell_handle_3>
    void set_nm_edge(Cell_handle_3 cell,int i1,int i2){
      non_manifold_features=true;
      nm_edges.set( CellInfo3::nm_edge_index[i1][i2] );
      if ( cell->vertex(i1)->info().v->info().on_contour && cell->vertex(i2)->info().v->info().on_contour )
        has_a_nm_edge_with_manifold_vertices=true;
    }
    
    void set_nm_edge_with_manifold_vertices(int i1,int i2){
      non_manifold_features=true;
      nm_edges.set( CellInfo3::nm_edge_index[i1][i2] );  
      has_a_nm_edge_with_manifold_vertices=true;
    }
    
    template <class Cell_handle_3>
    void set_nm_facet(Cell_handle_3 cell) {
      CGAL_precondition( type!=EDGE_EDGE );
      CGAL_precondition( f0->info().in_domain() );
      m_nmfacet=true;
      non_manifold_features=true;
      const int indices[3]={(i0+1)%4,(i0+2)%4,(i0+3)%4};
      const Vertex_handle_2 vertices_2[3]= { cell->vertex(indices[0])->info().v,
                                             cell->vertex(indices[1])->info().v,
                                             cell->vertex(indices[2])->info().v };
      const bool on_contour[3] = { vertices_2[0]->info().on_contour, vertices_2[1]->info().on_contour, vertices_2[2]->info().on_contour };
      

     //collect constrained edges
      std::set< std::pair<Vertex_handle_2,Vertex_handle_2> > constrained_edges;
      for (int i=0;i<3;++i)
        if ( f0->is_constrained(i) ) constrained_edges.insert( make_sorted_pair(f0->vertex((i+1)%3),f0->vertex((i+2)%3) ) );

      if (constrained_edges.size()==3){
        isolated_triangle=true;
        return;
      }
    
      //non-constrained edges are non-manifold
      //edge 01
      if ( on_contour[0] && on_contour[1] ){
        if ( constrained_edges.find( make_sorted_pair(vertices_2[0],vertices_2[1]) )==constrained_edges.end() )
          set_nm_edge_with_manifold_vertices(indices[0], indices[1]);
      }
      else
        set_nm_edge(cell,indices[0], indices[1]);
      //edge 02
      if ( on_contour[0] && on_contour[2] ){
        if ( constrained_edges.find( make_sorted_pair(vertices_2[0],vertices_2[2]) )==constrained_edges.end() )
          set_nm_edge_with_manifold_vertices(indices[0], indices[2]);
      }
      else
        set_nm_edge(cell,indices[0], indices[2]);
      //edge 21
      if ( on_contour[2] && on_contour[1] ){
        if ( constrained_edges.find( make_sorted_pair(vertices_2[2],vertices_2[1]) )==constrained_edges.end() )
          set_nm_edge_with_manifold_vertices(indices[2], indices[1]);
      }
      else
        set_nm_edge(cell,indices[2], indices[1]);

      //a vertex that is not on a contour is non-manifold
      if ( !on_contour[0] ) set_nm_vertices(indices[0]);
      if ( !on_contour[1] ) set_nm_vertices(indices[1]);
      if ( !on_contour[2] ) set_nm_vertices(indices[2]);
    }
    
    bool nm_facet() const {return m_nmfacet;}
    bool nm_edge(int i1,int i2) const { return nm_edges[ CellInfo3::nm_edge_index[i1][i2] ]; }
    #endif
    #warning replace volume, cc? and the array in facets_already_handled by a bitset! this will remove maps
    //used for orienting dual edges to eliminate edge_edge nodes
    //when a bit i is set, this indicates that there exists a solid path (crossing only facets)
    //going through Facet(cell,i) to reach the bottom or the top layer
    //(we follow the strategy from INRIA RR546 p19-21)
    std::bitset<4> out_face_status; 
  };

  struct VertexInfo3
  {
    Vertex_handle_2 v;
  };



  typedef CGAL::Triangulation_data_structure_3<CGAL::Triangulation_vertex_base_with_info_3<VertexInfo3,K>, CGAL::Triangulation_cell_base_with_info_3<CellInfo3, K> > Tds_3;
  typedef CGAL::Delaunay_triangulation_3<K, Tds_3> DT3;

  typedef typename DT3::Vertex_handle Vertex_handle_3;
  typedef typename DT3::Cell_handle Cell_handle_3;
  typedef typename DT3::Facet Facet_3;
  typedef typename DT3::Finite_cells_iterator Cell_iterator_3;
  typedef typename DT3::All_cells_iterator All_cells_iterator_3;
  typedef typename DT3::Finite_vertices_iterator Vertex_iterator_3;
  typedef typename DT3::Edge Edge_3;
  
  // Index the vertices for writing polyline files
  void
  index(const CDT2& cdt)
  {
    int i = 0;
    for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
      it->info().index = i++;
      it->info().poly_index = -1;
    }
  }


  void 
  mark_domains(const CDT2& ct, 
               Face_handle_2 start, 
               int index, 
               std::list<Edge_2>& border )
  {
    if(start->info().nesting_level != -1){
      return;
    }
    std::list<Face_handle_2> queue;
    queue.push_back(start);

    while(! queue.empty()){
      Face_handle_2 fh = queue.front();
      queue.pop_front();
      if(fh->info().nesting_level == -1){
        fh->info().nesting_level = index;
        for(int i = 0; i < 3; i++){
          Edge_2 e(fh,i);
          Face_handle_2 n = fh->neighbor(i);
          if(n->info().nesting_level == -1){
            if(ct.is_constrained(e)){
              border.push_back(e);
            } else {
              queue.push_back(n);
            }
          }
        }
      }
    }
  }


  void
  mark_domains(const CDT2& cdt)
  {
    for(All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
      it->info().nesting_level = -1;
    }

    int index = 0;
    std::list<Edge_2> border;
    mark_domains(cdt, cdt.infinite_face(), index++, border);
    while(! border.empty()){
      Edge_2 e = border.front();
      border.pop_front();
      Face_handle_2 n = e.first->neighbor(e.second);
      if(n->info().nesting_level == -1){
        mark_domains(cdt, n, e.first->info().nesting_level+1, border);
      }
    }
  }

  Point_3 centroid(const Point_3& p1,const Point_3& p2,const Point_3& p3) const{
    //TODO: check that the computed centroid is really inside fh
    To_exact to_exact;
    To_input to_input;
    return to_input(
      Exact_kernel::Construct_centroid_3() (to_exact(p1),to_exact(p2),to_exact(p3))
    );
  }
  
  Point_2 circumcenter(Face_handle_2 fh) const {
    To_exact to_exact;
    To_input to_input;
    return to_input(
      Exact_kernel::Construct_circumcenter_2() (
        to_exact(fh->vertex(0)->point()),
        to_exact(fh->vertex(1)->point()),
        to_exact(fh->vertex(2)->point())
      )
    );
  }

  Point_2 centroid(Face_handle_2 fh) const {
    //TODO: check that the computed centroid is really inside fh
    To_exact to_exact;
    To_input to_input;
    return to_input(
      Exact_kernel::Construct_centroid_2() (
        to_exact(fh->vertex(0)->point()),
        to_exact(fh->vertex(1)->point()),
        to_exact(fh->vertex(2)->point())
      )
    ); 
  }
  
  Point_2 midpoint(const Point_2& p1,const Point_2& p2) const {
    //TODO: check that the computed midpoint is really inside the edge (cmp distances)
    To_exact to_exact;
    To_input to_input;
    return to_input(
      Exact_kernel::Construct_midpoint_2() (
        to_exact(p1),
        to_exact(p2)
      )
    );
  }

  Point_3 midpoint(const Point_3& p1,const Point_3& p2) const {
    //TODO: check that the computed midpoint is really inside the edge (cmp distances)
    To_exact to_exact;
    To_input to_input;
    return to_input(
      Exact_kernel::Construct_midpoint_3() (
        to_exact(p1),
        to_exact(p2)
      )
    );
  }

  Point_2 midpoint(Vertex_handle_2 v1,Vertex_handle_2 v2) const {
    return midpoint(v1->point(),v2->point());
  }

  
  typedef std::map<std::pair<Vertex_handle_2,Vertex_handle_2>,std::list< std::pair<Point_2,Edge_2> > > Intersected_edges_map;
  
  //predicate used to sort points along a segment, when the points are not exactly
  //on the segment.
  struct Compare_distance_to_point{
    Point_2 ref;
    typedef bool result_type;
    Compare_distance_to_point(const Point_2& r):ref(r){}
    bool operator()(const Point_2& p1,const Point_2& p2) const{
      return CGAL::has_smaller_distance_to_point(ref,p1,p2);
    }
  };
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS  
  template <class RayOrSegment>
  Point_2 intersection_point(const RayOrSegment& s1,const Segment_2& s2) const {
    To_exact to_exact;
    CGAL::Object obj = CGAL::intersection(to_exact(s1),to_exact(s2));
    const Exact_kernel::Point_2* p = CGAL::object_cast<Exact_kernel::Point_2>(&obj);
    CGAL_assertion(p!=NULL);
    To_input to_input;
    return to_input(*p);
  }

  //turn around the vertex f_start->vertex(v_index) to find an edge intersected by object.
  //We know the segment/ray object passes through f_start->vertex(v_index).
  //This function returns true if f_start->vertex(v_index) lies on a constraint.
  //f_start and v_index are updated to indicate the edge intersected, as well as segment.
  template<class RayOrSegment>
  bool find_intersected_edge_in_vertex_link(const CDT2& cdt,Face_handle_2& f_start,int& v_index,Segment_2& segment,const RayOrSegment& object){
    int n_index = (v_index+1)%3; //turn around the vertex v_index
    Face_handle_2 n_f=f_start;
    int n_vi=v_index;
    
    if ( n_f->is_constrained(n_index) ) return true;
    
    do{
      boost::tie(n_f,n_index)=cdt.mirror_edge(Edge_2(n_f,n_index));
      n_vi=n_f->index( f_start->vertex(v_index) );
      n_index=(n_index+1)%3==n_vi?(n_index+2)%3:(n_index+1)%3;
      
      if ( n_f->is_constrained(n_index) ) return true;
      
      segment=Segment_2(n_f->vertex((n_vi+1)%3)->point(),n_f->vertex((n_vi+2)%3)->point());
      if ( CGAL::do_intersect(object,segment) ){
        f_start=n_f;
        v_index=n_vi;
        return false;
      }
      
    }while(n_f!=f_start);
    CGAL_assertion(n_f!=f_start);
    return false;
  }
  
  //start from a vertex inside a face and look for intersection points between a ray or a segment and triangulation edges.
  //we stop when we reach a constrained edge
  template <class RayOrSegment,class PointOutputIterator>
  void mark_edges_intersected(const CDT2& cdt,Face_handle_2 f_in,int opt_e_index,const RayOrSegment& object,bool target_on_constraint,Intersected_edges_map& edges_intersected,PointOutputIterator out){
    int k=0;
    Segment_2 segment;
    Point_2 ip;

    if ( opt_e_index ==-1 ){
      //the source of object is strictly inside a face
      for (;k<3;++k){
        segment=Segment_2(f_in->vertex((k+1)%3)->point(),f_in->vertex((k+2)%3)->point());
        if ( CGAL::do_intersect(object,segment) ) break;
      }
      CGAL_assertion(k!=3);
    }
    else{
      if(opt_e_index>2){
        //the source of object is on a vertex
        k=opt_e_index;
        do{
          k=k%3;
          if ( find_intersected_edge_in_vertex_link(cdt,f_in,k,segment,object) )
            return; //on a constrained edge
          ip=intersection_point(object,segment);
          
          if ( ip == f_in->vertex( (k+1)%3 )->point() ) k=(k+1)%3+3;
          else if ( ip == f_in->vertex( (k+2)%3 )->point() ) k=(k+2)%3+3;
          //here it should be notice that we cannot enter an infinite loop because in
          //find_intersected_edge_in_vertex_link we first change the current face, so
          //there is not risk to find an already encountered vertex on object.
        } while ( k>2 );
      }
      else{
        //the source of object is on an edge: we need to find which incident faces is the one
        //to consider to start the search.
        for (k=1;k<3;++k){
          segment=Segment_2(f_in->vertex((opt_e_index+k+1)%3)->point(),f_in->vertex((opt_e_index+k+2)%3)->point());
          if ( CGAL::do_intersect(object,segment) ) break;
        }
        if (k==3){
          k=0;
          cpp0x::tie(f_in,opt_e_index)=cdt.mirror_edge(Edge_2(f_in,opt_e_index));
          for (k=1;k<3;++k){
            segment=Segment_2(f_in->vertex((opt_e_index+k+1)%3)->point(),f_in->vertex((opt_e_index+k+2)%3)->point());
            if ( CGAL::do_intersect(object,segment) ) break;
          }
          CGAL_assertion(k!=3);
        }
        k=(opt_e_index+k)%3;
      }
    }

    do{
      CGAL_assertion(f_in->info().in_domain());
      ip=intersection_point(object,segment);
      while ( ip==segment.source() || ip==segment.target() ){
        //the segment/ray passes through a vertex, we turn around
        //while to find another edge or that the vertex lie on a constrained edge
        k = ip==f_in->vertex((k+1)%3)->point()?(k+1)%3:(k+2)%3;
        if ( find_intersected_edge_in_vertex_link(cdt,f_in,k,segment,object) )
          return; //on a constrained edge
        ip=intersection_point(object,segment);
      }
      if ( f_in->is_constrained(k) ) break;
      //do not insert point that are too close from a constraint or a point
      if ( min_distance_to_vertices(ip,f_in)   > m_min_point_squared_distance && 
           min_distance_to_constaints(ip,f_in) > m_min_point_squared_distance 
      ) *out++=ip;
      boost::tie(f_in,k)=cdt.mirror_edge(Edge_2(f_in,k));
      int i=1;
      for (;i!=3;++i){
        segment=Segment_2(f_in->vertex(k)->point(),
                          f_in->vertex((k+i)%3)->point());
        if ( CGAL::do_intersect(object,segment) ){
          k=(k+ (i==1?2:1) )%3;
          break;
        }
      }
      CGAL_assertion(i!=3);      
    }
    while(true);

    CGAL_assertion( f_in->is_constrained(k) );
    
    if (target_on_constraint) return;
    
    //insert the new point on the constraint only if not too close to an endpoint
    if ( min_distance_to_vertices(ip,f_in,k) > m_min_point_squared_distance ){
      typename Intersected_edges_map::iterator it_map=edges_intersected.insert(
        std::make_pair(make_sorted_pair(f_in->vertex((k+1)%3),f_in->vertex((k+2)%3)),
        std::list< std::pair<Point_2,Edge_2> >())
      ).first;
      it_map->second.push_back(std::make_pair(ip,Edge_2(f_in,k)));
    }
  }
#endif //DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS  
  
  //WARNING USE INTERVAL ARITHMETIC?
  double min_distance_to_vertices(const Point_2& p, Face_handle_2 fh) const {
    double d1=CGAL::squared_distance(p,fh->vertex(0)->point());
    double dmin=CGAL::squared_distance(p,fh->vertex(1)->point());
    dmin = (std::min) (d1,dmin);
    d1=CGAL::squared_distance(p,fh->vertex(2)->point());
    return (std::min) (d1,dmin);
  }

  double min_distance_to_vertices(const Point_2& p, Face_handle_2 fh,int i) const {
    return (std::min) (
      CGAL::squared_distance(p,fh->vertex((i+1)%3)->point()),
      CGAL::squared_distance(p,fh->vertex((i+2)%3)->point()) 
    );
  }
  
//  double min_distance_to_supporting_lines(const Point_2& p,Face_handle_2 fh) const {
//    double dmin=(std::min) (
//      CGAL::squared_distance(p,Line_2(fh->vertex(0)->point(),fh->vertex(1)->point())),
//      CGAL::squared_distance(p,Line_2(fh->vertex(1)->point(),fh->vertex(2)->point()))
//    );
//    return (std::min) ( dmin,
//      CGAL::squared_distance(p,Line_2(fh->vertex(0)->point(),fh->vertex(2)->point()))
//    );
//  }

  double distance_to_supporting_line(const Point_2& p,Face_handle_2 fh,int i) const {
    return CGAL::squared_distance(p,Line_2( fh->vertex((i+1)%3)->point(),
                                            fh->vertex((i+2)%3)->point()));
  }
  
  double min_distance_to_constaints(const Point_2& p, Face_handle_2 f) const {
    double min=std::numeric_limits<double>::max();
    for (int k=0;k<3;++k){
      if ( f->is_constrained(k) )
        min=(std::min) ( min,distance_to_supporting_line(p,f,k) );
    }
    return min;
  }

  #ifndef DO_NOT_FILTER_NOTCHES
  //used to collect points of the medial axis and their dual that should be inserted in cdtB
  typedef std::map<Face_handle_2,std::pair<Point_2,bool> > Graph_G;
  void vertex_to_add(const Point_2& center,FaceInfo2& A_info,Face_handle_2 A_face,CDT2& cdtB, Graph_G& graph_G,Intersected_edges_map& edges_intersected)
  #else
  #ifdef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
  void vertex_to_add(const Point_2& center,FaceInfo2&,Face_handle_2,CDT2& cdtB, std::list<Point_2>& points,Intersected_edges_map& edges_intersected)
  #else
  void vertex_to_add(const Point_2& center,FaceInfo2& A_info,Face_handle_2,CDT2& cdtB, std::list<Point_2>& points,Intersected_edges_map& edges_intersected)
  #endif  
  #endif  
  {
    Locate_type lt;
    int li;
    Face_handle_2 loc = cdtB.locate(center, lt, li);
    if( (lt == CDT2::FACE) && (loc->info().in_domain()) ){
      bool do_insert =  min_distance_to_vertices(center,loc)   > m_min_point_squared_distance && 
                        min_distance_to_constaints(center,loc) > m_min_point_squared_distance;
      #ifndef DO_NOT_FILTER_NOTCHES
      graph_G.insert(std::make_pair(A_face,std::make_pair(center,do_insert)));
      #else
      if ( do_insert ) points.push_back(center);
      #endif
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
      A_info.other_face = loc;
#endif
    }
    else{
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
      A_info.other_face=NULL;
#endif
      if (lt==CDT2::EDGE){
        if( loc->is_constrained(li) ){
          if ( min_distance_to_vertices(center,loc,li) > m_min_point_squared_distance )
          {
            typename Intersected_edges_map::iterator it_map=edges_intersected.insert(
                    std::make_pair(make_sorted_pair(loc->vertex((li+1)%3),loc->vertex((li+2)%3)),
                    std::list< std::pair<Point_2,Edge_2> >())
                  ).first;
            it_map->second.push_back(std::make_pair(center,Edge_2(loc,li)));          
          }
          #ifndef DO_NOT_FILTER_NOTCHES
          graph_G.insert(std::make_pair(A_face,std::make_pair(center,false)));
          #endif
          
          #ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
          A_info.dual_point_on_constraint=true;
          #endif
        }
        else{
          if (loc->info().in_domain() ){
            bool do_insert=min_distance_to_vertices(center,loc,li) > m_min_point_squared_distance;
            #ifndef DO_NOT_FILTER_NOTCHES
            graph_G.insert(std::make_pair(A_face,std::make_pair(center,do_insert)));
            #else
            if ( do_insert ) points.push_back(center);
            #endif            
            #ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
            A_info.other_face=loc;
            A_info.other_index=li;
            #endif
          }
        }
      }
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
      if (lt==CDT2::VERTEX){
        if ( loc->vertex(li)->info().on_contour ) A_info.dual_point_on_constraint=true;
        else{
          if (loc->info().in_domain()){
            A_info.other_face=loc;
            A_info.other_index=li+3;
          }
        }
      }
#endif
    }
  }
  
//look for outside voronoi vertices of cdtA that are inside cdtB
  void
  vertices_to_add(const CDT2& cdtA, CDT2& cdtB, std::list<Point_2>& points) 
  {
    //we use the following map to collect points that should be inserted on a constraint.
    Intersected_edges_map edges_intersected;
    
    #ifndef DO_NOT_FILTER_NOTCHES
    //used to collect points of the medial axis and their dual that should be inserted in cdtB
    typedef std::map<Face_handle_2,std::pair<Point_2,bool> > Graph_G;
    Graph_G graph_G;
    #endif
    
    for(Finite_faces_iterator it = cdtA.finite_faces_begin(); it != cdtA.finite_faces_end(); ++it){
      if(! it->info().in_domain()){
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
         it->info().dual_point_on_constraint=false;
#endif
        //exact construction of the circumcenter
        Point_2 center=circumcenter(it);
        #ifndef DO_NOT_FILTER_NOTCHES
        vertex_to_add(center,it->info(),it,cdtB,graph_G,edges_intersected);
        #else
        vertex_to_add(center,it->info(),it,cdtB,points,edges_intersected);
        #endif
//TODO: think about why the following is not improving the result, it correspond to a
//medial axis point! ---> I think the angle should be less than PI        
//do not use the following if DO_NOT_FILTER_NOTCHES is not defined
//        for (int k=0;k<3;++k){
//          if (it->is_constrained((k+1)%3) && it->is_constrained((k+2)%3)){
//            FaceInfo2 info;
//            vertex_to_add(it->vertex(k)->point(),info,NULL,cdtB,points);
//          }
//        }
          
        
      }
    }

    #ifndef DO_NOT_FILTER_NOTCHES
    //we first compute the bbox of cdtB
    CGAL::Bbox_2 bbox_cdtB=cdtB.finite_vertices_begin()->point().bbox();
    for (typename CDT2::Finite_vertices_iterator itpt=cpp0x::next(cdtB.finite_vertices_begin()),
                                                itpt_end=cdtB.finite_vertices_end();itpt!=itpt_end;++itpt
    ) bbox_cdtB=bbox_cdtB+itpt->point().bbox();

    //we recover connected component of the medial axis.
    //we filter only parts that are incident to the infinite facets and that are
    //small compare to contour in cdtB
    std::set<Face_handle_2> visited;
    for ( typename Graph_G::iterator ith=graph_G.begin(),ith_end=graph_G.end();ith!=ith_end;++ith ){
      bool is_connected_to_bottom_contour=false;
      if ( !visited.insert(ith->first).second ) continue; //already processed
      std::list< typename Graph_G::iterator > heap;
      heap.push_back(ith);
      std::list<typename Graph_G::iterator> cc_points;
      while(!heap.empty()){
        cc_points.splice(cc_points.end(),heap,heap.begin());
        if ( cc_points.back()->first->info().dual_point_on_constraint && 
             cc_points.back()->first->info().nesting_level==0   ) is_connected_to_bottom_contour=true;
        for (int i=0;i<3;++i){
          Face_handle_2 neighbor=cc_points.back()->first->neighbor(i);
          typename Graph_G::iterator itn=graph_G.find(neighbor);
          if ( itn==graph_G.end() ){ //the circumcenter is not inside bottom's contour
            if ( neighbor->info().nesting_level==0 ) is_connected_to_bottom_contour=true;
            continue;
          }
          if( !visited.insert(neighbor).second ) continue; //already processed
          heap.push_back(itn);
        }
      }
      
      CGAL_assertion(!cc_points.empty());
      CGAL::Bbox_2 bbox_ma=(*cc_points.begin())->second.first.bbox();
      for (typename std::list<typename Graph_G::iterator>::iterator itp=cc_points.begin(),itp_end=cc_points.end();itp!=itp_end;++itp)
        bbox_ma=bbox_ma+(*itp)->second.first.bbox();
     
      #warning Tune this parameter. Can't we use the diagonal length instead, or the area,...?
      if ( !is_connected_to_bottom_contour ||
        (std::max) (bbox_ma.xmax()-bbox_ma.xmin(),bbox_ma.ymax()-bbox_ma.ymin()) > 
        (std::min) (bbox_cdtB.xmax()-bbox_cdtB.xmin(),bbox_cdtB.ymax()-bbox_cdtB.ymin()) / m_bbox_ratio_ma_filtering
      )
        for (typename std::list<typename Graph_G::iterator>::iterator itp=cc_points.begin(),itp_end=cc_points.end();itp!=itp_end;++itp)
          if ( (*itp)->second.second ) 
            points.push_back((*itp)->second.first);
    }
    #endif

#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
    //We consider voronoi vertices of cdtA that when projected are inside the domain of cdtB.
    //In cdtA, we look at each of the facet's neighbors and consider those having its dual voronoi
    //vertex projected outside the domain of cdtB. Then we look for all intersection in cdtB with
    //the dual segment or ray (and stop when a constrained edge have been found)
    for(Finite_faces_iterator it = cdtA.finite_faces_begin(); it != cdtA.finite_faces_end(); ++it){
      if( !it->info().in_domain() && 
          it->info().other_face != NULL )
      {
        CGAL_assertion(!it->info().dual_point_on_constraint);
        for (int i=0;i<3;++i){
          if ( !cdtA.is_infinite(it->neighbor(i)) )
          {
            //handle segments here
            if ( !it->neighbor(i)->info().in_domain()  && //consider only outer medial axis
                 #ifndef DO_NOT_FILTER_NOTCHES
                 //TAG-IFMA I am no longer sure that points inserted inside the contour are helpful
                 !it->neighbor(i)->info().dual_point_on_constraint &&  //this point is not already on a constraint: this is only relevant to filter these elements
                                                                       //when we are not interested in intersection points of non-constrained edges with the exterior medial axis
                 #endif
                 it->neighbor(i)->info().other_face == NULL) //dual vertex is outside the domain of cdtB
            {
              Face_handle_2 f_in = it->info().other_face;
              CGAL_assertion(!cdtB.is_infinite(f_in));
              CGAL_assertion(f_in->info().in_domain());
              Point_2 c1=circumcenter(it);
              Point_2 c2=circumcenter(it->neighbor(i));
              CGAL_assertion(!it->info().in_domain());
              CGAL_assertion(!it->neighbor(i)->info().in_domain());
              Segment_2 seg(c1,c2);
              #ifdef DO_NOT_FILTER_NOTCHES
              mark_edges_intersected(cdtB,f_in,it->info().other_index,seg,it->neighbor(i)->info().dual_point_on_constraint,edges_intersected,std::back_inserter(points));
              #else
              //TAG-IFMA I am no longer sure that points inserted inside the contour are helpful, if yes then I need to
              //adapt the mechanism with graph_G
              mark_edges_intersected(cdtB,f_in,it->info().other_index,seg,false,edges_intersected,Emptyset_iterator());
              #endif
            }
          }
          else{
            //handle rays here
            Face_handle_2 f_in = it->info().other_face;
            CGAL_assertion(!cdtB.is_infinite(f_in));
            CGAL_assertion(f_in->info().in_domain());
            Point_2 c1=circumcenter(it);
            To_exact to_exact;
            To_input to_input;
            
            Exact_kernel::Vector_2 ortho_vector=
              (to_exact(it->vertex((i+2)%3)->point())-to_exact(it->vertex((i+1)%3)->point())).perpendicular(CGAL::CLOCKWISE);
            K::Ray_2 ray(c1,to_input(ortho_vector));
            #ifdef DO_NOT_FILTER_NOTCHES
            mark_edges_intersected(cdtB,f_in,it->info().other_index,ray,false,edges_intersected,std::back_inserter(points));
            #else
            //TAG-IFMA I am no longer sure that points inserted inside the contour are helpful, if yes then I need to
            //adapt the mechanism with graph_G              
            mark_edges_intersected(cdtB,f_in,it->info().other_index,ray,false,edges_intersected,Emptyset_iterator());
            #endif
          }
        }
      }
    }
#endif //DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
    
    //In the following we insert points that are on constraints. We first use stored face handles to unconstrain 
    //all edges. If several points are on the same edge, we sort them using the distance to one extremity of the edge,
    //insert sequentially the points and constrain the edges at the same time (wo do not use face handles stored that are
    //potentially invalid.
//    std::set<Face_handle_2> removed_facets;
    for (typename Intersected_edges_map::iterator it=edges_intersected.begin(),end=edges_intersected.end();it!=end;++it)
    {
//      Edge_2 e=it->second.front().second;
//      if (removed_facets.find(e.first)!=removed_facets.end()){
        Edge_2 e;
        CGAL_assertion_code( bool b= )
          cdtB.is_edge(it->first.first,it->first.second,e.first,e.second);
        CGAL_assertion(b);
//      }
      CGAL_assertion(e.first->is_constrained(e.second));
//      cdtB.remove_constrained_edge(e.first,e.second,std::inserter(removed_facets,removed_facets.begin()));
      cdtB.remove_constrained_edge(e.first,e.second);
    }
    
    for (typename Intersected_edges_map::iterator it=edges_intersected.begin(),end=edges_intersected.end();it!=end;++it)
    {
      if (CGAL::cpp0x::next(it->second.begin())==it->second.end()){
        Vertex_handle_2 vh=cdtB.insert(it->second.front().first);
        vh->info().z = it->first.first->info().z;
        cdtB.insert_constraint(it->first.first,vh);
        cdtB.insert_constraint(it->first.second,vh);
      }
      else{
        std::list< std::pair<Point_2,Edge_2> >& lst=it->second;
        std::vector<Point_2> on_edge;
        on_edge.reserve(lst.size());
        for (typename std::list< std::pair<Point_2,Edge_2> >::iterator itl=lst.begin();itl!=lst.end();++itl)
          on_edge.push_back(itl->first);
        std::sort(on_edge.begin(),on_edge.end(),Compare_distance_to_point(it->first.first->point()));
        Vertex_handle_2 prev=it->first.first;
        CGAL_assertion( CGAL::squared_distance(prev->point(),on_edge[0]) < CGAL::squared_distance(prev->point(),on_edge[1]) );
        #warning TODO do not insert point if too close from the previous one!
        for (std::size_t i=0;i<on_edge.size();++i)
        {
          Vertex_handle_2 vh=cdtB.insert(on_edge[i]);
          vh->info().z = it->first.first->info().z;
          cdtB.insert_constraint(prev,vh);
          prev=vh;
        }
        cdtB.insert_constraint(prev,it->first.second);
      }
    }
  }

  void
  add_vertices_inside(CDT2& cdt, const std::list<Point_2>& points)
  {
    double z = cdt.finite_vertices_begin()->info().z;
    for(std::list<Point_2>::const_iterator it = points.begin(); it != points.end(); ++it){
      Locate_type loc;
      int li;
      Face_handle_2 fh = cdt.locate(*it,loc,li);
      bool on_constrained_edge = (loc==CDT2::EDGE) &&
                                  fh->is_constrained(li);
      #warning:  TODO
      //TODO: This should be temporary. We should consider Graph_G and do the following for each connected component of the medial axis,
      //      build a graph and call a simplification algorithm (look at boost graph for example)
      if ( min_distance_to_vertices(*it,fh)   < m_min_point_squared_distance ) continue;
      
      Vertex_handle_2 vh = cdt.insert(*it,loc,fh,li);
      vh->info().z = z;
      vh->info().on_contour=on_constrained_edge;
    }
  }

  // After refining a CDT we have to set the z values for the new vertices 
  void
  update_z(CDT2& cdt)
  {
    double z = cdt.finite_vertices_begin()->info().z;
    for(typename CDT2::All_vertices_iterator it = cdt.all_vertices_begin();  it != cdt.all_vertices_end(); ++it){
      it->info().z = z;
    }
  }

  void 
  classify(Cell_handle_3 ch, const CDT2& top, const CDT2& bottom)
  {
    Vertex_handle_2 b[3], t[3];
    int ti=0, bi=0, tii=0, bii=0;
    double topz = top.finite_vertices_begin()->info().z;
    int nr_top = 0;
    for(int i=0; i<4; i++){
      if(ch->vertex(i)->point().z() == topz){
        t[ti++] = ch->vertex(i)->info().v;
        ++nr_top;
        tii = i;
      } else {
        b[bi++] = ch->vertex(i)->info().v;
        bii = i;
      }
    }
    if(nr_top == 3){
      ch->info().type = CellInfo3::TOP_TRIANGLE;
      Face_handle_2 fh;
      CGAL_assertion_code(bool found = )
        top.is_face(t[0],t[1],t[2], fh);
      CGAL_assertion(found);
      ch->info().f0 = fh;
      ch->info().i0 = bii; // i0 is the index of the vertex which is alone
    } else if(nr_top == 1){
      ch->info().type = CellInfo3::BOTTOM_TRIANGLE;
      Face_handle_2 fh;
      CGAL_assertion_code(bool found = )
        bottom.is_face(b[0],b[1],b[2], fh);
      CGAL_assertion(found);
      ch->info().f0 = fh;
      ch->info().i0 = tii;
    } else {
      CGAL_assertion(nr_top == 2);
      ch->info().type = CellInfo3::EDGE_EDGE;
      Face_handle_2 fh;
      int fi=0;
    //handle top
      bool found = top.is_edge(t[0],t[1], fh, fi);
      CGAL_assertion(found);
      ch->info().f0 = fh;
      ch->info().i0 = fi;
      //if the two incident top facets are in the domain, the pencil of tetrahedra
      //forms a solid component to reach the top layer
      if (fh->info().in_domain() && fh->neighbor(fi)->info().in_domain() ){
        ch->info().out_face_status.set(bii);
        do{
          bii=(bii+1)%4;
        }while(ch->vertex(bii)->point().z()==topz);
        ch->info().out_face_status.set(bii);
      }
    //handle bottom
      found = bottom.is_edge(b[0],b[1], fh, fi);
      CGAL_assertion(found);
      ch->info().f1 = fh;
      ch->info().i1 = fi;
      //if the two incident top facets are in the domain, the pencil of tetrahedra
      //forms a solid component to reach the top layer
      if (fh->info().in_domain() && fh->neighbor(fi)->info().in_domain() ){
        ch->info().out_face_status.set(tii);
        do{
          tii=(tii+1)%4;
        }while(ch->vertex(tii)->point().z()!=topz);
        ch->info().out_face_status.set(tii);
      }
    }
  }


  void
  create_tetrahedrization(const CDT2& top, const CDT2& bottom)
  {
    std::vector<Vertex_handle_2> vertices;

    CGAL_precondition(top.is_valid());
    CGAL_precondition(bottom.is_valid());
    
    vertices.reserve(top.number_of_vertices()+bottom.number_of_vertices());
    for(Finite_vertices_iterator it = top.vertices_begin();
        it != top.vertices_end();
        ++it){
      vertices.push_back(it);
    }
    for(Finite_vertices_iterator it = bottom.vertices_begin();
        it != bottom.vertices_end();
        ++it){
      vertices.push_back(it);
    }
    CGAL_assertion(vertices.size() == top.number_of_vertices()+bottom.number_of_vertices());
    random_shuffle(vertices.begin(), vertices.end());
    //TODO use the thing with property maps.......
    for(typename std::vector<Vertex_handle_2>::iterator it = vertices.begin();
        it!= vertices.end();
        ++it){
      Vertex_handle_3 vh = delaunay_3.insert(Point_3((*it)->point().x(), (*it)->point().y(), (*it)->info().z));
      vh->info().v = *it;
    }

    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin();
        it != delaunay_3.finite_cells_end();
        ++it){
      classify(it, top, bottom);
    }
  }


  void remove_cell(Cell_handle_3 it)
  {
    it->info().volume = false;
  }

  //the tetrahedron fstart.first is a tetrahedron which does not belong to the reconstructed
  //volume. The function turns around edge, crossing fstart to orient dual edges to indicate
  //that a solid path can be reach in the opposite direction (start being removed).
  //The turn-around loop is stopped either when a deleted tetrahedron is encountered
  //(remaining edges are or will be oriented the same way), or a non-deleted T1 or T2 tetrahedron is found.
  template <class CellOutputIterator>
  void orient_dual_edges(const Facet_3& fstart,const Edge_3& edge,CellOutputIterator out)
  {
    Facet_3 fcurr=fstart;
    Facet_3 mirror=delaunay_3.mirror_facet(fcurr);
    do{
      int opposite_index=mirror.second;
      int new_index=opposite_index;
      do{
        new_index=(new_index+1)%4;
      }
      while(mirror.first->vertex(new_index)==edge.first->vertex(edge.second) || 
            mirror.first->vertex(new_index)==edge.first->vertex(edge.third) );
      fcurr=Facet_3(mirror.first,new_index);
      //orient the dual edge (or delete it for the first loop)
      fcurr.first->info().out_face_status.reset(opposite_index);

      //stop if a deleted node is encountered
      if (!fcurr.first->info().volume){
        mirror=delaunay_3.mirror_facet(mirror);
        //indicate that the former facet cannot be crossed to find a solid path
        mirror.first->info().out_face_status.reset(mirror.second);//this one is not needed (since the cell is deleted)
        break;
      }
      
      //stop at a non deleted T1 or T2 node 
      if (fcurr.first->info().type!=CellInfo3::EDGE_EDGE) break;
      
      mirror=delaunay_3.mirror_facet(fcurr);
      //if the edge is already oriented the opposite way
      if( mirror.first->info().out_face_status[mirror.second] ){
        *out++=fcurr.first;
        remove_cell(fcurr.first);
        break; //the rest of the job will be done when handling fcurr.first
      }
      
      //set edge orientation
      fcurr.first->info().out_face_status.set(fcurr.second);
    }while(true);  
  }
  
  //add a function that can be used to filter tetrahedra with high slope.
  //it is not used for the moment.
  bool is_T22_too_slim(Cell_handle_3 cell,double threshold=0.009) const {
    CGAL_precondition(cell->info().type==CellInfo3::EDGE_EDGE); 
    //get the edge in P1 and the edge in P2
    int ea0=0,ea1,eb0,eb1;
    if ( cell->vertex(ea0)->point().z() == cell->vertex(1)->point().z() ){
      ea1=1;eb0=2;eb1=3;
    }
    else{
      if ( cell->vertex(ea0)->point().z() == cell->vertex(2)->point().z() ){
        ea1=2;eb0=1;eb1=3;
      }
      else{
        ea1=3;eb0=1;eb1=2;
      }
    }
    
    double h2=CGAL::square(cell->vertex(ea0)->point().z() - cell->vertex(eb0)->point().z());
    double d2_0= (CGAL::max)(
      CGAL::squared_distance( cell->vertex(ea0)->point(),cell->vertex(eb0)->point() ),
      CGAL::squared_distance( cell->vertex(ea0)->point(),cell->vertex(eb1)->point() )
    );
    double d2_1= (CGAL::max)(
      CGAL::squared_distance( cell->vertex(ea1)->point(),cell->vertex(eb0)->point() ),
      CGAL::squared_distance( cell->vertex(ea1)->point(),cell->vertex(eb1)->point() )
    );
    return h2/(CGAL::max)(d2_0,d2_1) < threshold;
  }

  //function to remove T31 cell that are too flat and connected parts of the surface
  //that are far apart. It is not used for the moment.
  bool is_T31_too_slim(Cell_handle_3 cell,double threshold=0.07) const {
    CGAL_precondition(cell->info().type!=CellInfo3::EDGE_EDGE); 
    int i0=cell->info().i0;
    double h2=abs(cell->vertex(i0)->point().z() - cell->vertex((i0+1)%4)->point().z());
    #if 0
    double d2= (CGAL::min)(
      CGAL::squared_distance( cell->vertex(i0)->point(),cell->vertex((i0+1)%4)->point() ),
      CGAL::squared_distance( cell->vertex(i0)->point(),cell->vertex((i0+2)%4)->point() )
    );
    d2= (CGAL::min)( d2,CGAL::squared_distance( cell->vertex(i0)->point(),cell->vertex((i0+3)%4)->point() ) );
    return h2/d2 < threshold;
    #else
    Point_2 center_2=centroid(cell->info().f0);
    Point_3 center_3(center_2.x(),center_2.y(),cell->vertex((i0+1)%4)->point().z());
    if (h2/CGAL::squared_distance(center_3,cell->vertex(i0)->point()) < threshold){
      std::cout << "killing one" << std::endl;
    }
    return  h2/sqrt(squared_distance(center_3,cell->vertex(i0)->point())) < threshold;
    #endif
  }
  
  void remove_cells()
  {
    std::list<Cell_handle_3> infinite_nodes;
    for(All_cells_iterator_3 it = delaunay_3.all_cells_begin(); it != delaunay_3.all_cells_end(); ++it){
      bool is_infinite=delaunay_3.is_infinite(it);
      it->info().volume = ! is_infinite;
      it->info().cc = false;
      if ( is_infinite ){
        infinite_nodes.push_back(it);
      }
    }

    std::list<Cell_handle_3> eliminated_nodes;
    
    // First pass: Remove cells with an edge or face outside the domain 
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it ){
      if(it->info().type == CellInfo3::EDGE_EDGE){
        if( ( (! it->info().f0->info().in_domain()) && (! it->info().f0->neighbor(it->info().i0)->info().in_domain()) ) || 
            ( (! it->info().f1->info().in_domain()) && (! it->info().f1->neighbor(it->info().i1)->info().in_domain()) ) ){
          remove_cell(it);
          eliminated_nodes.push_back(it);
        }
      } else {
        if(! it->info().f0->info().in_domain()
        //     || is_T31_too_slim(it,0.1) 
        )
        {
          remove_cell(it);
          eliminated_nodes.push_back(it);
        }
      }
    }

    //Second pass - step 1: orient dual edges from infinite tetrahedra that are adjacent to 
    //a edge-edge tetrahedron.
    for (typename std::list<Cell_handle_3>::iterator cit=infinite_nodes.begin();cit!=infinite_nodes.end();++cit)
    {
      Cell_handle_3 ch=*cit;
      int inf_index=ch->index(delaunay_3.infinite_vertex());
      int indices[3]={(inf_index+1)%4,(inf_index+2)%4,(inf_index+3)%4};
      double z=ch->vertex( indices[0] )->point().z();
      Edge_3 edge;
     
      if (z == ch->vertex( indices[1] )->point().z() ){
        if (z == ch->vertex( indices[2] )->point().z() ) continue; //non adjacent to a edge_edge tetra
        edge=Edge_3(ch,indices[0],indices[1]);
      }
      else{
        if (z == ch->vertex( indices[2] )->point().z() ){
          edge=Edge_3(ch,indices[0],indices[2]);
        }
        else{
          CGAL_assertion(ch->vertex( indices[1] )->point().z() == ch->vertex( indices[2] )->point().z());
          edge=Edge_3(ch,indices[1],indices[2]);
        }
      }
      //indicate that the solid path should go the other way
      orient_dual_edges(Facet_3(ch,inf_index),edge,std::back_inserter(eliminated_nodes));
    }
     
    //Second pass : we repeat step 2 and final step while we remove a node in final step
    //OPTI: maintain a list a non-remove edge_edge tetrahedra
    while(!eliminated_nodes.empty()){
    //Second pass - step 2: From the list of eliminated_nodes, orient dual edges to indicate
    //which way a facet should be crossed to follow a solid path to the top or bottom layer.      
      while(!eliminated_nodes.empty()){
        Cell_handle_3 node=eliminated_nodes.front();
        eliminated_nodes.pop_front();
        
        
        if (node->info().type==CellInfo3::EDGE_EDGE)
        {
          //get the edge in P1 and the edge in P2
          int ea0=0,ea1,eb0,eb1;
          if ( node->vertex(0)->point().z() == node->vertex(1)->point().z() ){
            ea1=1;eb0=2;eb1=3;
          }
          else{
            if ( node->vertex(0)->point().z() == node->vertex(2)->point().z() ){
              ea1=2;eb0=1;eb1=3;
            }
            else{
              ea1=3; eb0=1;eb1=2;
            }
          }
          
          CGAL_assertion(node->vertex(ea0)->point().z()==node->vertex(ea1)->point().z() && 
                         node->vertex(eb0)->point().z()==node->vertex(eb1)->point().z() );
          
          //delete dual edges (one way the otherway will be done in the loop)
          node->info().out_face_status.reset();
              
          //turn around first edge to set dual edges orientation
          orient_dual_edges(Facet_3(node,eb0),Edge_3(node,ea0,ea1),std::back_inserter(eliminated_nodes));
          orient_dual_edges(Facet_3(node,eb1),Edge_3(node,ea0,ea1),std::back_inserter(eliminated_nodes));
          //turn around second edge to set dual edges orientation
          orient_dual_edges(Facet_3(node,ea0),Edge_3(node,eb0,eb1),std::back_inserter(eliminated_nodes));
          orient_dual_edges(Facet_3(node,ea1),Edge_3(node,eb0,eb1),std::back_inserter(eliminated_nodes));
        }
        else{
          //this orientation of facet should be done only for edges on a contour, but the current implementation
          //should not be that expensive since it is stopped as soon as a deleted tetrahedron is found
          int indices[3]={(node->info().i0+1)%4,(node->info().i0+2)%4,(node->info().i0+3)%4};
          orient_dual_edges(Facet_3(node,indices[0]),Edge_3(node,indices[1],indices[2]),std::back_inserter(eliminated_nodes));
          orient_dual_edges(Facet_3(node,indices[2]),Edge_3(node,indices[1],indices[0]),std::back_inserter(eliminated_nodes));
          orient_dual_edges(Facet_3(node,indices[1]),Edge_3(node,indices[0],indices[2]),std::back_inserter(eliminated_nodes));
        }
      }
      
      // Second pass - final step: Remove edge_edge cells that does not have two facets to cross to
      //reach the top and the bottom layers along a solid path
      for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){
        if(it->info().type == CellInfo3::EDGE_EDGE && it->info().volume){
          
          //get the edge in P1 and the edge in P2
          int ea0=0,ea1,eb0,eb1;
          if ( it->vertex(0)->point().z() == it->vertex(1)->point().z() ){
            ea1=1;eb0=2;eb1=3;
          }
          else{
            if ( it->vertex(0)->point().z() == it->vertex(2)->point().z() ){
              ea1=2;eb0=1;eb1=3;
            }
            else{
              ea1=3; eb0=1;eb1=2;
            }
          }
          

          if ( ( it->info().out_face_status.test(ea0) || it->info().out_face_status.test(ea1) )
             &&
               ( it->info().out_face_status.test(eb0) || it->info().out_face_status.test(eb1) )
          ) 
          {
            continue;
          }
          remove_cell(it);
          eliminated_nodes.push_back(it);
        }
      }
    }
    #ifndef CGAL_NO_EDGE_EDGE_EXTRA_REMOVAL
    //remove edge_edge tetrahedra remaining that are connected to at most two other tetrahedra.
    //I notice that this kind of tetrahedron appears to be in curved parts and make the surface
    //less smooth.
    bool removed=false;
    do{
      int buried_edge[4],exposed_edge[4];
      removed=false;
      for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it ){
        int* b_ptr=buried_edge;
        int* e_ptr=exposed_edge;        
        if(it->info().type == CellInfo3::EDGE_EDGE && it->info().volume){
          int nbn=0;
          for (int k=0;k<4;++k)
            if (it->neighbor(k)->info().volume ){
              ++nbn;
              *b_ptr++=k;
            }
            else
              *e_ptr++=k;
          if (nbn<=2){
            CGAL_assertion(nbn==2);
            //buried_edge and exposed_edge contains indices of edges in the T22 cell which
            //are exposed and buried respectively (and not a contour edge).
            //The dihedral angle computed is the one around these two edges.            
            double buried_angle= Mesh_3::dihedral_angle(
                it->vertex(DT3::next_around_edge(buried_edge[0],buried_edge[1]))->point(),
                it->vertex(DT3::next_around_edge(buried_edge[1],buried_edge[0]))->point(),
                it->vertex(buried_edge[1])->point(),
                it->vertex(buried_edge[0])->point()
            );
            
            double exposed_angle= Mesh_3::dihedral_angle(
                it->vertex(DT3::next_around_edge(exposed_edge[0],exposed_edge[1]))->point(),
                it->vertex(DT3::next_around_edge(exposed_edge[1],exposed_edge[0]))->point(),
                it->vertex(exposed_edge[1])->point(),
                it->vertex(exposed_edge[0])->point()
            );
            //remove the cell only if we improve the smoothness of the surface
            if (buried_angle > min_smooth_dihedral_angle && exposed_angle < max_sharp_dihedral_angle){
              remove_cell(it);
              removed=true;
            }
          }
        }
      }
    }
    while(removed);
    #endif

    // Third pass: Remove connected components which do not have at least one edge_edge cell
    std::list<Cell_handle_3> cells_to_remove;
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it ) {
      if(it->info().type != CellInfo3::EDGE_EDGE){
        if(! it->info().cc){
          // explore the component
          bool found_edge_edge_cell = false;
          std::list<Cell_handle_3> component;
          std::list<Cell_handle_3> queue;
          it->info().cc = true;
          component.push_back(it);
          queue.push_back(it);
          while(! queue.empty()){
            Cell_handle_3 ch = queue.front();
            queue.pop_front();
            for(int i=0; i < 4; i++){
              Cell_handle_3 nh = ch->neighbor(i);
              if(nh->info().volume){
                if(! nh->info().cc){
                  if(nh->info().type == CellInfo3::EDGE_EDGE){
                    found_edge_edge_cell = true;
                  }
                  nh->info().cc = true;
                  component.push_back(nh);
                  queue.push_back(nh);
                }
              }
            }
          }
          if(! found_edge_edge_cell)
            cells_to_remove.splice(cells_to_remove.end(),component);
        }
      }
    }
    for(typename std::list<Cell_handle_3>::iterator it = cells_to_remove.begin(); it != cells_to_remove.end(); ++it){
      remove_cell(*it);
    }

  }

  struct Compare_sorted_pair_2D_vertices_from_3D{ //OUTPUT
    typedef bool result_type;
    typedef std::pair<Vertex_handle_3,Vertex_handle_3> Input_type;
    typedef std::pair<Vertex_handle_2,Vertex_handle_2> Sorting_type;
    
    Sorting_type get_key(const Input_type& input) const {
      Vertex_handle_2 v1=input.first->info().v;
      Vertex_handle_2 v2=input.second->info().v;
      if (v1 < v2) return std::make_pair(v1,v2);
      return std::make_pair(v2,v1);
    }
    
    bool operator()(const Input_type& p1,const Input_type& p2) const {
      return get_key(p1) < get_key(p2);
    }
  };
  
  template <class Iterator>
  struct Cmp_deref{
    typedef bool result_type;
    result_type operator()(Iterator it1,Iterator it2) const{
      return &(*it1)<&(*it2);
    }
  };
  
  //pair must be sorted before find and insert
  typedef std::set<std::pair<Vertex_handle_2,Vertex_handle_2> > Unique_incontour_edges;
  
  #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  struct DS_for_incontour_edge_handling{
    struct Indices{
      int opp_e; //index of the vertex opposite to the edge in the facet
      int mp_index; //midpoint index in the vector
      Indices(int oe,int v_):opp_e(oe),mp_index(v_){}
    };

    
    typedef std::map<std::pair<Cell_handle_3,int>,Indices> Non_planar_facets_to_split;
    //the map that associates an edge to the midpoint index in the vector
    typedef std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int> Edge_to_midpoint_index;
    typedef std::list< std::pair<Edge_to_midpoint_index,Cell_handle_3> > Planar_facets_to_split;
    typedef std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int,Compare_sorted_pair_2D_vertices_from_3D> Edge_to_vertex_index_map;

    void clear(){
      non_planar_facets.clear();
      planar_facets.clear();
      last_vertex_index=0;
      edge_to_vertex_index.clear();
      midpoint_indices.clear();
      triangular_isolated_contours.clear();
    }
    
    Non_planar_facets_to_split non_planar_facets;
    Planar_facets_to_split planar_facets;
    std::list<Cell_handle_3> triangular_isolated_contours;
    
    Edge_to_vertex_index_map edge_to_vertex_index;
    std::vector<int> midpoint_indices;
    int last_vertex_index;
    
    DS_for_incontour_edge_handling():last_vertex_index(0){}
    void set_non_planar_facet(const std::pair<Vertex_handle_3,Vertex_handle_3>& edge,Cell_handle_3 cell,int findex,int opp_e){
      std::pair<typename Edge_to_vertex_index_map::iterator,bool> res=edge_to_vertex_index.insert(std::make_pair(edge,last_vertex_index));
      if (res.second) ++last_vertex_index;
      non_planar_facets.insert(std::make_pair(std::make_pair(cell,findex),Indices(opp_e,res.first->second)));
    }
    
    
    void set_planar_facet(const std::list<std::pair<Vertex_handle_3,Vertex_handle_3> >& edges_to_split,Cell_handle_3 top_t31){
      planar_facets.push_back( std::make_pair(Edge_to_midpoint_index(),top_t31) );
      Edge_to_midpoint_index& edges_map=planar_facets.back().first;
      for (typename std::list<std::pair<Vertex_handle_3,Vertex_handle_3> >::const_iterator 
        it=edges_to_split.begin(),it_end=edges_to_split.end();it!=it_end;++it)
      {
        CGAL_assertion(it->first < it->second);
        std::pair<typename Edge_to_vertex_index_map::iterator,bool> res=edge_to_vertex_index.insert(std::make_pair(*it,last_vertex_index));
        if (res.second) ++last_vertex_index;
        edges_map.insert(std::make_pair(*it,res.first->second)); 
      }
    }
    
    void prepare(){
      CGAL_assertion(midpoint_indices.empty());
      midpoint_indices=std::vector<int>(last_vertex_index,-1);
    }
  } incontour_edge_ds;

  #if !defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT)
  template <class Facet_to_nm_indices>
  void write_planar_facets_with_non_manifold_incontour_edges(std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>& edges_to_split,Cell_handle_3 top_t31,Facet_to_nm_indices& planar_facets_nm_point)
  #else
  void write_planar_facets_with_non_manifold_incontour_edges(std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>& edges_to_split,Cell_handle_3 top_t31)
  #endif
  {
    int i=top_t31->info().i0;
    int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
    Vertex_handle_3 vh[3]={top_t31->vertex(indices[0]),top_t31->vertex(indices[1]),top_t31->vertex(indices[2])};
    #if !defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT)
    typename Facet_to_nm_indices::iterator fah_it=planar_facets_nm_point.find(Facet_3(top_t31,top_t31->info().i0));
    
    if ( fah_it!=planar_facets_nm_point.end() ){
      for (int a=0;a<3;++a){
        if ( fah_it->second[a]!=-1 ) indices[a]=fah_it->second[a];
        else indices[a]=get_vertex_index(vh[a]->info().v);
      }
      if (i%2==0){
        std::swap(indices[0],indices[1]);
        std::swap(vh[0],vh[1]);
      }
      planar_facets_nm_point.erase(fah_it); //remove it to not print it twice
    }
    else
    #endif
    {
      for (int a=0;a<3;++a) indices[a]=get_vertex_index(vh[a]->info().v);
      if (i%2==0){
        std::swap(indices[0],indices[1]);
        std::swap(vh[0],vh[1]);
      }
    }
    
    switch( edges_to_split.size() )
    {
      case 3:
      {
        int mi[3]={ get_vertex_index_on_edge(vh[0],vh[1],edges_to_split[make_sorted_pair(vh[0],vh[1])]),
                    get_vertex_index_on_edge(vh[1],vh[2],edges_to_split[make_sorted_pair(vh[1],vh[2])]),
                    get_vertex_index_on_edge(vh[2],vh[0],edges_to_split[make_sorted_pair(vh[2],vh[0])])
                  };
        for (int k=0;k<3;++k)
          slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[k],mi[k],mi[(k+2)%3]));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(mi[0],mi[1],mi[2]));
        break;
      }
      case 2:
      {
        int mi[2];
        int k;
        std::pair<Vertex_handle_3,Vertex_handle_3> vp=make_sorted_pair(vh[0],vh[1]);
        typename std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>::iterator itmi=
          edges_to_split.find(vp);
        if (itmi!=edges_to_split.end()){
          mi[0]=get_vertex_index_on_edge(vp.first,vp.second,itmi->second);
          vp=make_sorted_pair(vh[1],vh[2]);
          itmi=edges_to_split.find(vp);
          if (itmi!=edges_to_split.end()){
            mi[1]=get_vertex_index_on_edge(vp.first,vp.second,itmi->second);
            std::swap(mi[0],mi[1]);
            k=1; 
          }
          else{
            vp=make_sorted_pair(vh[2],vh[0]);
            mi[1]=get_vertex_index_on_edge(vp.first,vp.second,edges_to_split[vp]);
            k=0;
          }
        }
        else{
          vp=make_sorted_pair(vh[2],vh[0]);
          mi[0]=get_vertex_index_on_edge(vp.first,vp.second,edges_to_split[vp]);
          vp=make_sorted_pair(vh[1],vh[2]);
          mi[1]=get_vertex_index_on_edge(vp.first,vp.second,edges_to_split[vp]);
          k=2;
        }
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[k],mi[0],mi[1]));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(mi[0],indices[(k+1)%3],indices[(k+2)%3]));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(mi[0],indices[(k+2)%3],mi[1]));
        break;
      }
      default:
      {
        CGAL_assertion(edges_to_split.size()==1);
        int k=0,mi=0;
        for (;k<3;++k){
          std::pair<Vertex_handle_3,Vertex_handle_3> vp=make_sorted_pair(vh[k],vh[(k+1)%3]);
          typename std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>::iterator itmi=
            edges_to_split.find(vp);
          if (itmi!=edges_to_split.end()){
            mi=get_vertex_index_on_edge(vp.first,vp.second,itmi->second);
            break;
          }
        }
        CGAL_assertion(k!=3);
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[k],mi,indices[(k+2)%3]));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(mi,indices[(k+1)%3],indices[(k+2)%3]));
        break;
      }
    }
  }

  //function used to detect edges that have two contour points as endpoints and
  //that are inside the domain. These edges are possible locations where non-manifold
  //connections can happen
  void split_non_manifold_edges_between_contour_points(double topz,bool last_run){
    bool first_run=previous_layer_incontour_edges==NULL;
    
    incontour_edge_ds.clear();
    Unique_incontour_edges* current_layer_incontour_edges=last_run?NULL:new Unique_incontour_edges();
    
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it ) {
      if(it->info().type == CellInfo3::EDGE_EDGE && !it->info().volume){
        int top_i[2];
        int bottom_i[2];
        
        for (int i=0,* t_ptr=top_i,* b_ptr=bottom_i;i<4;++i)
          if ( it->vertex(i)->point().z()==topz ) *t_ptr++=i;
          else *b_ptr++=i;
        CGAL_assertion( it->vertex(top_i[0])->point().z()==topz &&
                        it->vertex(top_i[1])->point().z()==topz );
        CGAL_assertion(it->info().f0->vertex(it->info().i0)->info().z==topz);
        CGAL_assertion(it->info().f1->vertex(it->info().i1)->info().z==it->vertex(bottom_i[0])->point().z());
        
        //Looking for an edge in the top layer
        if ( it->vertex(top_i[0])->info().v->info().on_contour && 
             it->vertex(top_i[1])->info().v->info().on_contour &&
             !it->info().f0->is_constrained(it->info().i0)     && 
             it->info().f0->info().in_domain() )
        {
          CGAL_assertion( it->info().f0->neighbor(it->info().i0)->info().in_domain() );
          std::pair<Vertex_handle_2,Vertex_handle_2> vertex_pair=make_sorted_pair(it->vertex(top_i[0])->info().v,it->vertex(top_i[1])->info().v);
          if ( first_run || 
               previous_layer_incontour_edges->find(vertex_pair)!=previous_layer_incontour_edges->end() )
          {
            std::pair<Vertex_handle_3,Vertex_handle_3> vertex_pair_3(it->vertex(top_i[0]),it->vertex(top_i[1]));
            if (it->neighbor(bottom_i[0])->info().volume)
              incontour_edge_ds.set_non_planar_facet(vertex_pair_3,it,bottom_i[0],bottom_i[1]);
            if (it->neighbor(bottom_i[1])->info().volume)
              incontour_edge_ds.set_non_planar_facet(vertex_pair_3,it,bottom_i[1],bottom_i[0]);
          }
        }
        
        //Looking for an edge in the bottom layer
        if ( it->vertex(bottom_i[0])->info().v->info().on_contour && 
             it->vertex(bottom_i[1])->info().v->info().on_contour &&
             !it->info().f1->is_constrained(it->info().i1) &&
             it->info().f1->info().in_domain()  )
        {
          if (!last_run){
            std::pair<Vertex_handle_2,Vertex_handle_2> vertex_pair=
              make_sorted_pair(it->vertex(bottom_i[0])->info().v,it->vertex(bottom_i[1])->info().v);
            current_layer_incontour_edges->insert(vertex_pair);
          }
          else{
            std::pair<Vertex_handle_3,Vertex_handle_3> vertex_pair_3(it->vertex(bottom_i[0]),it->vertex(bottom_i[1]));
            if (it->neighbor(top_i[0])->info().volume)
              incontour_edge_ds.set_non_planar_facet(vertex_pair_3,it,top_i[0],top_i[1]);
            if (it->neighbor(top_i[1])->info().volume)
              incontour_edge_ds.set_non_planar_facet(vertex_pair_3,it,top_i[1],top_i[0]);
          }
        }
      }
    }
    
    // split incontour edges from planar facets that are non-manifold (that might not be the most efficient.
    // if we want to do better, upon insertion of an edge in current_layer_incontour_edges, take the two planar facets incident to the edge
    // and insert them in a list to be considered for the next layer.
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){
      if( (!last_run && it->info().type!=CellInfo3::TOP_TRIANGLE) ||
          (last_run && it->info().type==CellInfo3::EDGE_EDGE) ||
          it->info().volume ||
          !it->info().f0->info().in_domain()
      ) continue;

      int i=it->info().i0;
      int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};

      std::list<std::pair<Vertex_handle_3,Vertex_handle_3> > edges_to_split;
      Face_handle_2 fh=it->info().f0;
      
      bool has_only_constrained=true; //handle isolated contours made of only 3 edges
      for (int k=0;k<3;++k){
        if (!fh->is_constrained(k)) has_only_constrained=false;
        
        if (fh->vertex((k+1)%3)->info().on_contour && 
            fh->vertex((k+2)%3)->info().on_contour && 
            !fh->is_constrained(k))
        {
          int i1=0,i2=0;
          int j=0;
          for (;j<3;++j)
            if( fh->vertex((k+1)%3)==it->vertex(indices[j])->info().v ){
              i1=indices[j];
              break;
            }
          CGAL_assertion(j!=3);
          for (j=0;j<3;++j)
            if( fh->vertex((k+2)%3)==it->vertex(indices[j])->info().v ){
              i2=indices[j];
              break;
            }
          CGAL_assertion(j!=3);
          
          std::pair<Vertex_handle_2,Vertex_handle_2> vertex_pair=make_sorted_pair(it->vertex(i1)->info().v,it->vertex(i2)->info().v);
          if( first_run || (last_run && it->info().type==CellInfo3::BOTTOM_TRIANGLE) ||
              previous_layer_incontour_edges->find(vertex_pair)!=previous_layer_incontour_edges->end() )
          {
            std::pair<Vertex_handle_3,Vertex_handle_3> vertex_pair_3=make_sorted_pair(it->vertex(i1),it->vertex(i2));
            edges_to_split.push_back(vertex_pair_3);
          }
        }
      }
      


      //Write planar facets if non-manifold edges have been found
      if (!edges_to_split.empty()){
        planar_facets_already_handled.insert(it->info().f0); //do this to avoid printing those facets in save_surface_layer
        incontour_edge_ds.set_planar_facet(edges_to_split,it);
      }
      else
        if ( has_only_constrained &&( first_run || last_run || previous_layer_printed_planar_facets.find(it->info().f0)!=previous_layer_printed_planar_facets.end() ) ){
          incontour_edge_ds.triangular_isolated_contours.push_back(it);  //we cannot "render" it now, we must store it in order the central
                                                                         //point to have a correct z value
          planar_facets_already_handled.insert(it->info().f0); //do this to avoid printing those facets in save_surface_layer
        }
    }
    
    if ( !first_run ) delete previous_layer_incontour_edges;
    previous_layer_incontour_edges=current_layer_incontour_edges;
  }
  #endif
  
  void check_removed_cells(){
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it ) {
      //a edge_edge tetrahedron that is in the volume must have at least one neighbor in the volume
      if(it->info().type == CellInfo3::EDGE_EDGE && it->info().volume){
        if (it->info().volume){
          int i=0;
          for (;i<4;++i){
            Cell_handle_3 neighbor = it->neighbor(i);
            if (neighbor->info().volume) break;
          }
          
          if (i==4){
            std::cerr << "ERROR" << std::endl;
            for (i=0;i<4;++i)
              std::cerr << it->info().out_face_status[i];
            std::cerr << std::endl;
            for (i=0;i<4;++i){
              Facet_3 mirror=delaunay_3.mirror_facet(Facet_3(it,i));
              std::cerr << mirror.first->info().out_face_status[mirror.second];
            }
            std::cerr << std::endl;
            for (i=0;i<4;++i){
              Facet_3 mirror=delaunay_3.mirror_facet(Facet_3(it,i));
              std::cerr << (mirror.first->info().type==CellInfo3::EDGE_EDGE?"e":mirror.first->info().f0->info().in_domain()?"t-in":"t-out");
            }
            std::cerr << std::endl;
            for (i=0;i<4;++i)
              std::cerr << (it->vertex(0)->point().z()==it->vertex(i)->point().z());
            std::cerr << std::endl;  
            for (i=0;i<4;++i)
              std::cerr << &(*it->neighbor(i)) << " ";
            std::cerr << std::endl;            
              //std::cerr << it->vertex(i)->point() << std::endl;
          }
        }
      }
    }
  }

  void write_layer(std::ostream& os, const CDT2& cdt,bool in_domain=true)//DEBUG
  {
    double z =  cdt.finite_vertices_begin()->info().z;
    int number_of_domain_faces = 0;
    for(Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it){
      if(it->info().in_domain()==in_domain){
        ++number_of_domain_faces;
      }
    }
    os << "OFF\n" << cdt.number_of_vertices() << " " << number_of_domain_faces << " 0\n";
    for(Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); ++it){
      os << it->point() << " " << z << std::endl;
    }
    for(Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it){
      if(it->info().in_domain()==in_domain){
         os << "3 " << it->vertex(0)->info().index << " "  << it->vertex(1)->info().index << " "  << it->vertex(2)->info().index << std::endl;
      }
    }
  }

  
  void write_projected_circumcenters(std::ostream& os, const CDT2& cdt,double other_z)//DEBUG
  {
    for(Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it){
      if(!it->info().in_domain()){
        it->vertex(0)->point()=it->vertex(0)->point();
        Point_2 p=circumcenter(it);
        os << 2 << "\n" <<  Point_3(p.x(),p.y(),it->vertex(0)->info().z) << "\n" << Point_3(p.x(),p.y(),other_z) << "\n" << std::endl;
      }
    }
  }
  
  void write_voronoi(std::ostream& os, const CDT2& cdt)//DEBUG
  {
    double z =  cdt.finite_vertices_begin()->info().z;
    for(Finite_edges_iterator it = cdt.finite_edges_begin(); it != cdt.finite_edges_end(); ++it){
      Face_handle_2 fh, nh;
      int fi;
      boost::tie(fh,fi) = *it;
      nh = fh->neighbor(fi);
      
      //    if(fh->info().in_domain && nh->info().in_domain()){
      if(! cdt.is_infinite(fh) && ! cdt.is_infinite(nh)){
        os << "2 " << cdt.circumcenter(fh) << " " << z+1 <<  " "  << cdt.circumcenter(nh) << " " << z+1 << std::endl;;  
      }
    }  
  }
  
  void write_cells(std::ostream& os) //DEBUG
  { 
    int fc = 0; 
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){ 
      if(it->info().volume){ 
        for(int i=0; i < 4; i++){ 
          if( (it->info().type == CellInfo3::EDGE_EDGE) || (it->info().i0 != i)){ 
            if(! it->neighbor(i)->info().volume){
              ++fc; 
            } 
          } 
        } 
      } 
    } 
    os << "OFF\n" << delaunay_3.number_of_vertices() << " " << fc << " 0\n"; 

    int vc = 0; 
    std::map<Vertex_handle_3,int> vindices;
    for(Vertex_iterator_3 it = delaunay_3.finite_vertices_begin(); it != delaunay_3.finite_vertices_end(); ++it){ 
      vindices[it] = vc++; 
      os << it->point() << std::endl; 
    } 
     
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){    
      if(it->info().volume){  
        for(int i=0; i < 4; i++){ 
          if( (it->info().type == CellInfo3::EDGE_EDGE) || (it->info().i0 != i)){ 
            if(! it->neighbor(i)->info().volume){ 
              int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4}; 
              if (i%2!=0) 
                std::swap(indices[1],indices[2]); 

              os << "3 "  
                 << vindices[it->vertex(indices[0]) ]<<  " "   
                 << vindices[it->vertex(indices[1]) ] <<  " "   
                 << vindices[it->vertex(indices[2]) ] <<  std::endl; 
            }
          }
        }
      }
    }
  }

  void write_deleted_tetra()  //DEBUG
  { 
    static int nb = 0; 
            
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){ 
      if(!it->info().volume){
//        if (it->info().type == CellInfo3::EDGE_EDGE) continue;
//        if (!it->info().f0->info().in_domain()) continue;
        std::stringstream ss;
        ss << "deleted/";
        if (it->info().type == CellInfo3::EDGE_EDGE) ss << "T22-";
        else ss << "T31-";
        if (++nb<10) ss<<"0";
        if (nb<100) ss<<"0";
        ss << nb << ".off";
        std::string fname;
        ss >> fname;
        std::ofstream output(fname.c_str());
        output << "OFF 4 4 0\n";
        for(int i=0; i < 4; i++) output << it->vertex(i)->point() << "\n";
        output << "3 1 3 0\n3 2 3 1\n3 0 3 2\n3 2 1 0\n";
      }
    }
  }


  Point_3 get_point_3(const Point_3& coords) const {
    return Point_3( coords[x_index],coords[(x_index+1)%3],coords[(x_index+2)%3] );
  }  

  Point_3 get_point_3(double a,double b, double c) const {
    double coords[3]={a,b,c};
    return Point_3( coords[x_index],coords[(x_index+1)%3],coords[(x_index+2)%3] );
  }
  
  int get_vertex_index(Vertex_handle_2 vh){ //OUTPUT
    if (vh->info().poly_index==-1){
      double z=vh->info().z;
      slice_writer_ptr->surface_point_push_back( get_point_3(vh->point().x(),vh->point().y(),z) );
      vh->info().poly_index=slice_writer_ptr->last_point_index();
      }
    return vh->info().poly_index;    
  }
  
  #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
  Point_3 get_additional_point(Vertex_handle_2 vh){ //OUTPUT
    double z=vh->info().z;
    return get_point_3(vh->point().x(),vh->point().y(),z);
  }
  #endif

  int get_vertex_index(Vertex_handle_3 vh){ //OUTPUT
    return get_vertex_index(vh->info().v);
  }  
  
  #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  int get_vertex_index_on_edge(Vertex_handle_3 ve1,Vertex_handle_3 ve2,int mp_vector_index){//OUTPUT
    if (incontour_edge_ds.midpoint_indices[mp_vector_index]==-1){
      Point_2 mp2=midpoint(ve1->info().v,ve2->info().v);
      double z = ve1->point().z();
      slice_writer_ptr->surface_point_push_back( get_point_3(mp2.x(),mp2.y(),z) );
      incontour_edge_ds.midpoint_indices[mp_vector_index]=slice_writer_ptr->last_point_index();
    }
    return incontour_edge_ds.midpoint_indices[mp_vector_index];
  }
  #endif

  void write_first_planar_facets(const CDT2& top_layer,bool is_increasing_z) //OUTPUT
  {
    for (typename CDT2::Finite_faces_iterator fit=top_layer.finite_faces_begin(),
                                              end=top_layer.finite_faces_end();fit!=end;++fit)
    {
      if ( fit->info().in_domain() ){
        int ind[3] = {0,1,2};
        if (is_increasing_z) std::swap(ind[0],ind[1]);
        slice_writer_ptr->surface_indices_push_back(
            cpp0x::make_tuple(
              get_vertex_index(fit->vertex(ind[0])),
              get_vertex_index(fit->vertex(ind[1])),
              get_vertex_index(fit->vertex(ind[2])))
        );
      }
    }
  }
  
  void add_last_bottom_layer_facets_to_surface(CDT2& last_layer,bool is_increasing_z){ //OUTPUT
    
    #ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
    reinit_non_contour_vertices_poly_index(last_layer);
    #endif
    for (typename CDT2::Finite_faces_iterator fit=last_layer.finite_faces_begin(),
                                    end=last_layer.finite_faces_end();fit!=end;++fit)
    {
      if ( fit->info().in_domain() ){
        int ind[3]={0,1,2};
        if (!is_increasing_z) std::swap(ind[0],ind[1]);
        slice_writer_ptr->surface_indices_push_back(
          cpp0x::make_tuple(
            get_vertex_index( fit->vertex(ind[0]) ),
            get_vertex_index( fit->vertex(ind[1]) ),
            get_vertex_index( fit->vertex(ind[2]) ) )
        );
      }
    }
    slice_writer_ptr->finalize_layer();
  }
  
  struct Cell_int_map_inserter{
    Cell_int_map_inserter(std::map<Cell_handle_3,int>& p_map):m_map(p_map){}
    std::map<Cell_handle_3,int>& m_map;
    void operator()(Cell_handle_3 cell){
      m_map.insert(std::make_pair(cell,-1));
    }
  };
  
  #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  //same as cell->index(vertex) except that it returns -1 if the vertex
  //is not part of the cell
  int index_of_vertex_in_cell(Cell_handle_3 cell,Vertex_handle_3 vertex){
    if (cell->vertex(0)==vertex) return 0;
    if (cell->vertex(1)==vertex) return 1;
    if (cell->vertex(2)==vertex) return 2;
    if (cell->vertex(3)==vertex) return 3;
    return -1;
  }

  Point_3 build_new_point_on_edge(const Point_3& p1,const Point_3& p2){
      To_exact to_exact;
      To_input to_input;
      return to_input(
        Exact_kernel::Construct_barycenter_3() (to_exact(p1),to_exact(0.9),to_exact(p2))
      );
  }
  
  Point_3 build_new_point_in_facet(Vertex_handle_3 v0,Vertex_handle_3 v1,Vertex_handle_3 v2){
    To_exact to_exact;
    To_input to_input;
    typename Exact_kernel::Point_3 ev0=to_exact(v0->point());
    typename Exact_kernel::Vector_3 V1=to_exact(v1->point())-ev0;
    typename Exact_kernel::Vector_3 V2=to_exact(v2->point())-ev0;
    return to_input( ev0 + to_exact(0.1)* (V1+V2) );
  }
  
  template <class Map>
  int get_vertex_index_on_edge(Vertex_handle_3 v1,Vertex_handle_3 v2,Map& new_vertex_on_edge){
    typedef typename Map::iterator Iterator;
    std::pair<Iterator,bool> ires=
      new_vertex_on_edge.insert(
        std::make_pair( std::make_pair(v1,v2),slice_writer_ptr->last_point_index()+1 )
      );
    if (ires.second){
      const Point_3 p=build_new_point_on_edge(v1->point(),v2->point());
      slice_writer_ptr->surface_point_push_back( get_point_3(p) );
    }
    return ires.first->second;
  }


  cpp0x::tuple<Vertex_handle_3,Vertex_handle_3,Vertex_handle_3>
  make_key(Vertex_handle_3 a, Vertex_handle_3 b, Vertex_handle_3 c) const{
    if (a<b) std::swap(a,b);
    return cpp0x::make_tuple(a,b,c);
  }

  //point opposite to a non-manifold edge inside a facet
  template <class Map>
  int get_vertex_index_opposite_to_nm_edge(
      Vertex_handle_3 e1,Vertex_handle_3 e2, Vertex_handle_3 other,
      Map& new_opposite_nm_edge_vertex)
  {
    typedef typename Map::iterator Iterator;
    std::pair<Iterator,bool> ires=
      new_opposite_nm_edge_vertex.insert(
        std::make_pair( make_key(e1,e2,other),slice_writer_ptr->last_point_index()+1 )
      );
    if (ires.second){
      const Point_3 p=build_new_point_on_edge( midpoint(e1->point(),e2->point()),other->point() );
      slice_writer_ptr->surface_point_push_back( get_point_3(p) );
    }
    return ires.first->second;
  }
  
  //the point computed is the intersection of the segment opposite to i1 with the segment opposite to i2 in facet(cell,facet_index)
  //vertex_index is just used for the key
  template <class Map>
  int get_vertex_index_in_facet(Cell_handle_3& cell,int facet_index,int vertex_index,int i1, int i2,Map& new_vertex_in_facet){
    typedef typename Map::iterator Iterator;
    std::pair<typename Map::iterator,bool> ires=
      new_vertex_in_facet.insert(
        std::make_pair( std::make_pair(facet_index,vertex_index),slice_writer_ptr->last_point_index()+1 )
      );
    if (ires.second){
      CGAL_assertion(i1+i2+vertex_index+facet_index==6);
      const Point_3 p=build_new_point_in_facet(cell->vertex(vertex_index),cell->vertex(i1),cell->vertex(i2));
      slice_writer_ptr->surface_point_push_back( get_point_3(p) );
    }
    return ires.first->second;    
  }
  
  void display_surface_facet_with_a_nm_edge_and_manifold_endpoints(
    Cell_handle_3 cell,
    int e1,int e2,int other,
    std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>& new_vertex_on_edge,
    std::map<cpp0x::tuple<Vertex_handle_3,Vertex_handle_3,Vertex_handle_3>,int>& new_opposite_nm_edge_vertex)
  {
    int middle=get_vertex_index_opposite_to_nm_edge(cell->vertex(e1),cell->vertex(e2),cell->vertex(other),new_opposite_nm_edge_vertex);
    if ( cell->info().nm_vertices(other) ){
      int other_e1=get_vertex_index_on_edge(cell->vertex(other),cell->vertex(e1),new_vertex_on_edge);
      int other_e2=get_vertex_index_on_edge(cell->vertex(other),cell->vertex(e2),new_vertex_on_edge);
      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(e1)),middle,other_e1) );
      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(e2)),other_e2,middle) );
      
      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(middle,other_e2,other_e1) );
    }
    else{
      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(e1)),middle,get_vertex_index(cell->vertex(other))) );
      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(e2)),get_vertex_index(cell->vertex(other)),middle) );
    }
  }
  
  void display_surface_facet_in_cells_with_non_manifold_features(
    Cell_handle_3 cell,
    std::set<Face_handle_2>& previous_layer_printed_planar_facets,
    std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int>& new_vertex_on_edge,
    std::map<cpp0x::tuple<Vertex_handle_3,Vertex_handle_3,Vertex_handle_3>,int>& new_opposite_nm_edge_vertex)
  {
    CGAL_assertion(!cell->info().volume);

    //map storing vertex indices that are inside a facet (intersection of two edges in a boundary facet due to two non-manifold edge adjacent)
    std::map<std::pair<int,int>,int> new_vertex_in_facet;
    
    //1) Handle boundary facets of a cell with non-manifold features
    for (int i=0;i<4;++i){
      if ( cell->neighbor(i)->info().volume ){
        int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
        if (i%2==0) std::swap(indices[1],indices[2]);
        
        
        //special case to handle a non-manifold edge with two manifold endpoints
        if ( cell->info().has_a_nm_edge_with_manifold_vertices ){
          if ( cell->info().nm_edge(indices[0],indices[1]) && !cell->info().nm_vertices(indices[0]) && !cell->info().nm_vertices(indices[1]) )
          {
            display_surface_facet_with_a_nm_edge_and_manifold_endpoints(cell,indices[0],indices[1],indices[2],new_vertex_on_edge,new_opposite_nm_edge_vertex);
            continue;
          }
          if ( cell->info().nm_edge(indices[0],indices[2]) && !cell->info().nm_vertices(indices[0]) && !cell->info().nm_vertices(indices[2]) )
          {
            display_surface_facet_with_a_nm_edge_and_manifold_endpoints(cell,indices[2],indices[0],indices[1],new_vertex_on_edge,new_opposite_nm_edge_vertex);
            //special_function_call
            continue;
          }
          if ( cell->info().nm_edge(indices[1],indices[2]) && !cell->info().nm_vertices(indices[1]) && !cell->info().nm_vertices(indices[2]) )
          {
            display_surface_facet_with_a_nm_edge_and_manifold_endpoints(cell,indices[1],indices[2],indices[0],new_vertex_on_edge,new_opposite_nm_edge_vertex);
            continue;
          }
        }
        
        
        std::vector<int> ordered_indices;

        //turn around the facet and consider each vertex. We collect boundary points around the current vertex.
        //in ordered_indices we put points of the convex facet defined. In case of
        //a non-manifold edge with a manifold vertex we have a extra face displayed independantly
        for (int k=0;k<3;++k){
          const int i_curr = indices[k];
          const int i_before = indices[(k+2)%3];
          const int i_next = indices[(k+1)%3];
          
          //check if first edge is non-manifold
          if ( cell->info().nm_edge(i_curr,i_before) ){
            //check if the second edge is non-manifold
            if ( cell->info().nm_edge(i_curr,i_next) ){
              CGAL_assertion( cell->info().nm_vertices(i_curr) );//the vertex must be non-manifold
              ordered_indices.push_back( get_vertex_index_in_facet(cell,i,i_curr,i_before,i_next,new_vertex_in_facet) );
            }
            else{
              if ( cell->info().nm_vertices(i_curr) )
                ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_next),new_vertex_on_edge) );
              else
                ordered_indices.push_back(get_vertex_index(cell->vertex(i_curr)));
            }
            continue;
          }
          else{
            //check if the second edge is non-manifold
            if ( cell->info().nm_edge(i_curr,i_next) ){
              if ( cell->info().nm_vertices(i_curr) )
                ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_before),new_vertex_on_edge) );
              else
                ordered_indices.push_back(get_vertex_index(cell->vertex(i_curr)));
              continue;
            }
          }
          
          if ( cell->info().nm_vertices(i_curr) ){
            ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_before),new_vertex_on_edge) );
            ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_next),new_vertex_on_edge) );
          }
          else
            ordered_indices.push_back(get_vertex_index(cell->vertex(i_curr)));
        }
        
        std::size_t k=ordered_indices.size()-1;
        for (std::size_t i=1;i<k;++i)
          slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(ordered_indices[i],ordered_indices[i+1],ordered_indices[0]));
      }
      else{
        //case of a planar facet
        
        if ( cell->info().type!=CellInfo3::EDGE_EDGE && cell->info().i0==i && cell->info().f0->info().in_domain() )
        {
          CGAL_assertion( cell->info().f0->info().in_domain() );
          //collect facets that have been printed to detect non-manifold facet at the next step
          previous_layer_printed_planar_facets.insert(cell->info().f0);
          
          int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
          if (i%2==0) std::swap(indices[1],indices[2]);
          
          //special case of a isolated triangle
          if ( cell->info().isolated_triangle ){
            Point_3 pt=centroid( cell->vertex(indices[0])->point(), cell->vertex(indices[1])->point(), cell->vertex(indices[2])->point() );
            pt=build_new_point_on_edge(pt,cell->vertex(i)->point());
            slice_writer_ptr->surface_point_push_back( get_point_3(pt) );
            int central=slice_writer_ptr->last_point_index();
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( get_vertex_index(cell->vertex(indices[0])), get_vertex_index(cell->vertex(indices[1])), central ) );
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( get_vertex_index(cell->vertex(indices[1])), get_vertex_index(cell->vertex(indices[2])), central ) );
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( get_vertex_index(cell->vertex(indices[2])), get_vertex_index(cell->vertex(indices[0])), central ) );
            continue;
          }
          
          
          //special case to handle a non-manifold edge with two manifold endpoints
          if ( cell->info().has_a_nm_edge_with_manifold_vertices ){
            std::vector< int > edges; //store the first vertex (in ccw order) of a non-manifold edge with manifold endpoints

            if ( cell->info().nm_edge(indices[0],indices[1]) && !cell->info().nm_vertices(indices[0]) && !cell->info().nm_vertices(indices[1]) )
            { edges.push_back(0); }
            if ( cell->info().nm_edge(indices[1],indices[2]) && !cell->info().nm_vertices(indices[1]) && !cell->info().nm_vertices(indices[2]) )
            { edges.push_back(1); }
            if ( cell->info().nm_edge(indices[0],indices[2]) && !cell->info().nm_vertices(indices[0]) && !cell->info().nm_vertices(indices[2]) )
            { edges.push_back(2); }
            
            #warning think about what can happen if an edge incident to i0 is non-manifold. is it possible?
            
            switch (edges.size() ){
              case 3:
              {
                int e_01=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[0]),cell->vertex(indices[1]),cell->vertex(i),new_opposite_nm_edge_vertex);
                int e_12=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[1]),cell->vertex(indices[2]),cell->vertex(i),new_opposite_nm_edge_vertex);
                int e_20=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[2]),cell->vertex(indices[0]),cell->vertex(i),new_opposite_nm_edge_vertex);
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(e_01,e_12,e_20) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_12, e_01, get_vertex_index(cell->vertex(indices[1])) ) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_20, e_12, get_vertex_index(cell->vertex(indices[2])) ) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_01, e_20, get_vertex_index(cell->vertex(indices[0])) ) );
                continue;
              }
              case 2:
              {
                int other=3-edges[1]-edges[0];//the only edge which is manifold
                int e_p1=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[(other+1)%3]),cell->vertex(indices[(other+2)%3]),cell->vertex(i),new_opposite_nm_edge_vertex);
                int e_p2=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[(other+2)%3]),cell->vertex(indices[other]),cell->vertex(i),new_opposite_nm_edge_vertex);
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_p1, e_p2, get_vertex_index(cell->vertex(indices[(other+1)%3])) ) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_p2, e_p1, get_vertex_index(cell->vertex(indices[(other+2)%3])) ) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( e_p2, get_vertex_index(cell->vertex(indices[other])) , get_vertex_index(cell->vertex(indices[(other+1)%3])) ) );
                continue;
              }
              case 1:
              {
                //edges[0] is the only edge which is non-manifold with manifold endpoints
                int middle=get_vertex_index_opposite_to_nm_edge(cell->vertex(indices[edges[0]]),cell->vertex(indices[(edges[0]+1)%3]),cell->vertex(i),new_opposite_nm_edge_vertex);
                if (cell->info().nm_vertices( indices[(edges[0]+2)%3] )){
                  if ( cell->info().nm_facet() ){
                    int other=get_vertex_index_on_edge(cell->vertex( indices[(edges[0]+2)%3] ), cell->vertex(i), new_vertex_on_edge);
                    slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, other, get_vertex_index(cell->vertex(indices[edges[0]])) ) );
                    slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+1)%3])), other) );
                  }
                  else{
                    if (cell->info().nm_edge(indices[(edges[0]+1)%3],indices[(edges[0]+2)%3])){
                      int other=
                        (cell->info().nm_edge(indices[(edges[0]+2)%3],indices[edges[0]]))?
                          get_vertex_index_in_facet(cell,i,indices[(edges[0]+2)%3],indices[(edges[0]+1)%3],indices[edges[0]],new_vertex_in_facet)
                        :
                          get_vertex_index_on_edge(cell->vertex( indices[(edges[0]+2)%3] ), cell->vertex( indices[edges[0]]       ), new_vertex_on_edge);
                      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, other, get_vertex_index(cell->vertex(indices[edges[0]])) ) );
                      slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+1)%3])), other) );
                    }
                    else{
                      if (!cell->info().nm_edge(indices[(edges[0]+2)%3],indices[edges[0]])){
                        int other=get_vertex_index_on_edge(cell->vertex( indices[(edges[0]+2)%3] ), cell->vertex( indices[(edges[0]+1)%3] ), new_vertex_on_edge);
                        slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, other, get_vertex_index(cell->vertex(indices[edges[0]])) ) );
                        slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+1)%3])), other) );                        
                      }
                      else{
                        int e_p20=get_vertex_index_on_edge(cell->vertex( indices[(edges[0]+2)%3] ), cell->vertex( indices[edges[0]]       ), new_vertex_on_edge);
                        int e_p21=get_vertex_index_on_edge(cell->vertex( indices[(edges[0]+2)%3] ), cell->vertex( indices[(edges[0]+1)%3] ), new_vertex_on_edge);
                        slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, e_p20, get_vertex_index(cell->vertex(indices[edges[0]])) ) );
                        slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+1)%3])), e_p21 ) );
                        slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, e_p21, e_p20 ) );                                          
                      }
                    }
                  }
                }
                else{
                  slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+2)%3])), get_vertex_index(cell->vertex(indices[edges[0]])) ) );
                  slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple( middle, get_vertex_index(cell->vertex(indices[(edges[0]+1)%3])), get_vertex_index(cell->vertex(indices[(edges[0]+2)%3])) ) );
                }
                continue;
              }
              default:
                CGAL_assertion( edges.size()==0 );
            }
          }
          
          //case when no non-manifold edge with manifold vertices has been found
          std::vector<int> ordered_indices;
          
          for (int k=0;k<3;++k){
            const int i_curr = indices[k];
            const int i_before = indices[(k+2)%3];
            const int i_next = indices[(k+1)%3];
            
            //check is the edge(i,i_curr) is non-manifold edge 
            if ( cell->info().nm_edge(i_curr,i) ){
              ordered_indices.push_back( get_vertex_index_in_facet(cell,i_next,i_curr,i,i_before,new_vertex_in_facet) );
              ordered_indices.push_back( get_vertex_index_in_facet(cell,i_before,i_curr,i,i_next,new_vertex_in_facet) );
            }
            else{
              if ( cell->info().nm_vertices(i_curr) ){
                if ( cell->info().nm_facet() )
                  ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i),new_vertex_on_edge) );
                else{
                  if ( cell->info().nm_edge(i_curr,i_before) ){
                    if ( cell->info().nm_edge(i_curr,i_next) ){
                      ordered_indices.push_back( get_vertex_index_in_facet(cell,i,i_curr,i_before,i_next,new_vertex_in_facet) );
                    }
                    else{
                      ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_next),new_vertex_on_edge) );
                    }
                  }
                  else
                  {
                    if ( cell->info().nm_edge(i_curr,i_next) ){
                      ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_before),new_vertex_on_edge) );
                    }
                    else{
                      ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_before),new_vertex_on_edge) );
                      ordered_indices.push_back( get_vertex_index_on_edge(cell->vertex(i_curr),cell->vertex(i_next),new_vertex_on_edge) );
                    }
                  }
                }
              }
              else
                ordered_indices.push_back( get_vertex_index(cell->vertex(i_curr) ) );
            }
          }
          std::size_t k=ordered_indices.size()-1;
          for (std::size_t i=1;i<k;++i)
            slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(ordered_indices[i],ordered_indices[i+1],ordered_indices[0]));
        }
      }
    }

    //2) handle parts added as non-manifold edges
    std::map< std::pair<int,int>, std::pair<int,int> > indices_for_edges;
    
    //consider each vertex and look for its non-manifold incident edges and add parts of the cut plane
    for (int i=0;i<4;++i)
    {
      if (cell->info().nm_vertices(i)){
        //cout the number of incident non-manifold edges
        std::vector<int> edges_nm;
        for(int k=1;k<4;++k){
          if ( cell->info().nm_edge(i,(i+k)%4) ){
            edges_nm.push_back((i+k)%4);
          }
        }
        
        switch(edges_nm.size()){
          case 0:
          {
            if (cell->info().type!=CellInfo3::EDGE_EDGE && i!=cell->info().i0 && cell->info().nm_facet()) continue;
            
            int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
            if (i%2!=0) std::swap(indices[1],indices[2]);
            slice_writer_ptr->surface_indices_push_back(
              cpp0x::make_tuple(
                get_vertex_index_on_edge(cell->vertex(i),cell->vertex(indices[0]),new_vertex_on_edge),
                get_vertex_index_on_edge(cell->vertex(i),cell->vertex(indices[1]),new_vertex_on_edge),
                get_vertex_index_on_edge(cell->vertex(i),cell->vertex(indices[2]),new_vertex_on_edge)
              )
            );
            break;              
          }
          case 1:
          {
            const int e1=i;
            const int e2=edges_nm[0];
            int indices[3]={(e2+1)%4,(e2+2)%4,(e2+3)%4};
            if (e2%2!=0) std::swap(indices[1],indices[2]);              
            for (int k=0;k<3;++k)
              if ( indices[k]==e1 ){
                  int i1=get_vertex_index_on_edge(cell->vertex(e1),cell->vertex(indices[(k+2)%3]),new_vertex_on_edge);
                  int i2=get_vertex_index_on_edge(cell->vertex(e1),cell->vertex(indices[(k+1)%3]),new_vertex_on_edge);
                  indices_for_edges[std::make_pair(e1,e2)]=std::make_pair(i1,i2);
                break;
              }
            break; 
          }
          case 2:
          {
            if ( cell->info().type!=CellInfo3::EDGE_EDGE && cell->info().nm_facet() ){
              if ( edges_nm[0]+edges_nm[1]+i+cell->info().i0== 6 ) continue; //the 3 vertices corresponding to non-manifold edges are in a non-manifold facet
            }
            
            const int e1=i;
            const int e2=edges_nm[0];
            
            int indices[3]={(e2+1)%4,(e2+2)%4,(e2+3)%4};
            if (e2%2!=0) std::swap(indices[1],indices[2]);
            for (int k=0;k<3;++k)
              if ( indices[k]==e1 ){
                  int i1,i2;
                  if ( indices[(k+2)%3] == edges_nm[1] ){
                    i1=get_vertex_index_in_facet(cell,indices[(k+1)%3],e1,e2,indices[(k+2)%3],new_vertex_in_facet);
                    i2=get_vertex_index_on_edge(cell->vertex(e1),cell->vertex(indices[(k+1)%3]),new_vertex_on_edge);
                  }
                  else{
                    i1=get_vertex_index_on_edge(cell->vertex(e1),cell->vertex(indices[(k+2)%3]),new_vertex_on_edge);
                    i2=get_vertex_index_in_facet(cell,indices[(k+2)%3],e1,e2,indices[(k+1)%3],new_vertex_in_facet);
                  }
                  indices_for_edges[std::make_pair(e1,e2)]=std::make_pair(i1,i2);
                  indices_for_edges[std::make_pair(e1,edges_nm[1])]=std::make_pair(i2,i1);
                break;
              }
            break; 
          }
          default:
          {
            CGAL_assertion(edges_nm.size()==3);
            int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
            if (i%2!=0) std::swap(indices[1],indices[2]);     
            int poly_indices[3]={-1,-1.-1};
            for (int k=0;k<3;++k)
              poly_indices[k]=get_vertex_index_in_facet(cell,indices[k],i,indices[(k+1)%3],indices[(k+2)%3],new_vertex_in_facet);
            for (int k=0;k<3;++k)
              indices_for_edges[std::make_pair(i,indices[k])]=std::make_pair(poly_indices[(k+2)%3],poly_indices[(k+1)%3]);
            //add the top part
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(poly_indices[0],poly_indices[1],poly_indices[2]) );
          }
        }
      }
    }

    //print a quad per non-manifold edge with two non-manifold vertices
    //print a triangle per non-manifold edge with one non-manifold vertex
    //print two triangle per non-manifold edge with manifold vertices
    for (int i=0;i<4;++i)
      for(int j=i+1;j<4;++j)
        if (
            (!cell->info().nm_facet() || i==cell->info().i0 || j==cell->info().i0)  //do not draw for edges of a non-manifold facet
            && cell->info().nm_edge(i,j) 
        ){
          if ( cell->info().nm_vertices(i) ){
            if ( cell->info().nm_vertices(j) ){
              typename std::map< std::pair<int,int>, std::pair<int,int> >::iterator it_ij=indices_for_edges.find(std::make_pair(i,j));
              CGAL_assertion( it_ij!=indices_for_edges.end() );
              typename std::map< std::pair<int,int>, std::pair<int,int> >::iterator it_ji=indices_for_edges.find(std::make_pair(j,i));
              CGAL_assertion( it_ji!=indices_for_edges.end() );
              
              std::pair<int,int>& indices_i=it_ij->second;
              std::pair<int,int>& indices_j=it_ji->second;
              slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(indices_j.first,indices_i.first,indices_i.second ) );
              slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(indices_i.first,indices_j.first,indices_j.second ) );            
            }
            else{
              typename std::map< std::pair<int,int>, std::pair<int,int> >::iterator it_ij=indices_for_edges.find(std::make_pair(i,j));
              CGAL_assertion( it_ij!=indices_for_edges.end() );

              std::pair<int,int>& indices_i=it_ij->second;
              slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(j)),indices_i.first,indices_i.second ) );
            }
          }
          else{
            if ( cell->info().nm_vertices(j) ){
              typename std::map< std::pair<int,int>, std::pair<int,int> >::iterator it_ji=indices_for_edges.find(std::make_pair(j,i));
              CGAL_assertion( it_ji!=indices_for_edges.end() );
              std::pair<int,int>& indices_j=it_ji->second;
              slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(i)),indices_j.first,indices_j.second ) );
            }
            else{
              //case of a non-manifold edge with manifold endpoints. Nothing need to be done if the facet is planar
              if ( cell->info().type==CellInfo3::EDGE_EDGE ){
                int k=0;
                int o1,o2;
                while ( (i+k+1)%4!=j ) ++k;
                o1=( i + 1 + (k+1)%3 )%4;
                o2=( i + 1 + (k+2)%3 )%4;
                if (i%2==0) std::swap(o1,o2);
                //the vertex i sees j,o1,o2 in ccw order
                CGAL_assertion(i+j+o1+o2==6);
                int ij1=get_vertex_index_opposite_to_nm_edge(cell->vertex(i),cell->vertex(j),cell->vertex(o1),new_opposite_nm_edge_vertex);
                int ij2=get_vertex_index_opposite_to_nm_edge(cell->vertex(i),cell->vertex(j),cell->vertex(o2),new_opposite_nm_edge_vertex);
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(i)),ij1,ij2) );
                slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(get_vertex_index(cell->vertex(j)),ij2,ij1) );
              }
            }
          }
        }

    //3) handle facets that are did not contribute to the boundary, but since we resolved a non-manifold situation
    //   parts of that facet should be displayed (those that are on the boundary of the cell).
    if ( cell->info().type==CellInfo3::EDGE_EDGE ){
      for (int i=0;i<4;++i){
        for (int j=i+1;j<4;++j)
          if ( cell->info().nm_edge(i,j) &&
               cell->vertex(i)->point().z()==cell->vertex(j)->point().z()
          )
          {
            //look whether in the neighbor whether vertex i is non-manifold
            int index_in_neigh=cell->neighbor(j)->index( cell->vertex(i) );
            if ( cell->info().nm_vertices(i) &&
                 !cell->neighbor(j)->info().volume && 
                 !cell->neighbor(j)->info().nm_vertices(index_in_neigh)  ){
              //CGAL_assertion( cell->vertex(i)->info().v->info().on_contour );
              int indices[3]={(j+1)%4,(j+2)%4,(j+3)%4};
              if (j%2!=0) std::swap(indices[1],indices[2]);
              int index_of_i = 0;
              while (indices[index_of_i]!=i) ++index_of_i;
              slice_writer_ptr->surface_indices_push_back( 
                cpp0x::make_tuple(
                  get_vertex_index( cell->vertex(i) ),
                  get_vertex_index_on_edge(cell->vertex(i),cell->vertex(indices[(index_of_i+1)%3]),new_vertex_on_edge),
                  get_vertex_index_on_edge(cell->vertex(i),cell->vertex(indices[(index_of_i+2)%3]),new_vertex_on_edge)
                )
              );    
            }
            
            //look whether in the neighbor whether vertex j is non-manifold
            index_in_neigh=cell->neighbor(i)->index( cell->vertex(j) );
            if ( cell->info().nm_vertices(j) &&
                 !cell->neighbor(i)->info().volume && 
                 !cell->neighbor(i)->info().nm_vertices(index_in_neigh)  ){
              //CGAL_assertion( cell->vertex(j)->info().v->info().on_contour );
              int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
              if (i%2!=0) std::swap(indices[1],indices[2]);
              int index_of_j = 0;
              while (indices[index_of_j]!=j) ++index_of_j;
              slice_writer_ptr->surface_indices_push_back( 
                cpp0x::make_tuple(
                  get_vertex_index( cell->vertex(j) ),
                  get_vertex_index_on_edge(cell->vertex(j),cell->vertex(indices[(index_of_j+1)%3]),new_vertex_on_edge),
                  get_vertex_index_on_edge(cell->vertex(j),cell->vertex(indices[(index_of_j+2)%3]),new_vertex_on_edge)
                )
              );    
            }
          }
      }
    }
    else{
      if ( cell->info().nm_facet() ){
        const int i0=cell->info().i0;
        for (int i=0;i<3;++i){
          Cell_handle_3 neighbor=cell->neighbor( (i0+i+1)%4 );
          if ( !neighbor->info().volume && 
               !neighbor->info().nm_facet() &&
               cell->info().nm_vertices((i0+1+(i+1)%3)%4) && cell->info().nm_vertices((i0+1+(i+2)%3)%4)
          ){
            if ( neighbor->info().non_manifold_features )
            {
              int i1 = neighbor->index( cell->vertex( (i0+1+(i+1)%3)%4 ) );
              int i2 = neighbor->index( cell->vertex( (i0+1+(i+2)%3)%4 ) );
              CGAL_assertion( (i0+i+1)%4+i0+(i0+1+(i+1)%3)%4+(i0+1+(i+2)%3)%4==6 );
              if ( neighbor->info().nm_edge(i1,i2) ) continue;
            }
            
            int i1,i2;
            //display a quad
            const int vi=(i0+i+1)%4;
            for (int k=0;k<3;++k)
              if ( (vi+k+1)%4==i0 ){
                i1=( vi+1+(k+1)%3 )%4;
                i2=( vi+1+(k+2)%3 )%4;
                break;
              }
            if (vi%2!=0) std::swap(i1,i2);

            CGAL_assertion(vi+i0+i1+i2==6);

            int i3=get_vertex_index_on_edge(cell->vertex(i2),cell->vertex(i0),new_vertex_on_edge);
            int i4=get_vertex_index_on_edge(cell->vertex(i1),cell->vertex(i0),new_vertex_on_edge);
            i1=get_vertex_index( cell->vertex(i1) );
            i2=get_vertex_index( cell->vertex(i2) );
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(i1,i2,i3) );  
            slice_writer_ptr->surface_indices_push_back( cpp0x::make_tuple(i1,i3,i4) );  
          }
        }
      }
    }
    
  }
  #endif
  
  void save_surface_layer(bool last_run){ //OUTPUT
   
    #ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
    
    #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    double top_z = top_ptr->finite_vertices_begin()->info().z;
    bool first_run = previous_bottom_incontour_nm_vertices==NULL;
    #endif

    //POSSIBLE OPTI: One way to directly compute incident facets is to make a pass on all cells and use a vector
    //  of size number_of_vertices containing the list of incident facets per vertex. We need to have an index in
    //  the info of the vertex for that
    
  //handle facets incident to non-manifold contour and non-contour points
    
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    //this map type is needed because a facet can have several non-manifold vertex to duplicate
    typedef std::map<Facet_3,cpp0x::array<int,3> >  Facet_to_nm_indices;
    //map non-planar facets that have at least one non-manifold vertex to duplicate to duplicated index of vertices or default otherwise
    Facet_to_nm_indices facets_already_handled;
    //map planar facets that have at least one non-manifold vertex to duplicate to duplicated index of vertices or default otherwise.
    Facet_to_nm_indices planar_facets_nm_point;    
    #endif

    
    //iterate over all vertices and look for their one neighbors. We look at the connected components
    //of tetrahedra that are outside the volume. If there is more than one such connected components, then
    //the vertex is non-manifold. For a contour vertex, among the connected components one contains facets incident to the contour edges 
    //the non-manifold vertex is indicent to. Other connected compondents are called umbrella and form a kind of pyramid
    //with apex the non-manifold vertex. For each umbrella, we need to duplicate the non manifold vertex and 
    //associate it to facets in the umbrella. The same procedure is working for non-contour points, but the connected component
    //that is left as is is the first one.
    for (typename DT3::Finite_vertices_iterator vit=delaunay_3.finite_vertices_begin(),
                                                vit_end=delaunay_3.finite_vertices_end();vit!=vit_end;++vit)
    {
      const bool vertex_on_contour = vit->info().v->info().on_contour;
      std::map<Cell_handle_3,int> incident_cells;
      //collect tetrahedra incident to vit
      delaunay_3.incident_cells(vit,
        boost::make_function_output_iterator(Cell_int_map_inserter(incident_cells)));
      
      //first pass to mark infinite cells incident to a T31;
      //in this pass we also collect the contour vertices incident to vit if vit is a contour vertex
      bool found_incident_vertices=!vertex_on_contour;
      Vertex_handle_2 vhandles[2];
      for (typename std::map<Cell_handle_3,int>::iterator it=incident_cells.begin(),
                                                          it_end=incident_cells.end();it!=it_end;++it)
      {
        if(!delaunay_3.is_infinite(it->first) && it->first->info().type != CellInfo3::EDGE_EDGE ){
          typename std::map<Cell_handle_3,int>::iterator rf=incident_cells.find(it->first->neighbor(it->first->info().i0));
          if( rf!=incident_cells.end() ){
            CGAL_assertion(it->first->vertex(it->first->info().i0) != Vertex_handle_3(vit));
            rf->second=-2;
            if ( !found_incident_vertices ){
              CGAL_assertion( it->first->info().f0->vertex(0)->info().z==vit->point().z() );
            //  turn around vit to find the two incident constrained edges 
              Face_handle_2 fh=it->first->info().f0;
              Vertex_handle_2 vh=vit->info().v;
              int vertex_index=fh->index(vh);
              int turn_index=(vertex_index+1)%3;
              int cindex=-1;
              CGAL_assertion_code(Face_handle_2 start=fh);
              do{
                if(fh->info().in_domain()){
                  if ( fh->is_constrained((vertex_index+1)%3) )
                    vhandles[++cindex]=fh->vertex( (vertex_index+2)%3 );
                  if ( fh->is_constrained((vertex_index+2)%3) )
                    vhandles[++cindex]=fh->vertex( (vertex_index+1)%3 );
                }
                
                Face_handle_2 new_fh=fh->neighbor(turn_index);
                turn_index=new_fh->index(fh);
                if ( new_fh->vertex((turn_index+1)%3)!=vh )   turn_index=(turn_index+1)%3;
                else turn_index=(turn_index+2)%3;
                fh=new_fh;
                vertex_index=fh->index(vh);
                CGAL_assertion(fh!=start|| cindex==1);
              }while(cindex!=1);
              found_incident_vertices=true;
            }
          }
        }
      }
      
      CGAL_precondition(found_incident_vertices);
      
      #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
      std::vector< std::vector<Cell_handle_3> > incident_cells_per_cc;
      #else
      std::vector< std::list<Facet_3> > boundary_facets_per_cc;
      std::vector< std::list<Facet_3> > boundary_planar_facets_per_cc;
      #endif
      //in practice it will always be 2 but put 5 to be on the safe side
      #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
      incident_cells_per_cc.reserve(5);
      #else
      boundary_facets_per_cc.reserve(5);
      boundary_planar_facets_per_cc.reserve(5);
      #endif
     
      #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
      //used to detect non-manifold non-planar edges: if we cut such an edge, it remains non-manifold
      std::set<Vertex_handle_3> incident_vertices;
      std::vector<Vertex_handle_3> non_manifold_edges;
      #endif
      
      //identify surface connected components containing vit 
      int cc_number=0;
      int cc_to_keep=-1;
      for (typename std::map<Cell_handle_3,int>::iterator it=incident_cells.begin(),
                                                          it_end=incident_cells.end();it!=it_end;++it)
      {
        if( it->first->info().volume || it->second!=-1 ) continue;

        #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
        //make a copy to test the current cc against the other
        std::set<Vertex_handle_3> last_cc_incident_vertices (incident_vertices);
        //start a new connected component
        incident_cells_per_cc.push_back(std::vector<Cell_handle_3>());
        #else
        boundary_facets_per_cc.push_back(std::list<Facet_3>());
        boundary_planar_facets_per_cc.push_back(std::list<Facet_3>());
        #endif
        
        std::list<Cell_handle_3> stack;
        stack.push_back(it->first);
        #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
        incident_cells_per_cc.back().push_back(it->first);
        #endif
        it->second=cc_number;
        bool has_neighbor_vertex=!vertex_on_contour;
        while ( !stack.empty() ){
          Cell_handle_3 cell=stack.back();
          stack.pop_back();
          
          #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
          for (int iv=0;iv<4;++iv){
            Vertex_handle_3 vh=cell->vertex(iv);
            if (vh==Vertex_handle_3(vit)) continue;
            typename std::set<Vertex_handle_3>::iterator it_value=last_cc_incident_vertices.find(vh);
            if ( it_value!=last_cc_incident_vertices.end() ){
//              std::cout << "Found a non-manifold edge\n";
//              std::cout << vit->point() << " " << vh->point() << std::endl;
              non_manifold_edges.push_back(vh);
              incident_vertices.erase(vh);
              last_cc_incident_vertices.erase(it_value);
            }
            else
              incident_vertices.insert(vh);
          }
          #endif
          
          int index_vit=cell->index(vit);
          for (int k=0;k<4;++k){
            if (k==index_vit) continue;
            Cell_handle_3 neighbor=cell->neighbor(k);
            if ( neighbor->info().volume ){
              CGAL_assertion(!cell->info().volume);
              #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
              boundary_facets_per_cc.back().push_back(Facet_3(cell,k));
              #endif
              //this is enough to look for one vertex only, the second will be inside too
              if ( !has_neighbor_vertex ){
                for (int j=1;j<4;++j){
                  Vertex_handle_3 vh3=cell->vertex((k+j)%4);
                  //we need to compare to vertex handles and cannot use contrained edges as the cells are
                  //outside the volume and we cannot necessarily access the corresponding T31
                  if (vh3->info().v==vhandles[0] || vh3->info().v==vhandles[1]){
                    has_neighbor_vertex=true;
                    break;
                  }
                }
              }
              continue;
            }
            typename std::map<Cell_handle_3,int>::iterator rf=incident_cells.find(neighbor);
            if (rf == incident_cells.end() || rf->second!=-1) continue;
            rf->second=cc_number;
            stack.push_back(neighbor);
            #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
            incident_cells_per_cc.back().push_back(neighbor);
            #endif
          }
          #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
          if (cell->info().type!=CellInfo3::EDGE_EDGE && cell->info().i0!=index_vit)
            boundary_planar_facets_per_cc.back().push_back(Facet_3(cell,cell->info().i0));
          #endif
        }
        if (vertex_on_contour && has_neighbor_vertex){
          CGAL_assertion( cc_to_keep==-1 );
          cc_to_keep=cc_number;
        }
        ++cc_number;
      }

      if (cc_number > 1)
      {
        CGAL_assertion(!vit->info().v->info().on_contour || cc_to_keep!=-1);
        //In case the vertex is not on the contour, the first component is the one kept
        if( !vit->info().v->info().on_contour ) cc_to_keep=0;
        
        #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
        //in case a non-manifold edge is found, all the components are modified because we have no guarantee that when
        //the second endpoint of the edge will be handled, the same component will be selected as the one "kept".
        if ( !vit->info().v->info().on_contour && !non_manifold_edges.empty() ) cc_to_keep=-1;
        for (int ccn=0;ccn<cc_number;++ccn){
          if (ccn==cc_to_keep) continue;

          std::size_t nb_mne=non_manifold_edges.size();
          //make each incident cell having non-manifold features (vertices and non-planar edges)
          for( typename std::vector<Cell_handle_3>::iterator itcell=incident_cells_per_cc[ccn].begin(),
                                                             itcell_end=incident_cells_per_cc[ccn].end();itcell!=itcell_end;++itcell)
          {
            (*itcell)->info().non_manifold_features=true;
            int index_of_vit=(*itcell)->index(vit);
            (*itcell)->info().set_nm_vertices(index_of_vit);
            //the loop is only done when a non-planar non-manifold edge has been found
            for (std::size_t kk=0;kk<nb_mne;++kk){
              Vertex_handle_3 vh=non_manifold_edges[kk];
              int index_of_vh=index_of_vertex_in_cell(*itcell,vh);
              if ( index_of_vh!=-1 )//the non-manifold edge is part of that cell
              {
                (*itcell)->info().set_nm_edge(*itcell,index_of_vit,index_of_vh);
              }
            }
          }
        }
        #else
        for (int ccn=0;ccn<cc_number;++ccn){
          if (ccn==cc_to_keep) continue;

          std::list<Facet_3>& boundary_facets=boundary_facets_per_cc[ccn];
          std::list<Facet_3>& boundary_planar_facets=boundary_planar_facets_per_cc[ccn];
          
          slice_writer_ptr->surface_point_push_back( get_point_3(vit->point()) );
          int vit_index=slice_writer_ptr->last_point_index();

          //handle non-planar facets in the component
          const cpp0x::array<int,3> init_array = {{-1,-1,-1}};
          for(typename std::list<Facet_3>::const_iterator bf_it=boundary_facets.begin(),
                                                          bf_it_end=boundary_facets.end();bf_it!=bf_it_end;++bf_it)
          {
            typename Facet_to_nm_indices::iterator fah_it=facets_already_handled.insert(std::make_pair(*bf_it,init_array)).first;
            int ii=0;
            for (;ii<3;++ii){
              if ( bf_it->first->vertex((bf_it->second+ii+1)%4)==Vertex_handle_3(vit) ){
                fah_it->second[ii]=vit_index;
                break;
              }
            }
            CGAL_assertion(ii!=3);
          }
         
          //handle planar facets in the component
          for(typename std::list<Facet_3>::iterator bf_it=boundary_planar_facets.begin(),
                                                    bf_it_end=boundary_planar_facets.end();bf_it!=bf_it_end;++bf_it)
          {
            //if (!vertex_on_contour) CGAL_assertion(!"We should get here"); this can also happen for vertex inside a contour!
            CGAL_precondition(!bf_it->first->info().volume);
            //only planar facet in the domain can be in such componant (otherwise, the out-of-domain facet can be incident to a contour)
            CGAL_precondition(bf_it->first->info().f0->info().in_domain()); 
            planar_facets_already_handled.insert(bf_it->first->info().f0);

            
            #warning this is dirty and should not be done like that: do something at the cell level
            Facet_3 opposite_facet(bf_it->first,bf_it->first->index(Vertex_handle_3(vit)));
            if ( opposite_facet.first->neighbor(opposite_facet.second)->info().volume ){
              typename Facet_to_nm_indices::iterator fah_it_2=facets_already_handled.insert(std::make_pair(opposite_facet,init_array)).first;
              for (int kkk=0;kkk<3;++kkk)
                if ( bf_it->first->vertex((opposite_facet.second+kkk+1)%4)->point().z()==vit->point().z() ) fah_it_2->second[kkk]=33;
            }
            
            typename Facet_to_nm_indices::iterator fah_it=planar_facets_nm_point.insert(std::make_pair(*bf_it,init_array)).first;
            int ii=0;
            for (;ii<3;++ii){
              if ( bf_it->first->vertex((bf_it->second+ii+1)%4)==Vertex_handle_3(vit) ){
                fah_it->second[ii]=vit_index;
                break;
              }
            }
            CGAL_assertion(ii!=3);
          }
        }        
        #endif
      }
    }
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    reinit_non_contour_vertices_poly_index(*top_ptr);
    #endif
    #endif //DO_NOT_HANDLE_NON_MANIFOLD_POINT
    
    #if !defined(CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS)
    if ( split_non_manifold_incontour_edges ){
      incontour_edge_ds.prepare();//init vertex indices to -1;
    
      //write planar facets having at least one non-manifold incontour edge
      for (typename DS_for_incontour_edge_handling::Planar_facets_to_split::iterator 
        it=incontour_edge_ds.planar_facets.begin(),it_end=incontour_edge_ds.planar_facets.end();it!=it_end;++it)
      {
        #ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
        write_planar_facets_with_non_manifold_incontour_edges(it->first,it->second,planar_facets_nm_point);
        #else
        write_planar_facets_with_non_manifold_incontour_edges(it->first,it->second);
        #endif
      }
    
      //write the special case when a contour is made of only 3 edges and that it is isolated (both incident
      //tetrahedra have been eliminated
      for (typename std::list<Cell_handle_3>::iterator cit=incontour_edge_ds.triangular_isolated_contours.begin();
                                                       cit!=incontour_edge_ds.triangular_isolated_contours.end(); ++cit)
      {
        int i=(*cit)->info().i0;
        int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
        if (i%2==0) std::swap(indices[0],indices[1]);
        for (int k=0;k<3;++k) indices[k]=get_vertex_index( (*cit)->vertex(indices[k]) );

        Point_2 pt_c2=centroid((*cit)->info().f0);
        double z = (*cit)->info().f0->vertex(0)->info().z;
        slice_writer_ptr->surface_point_push_back( get_point_3(pt_c2.x(),pt_c2.y(),z) );

        int cindex=slice_writer_ptr->last_point_index();
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[0],indices[1],cindex));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[1],indices[2],cindex));
        slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[2],indices[0],cindex));
      }
    }

    #ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
    //this should be done after calling write_planar_facets_with_non_manifold_incontour_edges
    // to delete planar facet with both non-manifold edges and vertices
    //write planar facets with at least one non-manifold vertex
    for (typename Facet_to_nm_indices::iterator it=planar_facets_nm_point.begin(),
                                       it_end=planar_facets_nm_point.end();it!=it_end;++it)
    {
      CGAL_assertion(it->first.first->info().i0==it->first.second);
      int i=it->first.first->info().i0;
      int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};

      for (int a=0;a<3;++a){
        if ( it->second[a]!=-1 ) indices[a]=it->second[a];
        else indices[a]=get_vertex_index(it->first.first->vertex(indices[a])->info().v);
      }

      if (i%2==0) std::swap(indices[1],indices[2]);
      slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[0],indices[1],indices[2]));                
      
      //TMP DEBUG
      #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
      additional_surface_points.push_back(get_additional_point(it->first.first->vertex((i+1)%4)->info().v));
      additional_surface_points.push_back(get_additional_point(it->first.first->vertex((i+2)%4)->info().v));
      additional_surface_points.push_back(get_additional_point(it->first.first->vertex((i+3)%4)->info().v));
      additional_surface_indices.push_back(cpp0x::make_tuple(additional_index+1,additional_index+2,additional_index+3));                
      additional_index+=3;
      #endif
    }
    #endif //!DO_NOT_HANDLE_NON_MANIFOLD_POINT
    #endif //CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    
    #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    //store indices of new points added on edges incident to a non-manifold vertex
    std::map<std::pair<Vertex_handle_3,Vertex_handle_3>,int> new_vertex_on_edge;
    //store indices of new points added in a facet, opposite to a non-manifold edge with manifold endpoints
    //two first vertices are sorted edge endpoint, third vertex is the last one in facet
    std::map< cpp0x::tuple<Vertex_handle_3,Vertex_handle_3,Vertex_handle_3>,int > new_opposite_nm_edge_vertex;
    
    //sets to detect non-manifold vertices and edges
    std::set<Vertex_handle_2>* current_bottom_incontour_nm_vertices=NULL;
    if (!last_run) current_bottom_incontour_nm_vertices=new std::set<Vertex_handle_2>();
    std::set< std::pair<Vertex_handle_2,Vertex_handle_2> >* current_bottom_incontour_nm_edges=NULL;
    if (!last_run) current_bottom_incontour_nm_edges=new std::set<std::pair<Vertex_handle_2,Vertex_handle_2> >();
    
    //mark planar non-manifold edges inside the contour as well as in contour vertices and non-manifold facets
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it)
      if ( !it->info().volume ){
        switch( it->info().type ){
          case CellInfo3::EDGE_EDGE:
          {
            //get indices of top and bottom indices
            int bi[2],ti[2];
            int *bi_ptr=bi, *ti_ptr=ti;
            for (int i=0;i<4;++i)
              if (it->vertex(i)->point().z()==top_z) 
                *ti_ptr++=i;
              else
                *bi_ptr++=i;
            CGAL_assertion(ti_ptr==ti+2 && bi_ptr==bi+2);
            CGAL_assertion(ti[0]+ti[1]+bi[0]+bi[1]==6);
            //handle the simplices in the top layer
            if ( it->info().f0->info().in_domain() && it->info().f0->neighbor(it->info().i0)->info().in_domain() ){  //the edge is inside the domain and not on the boundary
              //handle the edge
              if ( first_run || 
                previous_bottom_incontour_nm_edges->find( make_sorted_pair(it->vertex(ti[0])->info().v,it->vertex(ti[1])->info().v) )
                  !=previous_bottom_incontour_nm_edges->end() ||
                previous_layer_printed_planar_facets.find(it->info().f0)!=previous_layer_printed_planar_facets.end() ||
                previous_layer_printed_planar_facets.find(it->info().f0->neighbor(it->info().i0))!=previous_layer_printed_planar_facets.end()
              ) {  it->info().set_nm_edge(it,ti[0],ti[1]); }
              //handle the vertices
                if (!it->vertex(ti[0])->info().v->info().on_contour && (first_run || previous_bottom_incontour_nm_vertices->find(it->vertex(ti[0])->info().v)!=previous_bottom_incontour_nm_vertices->end()) )
                  it->info().set_nm_vertices(ti[0]);
                if (!it->vertex(ti[1])->info().v->info().on_contour && (first_run || previous_bottom_incontour_nm_vertices->find(it->vertex(ti[1])->info().v)!=previous_bottom_incontour_nm_vertices->end()) )
                  it->info().set_nm_vertices(ti[1]);
            }


            //handle the simplices in bottom layer
            if ( it->info().f1->info().in_domain() && it->info().f1->neighbor(it->info().i1)->info().in_domain() )
            {
              //the edge is inside the domain and not on the boundary
              if (last_run){
                it->info().set_nm_edge(it,bi[0],bi[1]);
                if ( !it->vertex(bi[0])->info().v->info().on_contour ) it->info().set_nm_vertices(bi[0]);
                if ( !it->vertex(bi[1])->info().v->info().on_contour ) it->info().set_nm_vertices(bi[1]);
              }
              else{
                current_bottom_incontour_nm_edges->insert( make_sorted_pair(it->vertex(bi[0])->info().v,it->vertex(bi[1])->info().v) );
                if ( !it->vertex(bi[0])->info().v->info().on_contour ) current_bottom_incontour_nm_vertices->insert(it->vertex(bi[0])->info().v);
                if ( !it->vertex(bi[1])->info().v->info().on_contour ) current_bottom_incontour_nm_vertices->insert(it->vertex(bi[1])->info().v);
              }
            }
            break;
          }
          case CellInfo3::TOP_TRIANGLE:
          {
            if ( it->info().f0->info().in_domain() ){
              //handle the facet
              if( first_run || previous_layer_printed_planar_facets.find(it->info().f0)!=previous_layer_printed_planar_facets.end() )
                it->info().set_nm_facet(it);
              else{
                int indices[3]={ (it->info().i0+1)%4, (it->info().i0+2)%4, (it->info().i0+3)%4 };
                for (int i=0;i<3;++i){
                  //handle edge
                  if (first_run || previous_bottom_incontour_nm_edges->find( make_sorted_pair(it->vertex(indices[i])->info().v,it->vertex(indices[(i+1)%3])->info().v) )
                        !=previous_bottom_incontour_nm_edges->end() )
                  {
                    it->info().set_nm_edge(it,indices[i],indices[(i+1)%3]);
                  }
                  //handle vertex
                  if ( !it->vertex(indices[i])->info().v->info().on_contour &&
                       ( first_run || previous_bottom_incontour_nm_vertices->find(it->vertex(indices[i])->info().v)!=previous_bottom_incontour_nm_vertices->end()) )
                  {
                    it->info().set_nm_vertices(indices[i]);
                  }
                }
              }
            }
            //handle the vertex isolated in the bottom layer
            if (!it->vertex(it->info().i0)->info().v->info().on_contour){
              if (last_run)
                it->info().set_nm_vertices(it->info().i0);
              else
                current_bottom_incontour_nm_vertices->insert(it->vertex(it->info().i0)->info().v);
            }
            break;
          }
          case CellInfo3::BOTTOM_TRIANGLE:
          {
            //handle vertex isolated in the top layer
            if (!it->vertex(it->info().i0)->info().v->info().on_contour){
              if (first_run || previous_bottom_incontour_nm_vertices->find(it->vertex(it->info().i0)->info().v)!=previous_bottom_incontour_nm_vertices->end() )
                it->info().set_nm_vertices(it->info().i0);
            }
            //handle facet
            if( it->info().f0->info().in_domain() ){
              if (last_run)
                it->info().set_nm_facet(it);
              else{
                const int indices[3]={ (it->info().i0+1)%4, (it->info().i0+2)%4, (it->info().i0+3)%4 };
                current_bottom_incontour_nm_vertices->insert(it->vertex(indices[0])->info().v);
                current_bottom_incontour_nm_vertices->insert(it->vertex(indices[1])->info().v);
                current_bottom_incontour_nm_vertices->insert(it->vertex(indices[2])->info().v);
              }
            }
          }
        }
      }

    if (!first_run){
      delete previous_bottom_incontour_nm_edges;
      delete previous_bottom_incontour_nm_vertices;
    }
    std::swap(previous_bottom_incontour_nm_edges,current_bottom_incontour_nm_edges);
    std::swap(previous_bottom_incontour_nm_vertices,current_bottom_incontour_nm_vertices);
    previous_layer_printed_planar_facets.clear();
    #endif
    
    //write surface facets for the current slice
    for(Cell_iterator_3 it = delaunay_3.finite_cells_begin(); it != delaunay_3.finite_cells_end(); ++it){
      #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
      if ( it->info().non_manifold_features ){
        display_surface_facet_in_cells_with_non_manifold_features(it,previous_layer_printed_planar_facets,new_vertex_on_edge,new_opposite_nm_edge_vertex);
        continue;
      }
      #endif
      
      //Handle regular facets on the boundary of in volume tetrahedron
      if(it->info().volume){
        for(int i=0; i < 4; ++i){
          if( (it->info().type == CellInfo3::EDGE_EDGE) || (it->info().i0 != i)){
            
            #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
            if(! it->neighbor(i)->info().volume && !it->neighbor(i)->info().non_manifold_features)
            #else
            if(! it->neighbor(i)->info().volume)
            #endif
            {
              Facet_3 mirror_facet=delaunay_3.mirror_facet(Facet_3(it,i));
              int indices[3]={(mirror_facet.second+1)%4,(mirror_facet.second+2)%4,(mirror_facet.second+3)%4};

              //catch by non-manifold vertex
              #if !defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT) && !defined(CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS)
              typename Facet_to_nm_indices::iterator fah_it= facets_already_handled.find(mirror_facet); 
              if ( fah_it!=facets_already_handled.end() ){
                //TMP DEBUG
                #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
                additional_surface_points.push_back(get_additional_point(it->vertex((i+1)%4)->info().v));
                additional_surface_points.push_back(get_additional_point(it->vertex((i+2)%4)->info().v));
                additional_surface_points.push_back(get_additional_point(it->vertex((i+3)%4)->info().v));
                additional_surface_indices.push_back(cpp0x::make_tuple(additional_index+1,additional_index+2,additional_index+3));                
                additional_index+=3;                
                #endif
                for (int a=0;a<3;++a){
                  if ( fah_it->second[a]!=-1 ) indices[a]=fah_it->second[a];
                  else indices[a]=get_vertex_index(mirror_facet.first->vertex(indices[a])->info().v);
                }                
              }
              else
              #endif
              {
                indices[0]=get_vertex_index( mirror_facet.first->vertex(indices[0]) );
                indices[1]=get_vertex_index( mirror_facet.first->vertex(indices[1]) );
                indices[2]=get_vertex_index( mirror_facet.first->vertex(indices[2]) );
              }
              if (mirror_facet.second%2==0) std::swap(indices[1],indices[2]);
             
              #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
              //split incontour non-manifold edges
              if (split_non_manifold_incontour_edges){
                typename DS_for_incontour_edge_handling::Non_planar_facets_to_split::iterator it_np=
                  incontour_edge_ds.non_planar_facets.find(mirror_facet);
                if ( it_np!=incontour_edge_ds.non_planar_facets.end() ){
                  int k=mirror_facet.second;
                  int k_indices[3]={(k+1)%4,(k+2)%4,(k+3)%4};
                  if (k%2==0) std::swap(k_indices[1],k_indices[2]);
                  //set k so that k_indices[k] is the index of the vertex opposite to the edge to be split
                  for (k=0;k<3;++k)
                    if ( k_indices[k]==it_np->second.opp_e ) break;
                  CGAL_assertion(k!=3);
                  Vertex_handle_3 ve1=mirror_facet.first->vertex(k_indices[(k+1)%3]);
                  Vertex_handle_3 ve2=mirror_facet.first->vertex(k_indices[(k+2)%3]);
                  //get index of the midpoint
                  int mp_index=get_vertex_index_on_edge(ve1,ve2,it_np->second.mp_index);
                  //split one facet in two facets
                  slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[(k+1)%3],mp_index,indices[k]));
                  slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(mp_index,indices[(k+2)%3],indices[k]));
                  continue;
                }
              }
              #endif
              slice_writer_ptr->surface_indices_push_back(cpp0x::make_tuple(indices[0],indices[1],indices[2]));
            }
          }
        }
      }
      else{
        //handle planar facets
        if( it->info().type==CellInfo3::EDGE_EDGE  || !it->info().f0->info().in_domain() ) continue;
        int i=it->info().i0;
        int indices[3]={(i+1)%4,(i+2)%4,(i+3)%4};
        if (i%2==0)
          std::swap(indices[1],indices[2]);
        
        if (it->info().type==CellInfo3::TOP_TRIANGLE){
          #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
          if (planar_facets_already_handled.find(it->info().f0) == planar_facets_already_handled.end() )
          #endif
              slice_writer_ptr->surface_indices_push_back(
                cpp0x::make_tuple(
                  get_vertex_index(it->vertex(indices[0])),
                  get_vertex_index(it->vertex(indices[1])),
                  get_vertex_index(it->vertex(indices[2]))
                )
              );
        }
        else{
          CGAL_assertion(it->info().type==CellInfo3::BOTTOM_TRIANGLE);
          
          #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
          if ( planar_facets_already_handled.find(it->info().f0)==planar_facets_already_handled.end() )
          #endif
            slice_writer_ptr->surface_indices_push_back(
              cpp0x::make_tuple(
                  get_vertex_index(it->vertex(indices[0])),
                  get_vertex_index(it->vertex(indices[1])),
                  get_vertex_index(it->vertex(indices[2]))
              )
            );

          #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
          if ( !last_run )
            //collect facets that have been printed to detect non-manifold facet at the next step
            previous_layer_printed_planar_facets.insert(it->info().f0);
          #else
          if ( split_non_manifold_incontour_edges && !last_run){
            //collect facets that have been printed to detect non-manifold facet at the next step
            previous_layer_printed_planar_facets.insert(it->info().f0);
            
            //add incontour edges that are in planar facets in the potental list of non-manifold edges
            Face_handle_2 fh=it->info().f0;
            for (int k=0;k<3;++k){
              if (fh->vertex((k+1)%3)->info().on_contour && 
                  fh->vertex((k+2)%3)->info().on_contour && 
                  !fh->is_constrained(k)
                 )
              {
                Vertex_handle_3 v1,v2;
                int j=0;
                for (;j<3;++j)
                  if( fh->vertex((k+1)%3)==it->vertex(indices[j])->info().v ){
                    v1=it->vertex(indices[j]);
                    break;
                  }
                CGAL_assertion(j!=3);
                for (j=0;j<3;++j)
                  if( fh->vertex((k+2)%3)==it->vertex(indices[j])->info().v ){
                    v2=it->vertex(indices[j]);
                    break;
                  }
                CGAL_assertion(j!=3);
                std::pair<Vertex_handle_2,Vertex_handle_2> vertex_pair=make_sorted_pair(v1->info().v,v2->info().v);               
                previous_layer_incontour_edges->insert(vertex_pair);
              }
            }
          }
          #endif
        }
      }
    }

    if (!last_run) slice_writer_ptr->finalize_layer();
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    planar_facets_already_handled.clear();
    #endif
  }

  #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
  template <class Points,class Indices>
  void write_final_model(std::ostream& output,const Points& points,const Indices& indices){ //OUTPUT
    output << "OFF " << points.size() << " " << indices.size() << " 0\n";
    std::copy(points.begin(),points.end(),std::ostream_iterator<Point_3>(output,"\n"));
    for( typename Indices::const_iterator it=indices.begin();it!=indices.end();++it)
      output << "3 " << cpp0x::get<0>(*it) << " " << cpp0x::get<1>(*it) << " " << cpp0x::get<2>(*it) << "\n";
  }

  void write_additional_facets(std::ostream& output){ //OUTPUT
    write_final_model(output,additional_surface_points,additional_surface_indices);
  }
  #endif
  
  #ifndef DO_NOT_HANDLE_NON_MANIFOLD_POINT
  //vertices that are inside the contour are duplicated in each reconstructed slice
  void reinit_non_contour_vertices_poly_index(CDT2& cdt){
    for (typename CDT2::Finite_vertices_iterator vit=cdt.finite_vertices_begin(),end=cdt.finite_vertices_end();vit!=end;++vit){
      if ( !vit->info().on_contour )
        vit->info().poly_index=-1;
    }
  }
  #endif
  
  void run(bool last_run){
    if (!last_run){
      CGAL::make_conforming_Gabriel_2(*next_ptr);
      mark_domains(*next_ptr);
    }
  //insert exterior Voronoi vertices of top and next that are inside the contour of bottom
    std::list<Point_2> add_to_bottom;
    vertices_to_add(*top_ptr, *bottom_ptr, add_to_bottom);
    if (!last_run){
#ifndef DO_NOT_INTERSECT_CONTOURS_WITH_MEDIAL_AXIS
      mark_domains(*bottom_ptr);
#endif
      vertices_to_add(*next_ptr, *bottom_ptr, add_to_bottom);
    }
  //insert selected vertices into constrainted triangulations
    add_vertices_inside(*bottom_ptr, add_to_bottom);
  //make both constrained triangulations conforming Gabriel after point insertion
    CGAL::make_conforming_Gabriel_2(*bottom_ptr);
  //update triangles inside the contours
    mark_domains(*bottom_ptr);
  //after refining a CDT we have to set the z values for the new vertices 
    update_z(*bottom_ptr);

    //create a 3D Delaunay triangulation from vertices in top and bottom and
  //classify each tetrahedron:either having one face in P1, or one face in P2 or one edge in P1 and one edge in P2
    delaunay_3.clear();
    create_tetrahedrization(*top_ptr, *bottom_ptr);
  //remove tetrahedra that should not be in the final volumic slice
    remove_cells();
//    check_removed_cells();//DEBUG
    #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
  //DEBUG OUTPUT: assign an index to each finite vertex  
    index(*bottom_ptr);
  //DEBUG OUTPUT: write CDTs into off format
    std::ofstream output_1("data/bottom.off");
    std::ofstream output_2("data/top.off");
    write_layer(output_1,*bottom_ptr);
    write_layer(output_2, *top_ptr);
    output_1.close();
    output_2.close();
    output_1.open("data/bottom-out.off");
    output_2.open("data/top-out.off");
    write_layer(output_1,*bottom_ptr,false);
    write_layer(output_2, *top_ptr,false);
//DEBUG code to display projected circumcenters
//    output_1.close();
//    output_2.close();
//    output_1.open("data/bottom-out.cgal");
//    output_2.open("data/top-out.cgal");   
//    write_projected_circumcenters(output_1,*bottom_ptr,top_ptr->finite_vertices_begin()->info().z);
//    write_projected_circumcenters(output_2,*top_ptr,bottom_ptr->finite_vertices_begin()->info().z);
    #endif
    //write_voronoi(std::ofstream("bottomVoronoi.cgal"),bottom);
    //write_voronoi(std::ofstream("topVoronoi.cgal"), top);
    
//  std::ofstream output("last-graph.off");//DEBUG
//  write_cells(output);//DEBUG
//  write_deleted_tetra();//DEBUG
//    bool first_run=previous_layer_incontour_edges==NULL;
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    if( split_non_manifold_incontour_edges ) 
      split_non_manifold_edges_between_contour_points(top_ptr->finite_vertices_begin()->info().z,last_run);
    #endif
    save_surface_layer(last_run);
    #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
    std::cout << "X" << std::flush;
    #endif
  }

  template <class PolygonContainer>
  void read(PolygonContainer& polygon_reader, CDT2*& cdt_ptr)
  {
    Point_2 p2, q2;
    Vertex_handle_2 vp, vq;
    cdt_ptr=new CDT2();
    double z;
    do{
      bool contour_made_of_one_point=true;
      int n=polygon_reader.number_of_points();
      p2 = polygon_reader.get_point(z);
      vp = cdt_ptr->insert(p2);
      vp->info().z = z;
      for(int i=1; i < n; i++){
        q2 = polygon_reader.get_point();

        if (contour_made_of_one_point && p2==q2) continue;
        else contour_made_of_one_point=false;
        
        vq = cdt_ptr->insert(q2, vp->face());
        vq->info().z = z;
        if (vp!=vq)
          cdt_ptr->insert_constraint(vp, vq);

        vp = vq;
        p2 = q2;
      }
      if (contour_made_of_one_point) vp->info().on_contour=false;
      if (polygon_reader.has_another_component()) polygon_reader.next_polygon();
      else break;
    }while(true);
  }
  
public:

  Reconstruction_from_parallel_slices_3():
    top_ptr(NULL),bottom_ptr(NULL),next_ptr(NULL)
  #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
  ,additional_index(-1)
  #endif
  #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  ,previous_bottom_incontour_nm_vertices(NULL)
  ,previous_bottom_incontour_nm_edges(NULL)
  #endif

  {}

  template <class PolygonContainer>
  void run(PolygonContainer& polygon_reader,Slice_writer& slice_writer, unsigned int constant_coordinate){
    slice_writer_ptr=&slice_writer;
    x_index=2-constant_coordinate; //index of x
    CGAL_assertion(!polygon_reader.empty() || !"No contour has been found");
    read(polygon_reader, bottom_ptr);
    CGAL_assertion(polygon_reader.has_next_planar_contour() || !"At least two contours should be available");
    polygon_reader.next_polygon();
    read(polygon_reader, next_ptr);
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  //init container to check no two opposite planar facets  are in the final volume
    previous_layer_incontour_edges=NULL;
    #endif
  //initialization
    CGAL::make_conforming_Gabriel_2(*bottom_ptr);
    CGAL::make_conforming_Gabriel_2(*next_ptr);
  //mark triangles that are inside and outside domains
    mark_domains(*bottom_ptr);
    mark_domains(*next_ptr);
  //insert exterior Voronoi vertices of next that are inside the contour of bottom
    std::list<Point_2> add_to_bottom;
    vertices_to_add(*next_ptr, *bottom_ptr, add_to_bottom);
  //insert selected vertices into constrainted triangulation
    add_vertices_inside(*bottom_ptr, add_to_bottom);
  //make conforming Gabriel after point insertion
    CGAL::make_conforming_Gabriel_2(*bottom_ptr);
  //update triangles inside the contours
    mark_domains(*bottom_ptr);
    update_z(*bottom_ptr);
    index(*bottom_ptr);

    bool is_increasing_z= bottom_ptr->vertices_begin()->info().z < 
                          next_ptr->vertices_begin()->info().z;
    //write the planar facets from the top layer
    write_first_planar_facets(*bottom_ptr,is_increasing_z);
    do{
      if (top_ptr) delete top_ptr;
      top_ptr=bottom_ptr;
      bottom_ptr=next_ptr;
      next_ptr=NULL;

      if (polygon_reader.has_next_planar_contour()){
        polygon_reader.next_polygon();
        read(polygon_reader, next_ptr);
      }
      run(next_ptr==NULL);
    }while (next_ptr!=NULL);
    #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
    std::cout << std::endl;
    #endif
    add_last_bottom_layer_facets_to_surface(*bottom_ptr,is_increasing_z);
    delete bottom_ptr;
    delete top_ptr;
    #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
    if( split_non_manifold_incontour_edges && previous_layer_incontour_edges!=NULL)
      delete previous_layer_incontour_edges;
    #endif
    slice_writer_ptr->finalize();
    #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
    std::ofstream other("additional.off");
    write_additional_facets(other);
    #endif
  }
  
private:
  //members
  DT3 delaunay_3;
  CDT2* top_ptr,*bottom_ptr,*next_ptr;

  #ifdef CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3_DEBUG
  int additional_index;
  std::list<cpp0x::tuple<int,int,int> > additional_surface_indices;
  std::list<Point_3> additional_surface_points;
  #endif

  Slice_writer* slice_writer_ptr;

  #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  //set of planar facets with at least one non-manifold vertex or non-manifold edge
  std::set<Face_handle_2> planar_facets_already_handled;
  #endif

  //set of planar facets that have been displayed in the previous bottom layer (which is the top layer in the current run)
  std::set<Face_handle_2> previous_layer_printed_planar_facets;

  #ifndef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  Unique_incontour_edges* previous_layer_incontour_edges;
  #endif

  unsigned int x_index;
  
  #ifdef CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS
  //use to track in-contour non-manifold points
  std::set<Vertex_handle_2>* previous_bottom_incontour_nm_vertices;
  //used to track planar non-manifold edges in T22 that are inside to the contour
  std::set< std::pair<Vertex_handle_2,Vertex_handle_2> >* previous_bottom_incontour_nm_edges;
  #endif
};

#if !defined(DO_NOT_HANDLE_NON_MANIFOLD_POINT) && defined(CGAL_ADD_VOLUME_TO_REMOVE_NON_MANIFOLDNESS)
template <class Slice_writer>
const int Reconstruction_from_parallel_slices_3<Slice_writer>::CellInfo3::nm_edge_index[4][4] = {
  {-1,0,1,2}, //0
  {0,-1,3,4}, //1
  {1,3,-1,5}, //2
  {2,4,5,-1} }; //3
#endif

} //namespace CGAL

#endif //CGAL_RECONSTRUCTION_FROM_PARALLEL_SLICES_3
