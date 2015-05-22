// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_PARAMETERIZATION_POLYHEDRON_ADAPTOR3_H
#define CGAL_PARAMETERIZATION_POLYHEDRON_ADAPTOR3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Convertible_iterator_project.h>
#include <CGAL/Convertible_circulator_project.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/foreach.hpp>
#include <list>
#include <set>

/// \file Parameterization_polyhedron_adaptor_3.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationMesh
///
/// `Parameterization_polyhedron_adaptor_3` is an adaptor class to access to a Polyhedron
/// 3D mesh using the `ParameterizationPatchableMesh_3` interface.
/// Among other things, this concept defines the accessor to the (u,v) values
/// computed by parameterizations methods.
///
/// Note that these interfaces are decorators that add "on the fly"
/// the necessary fields to unmodified CGAL data structures (using STL maps).
/// For performance reasons, it is recommended to use CGAL data structures
/// enriched with the proper fields.
///
/// A `ParameterizationMesh_3` surface consists of vertices,
/// facets and an incidence relation on them.
/// No notion of edge is requested.
///
/// `ParameterizationMesh_3` meshes can have any genus, arity or number of components.
///
/// It can have have any number of borders. Its "main border"
/// will be the mesh's longest border (if there is at least one border).
///
/// It has also the ability to support patches and virtual seams.
/// <i>Patches</i> are a subset of a 3D mesh. "Virtual seams" are the ability
/// to behave exactly as if the surface was cut following a certain path.
///
/// \cgalModels `ParameterizationPatchableMesh_3`
///

  template<class Polyhedron_3_, class VertexIndexMap, class HalfedgeIndexMap, class HalfedgeUvMap>
class Parameterization_polyhedron_adaptor_3
{
// Forward references
public:
    class                                   Halfedge_info;
    class                                   Vertex_info;

// Private types
private:

  typedef typename Polyhedron_3_ TriangleMesh;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator face_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_iterator halfedge_iterator;
 
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Halfedge_around_target_circulator<TriangleMesh> halfedge_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

    /// Additional info attached to halfedges
    typedef typename std::vector<Halfedge_info>         Halfedge_info_map;
    /// Additional info attached to vertices
    typedef typename std::vector<Vertex_info>         Vertex_info_map;

  typedef Parameterization_polyhedron_adaptor_3<Polyhedron_3_, VertexIndexMap,HalfedgeIndexMap, HalfedgeUvMap> Adaptor;

// Public types
public:

    //*******************************************************************
    // @name Interface specific to Parameterization_polyhedron_adaptor_3
    //*******************************************************************
    //@{

    /// Export template parameter.
    typedef Polyhedron_3_                  Polyhedron;

    /// \cond SKIP_IN_MANUAL

    /// Additional info attached to halfedges.
    class Halfedge_info
    {
    private:
        int m_tag;                  ///< general purpose tag
        bool m_is_parameterized;    ///< is parameterized?
        int m_seaming;              ///< seaming status
        int m_index;                ///< unique index
        int m_vindex;               /// AF
        int m_on_border;            /// AF

    public:
        /// Default constructor.
        Halfedge_info()
        {
            m_tag = -1;             // uninitialized
            m_index = -1;           // uninitialized
            m_seaming = 0;         // uninitialized
            m_is_parameterized = false;
            m_vindex = -1;
            m_on_border = 0;
        }

        // Default destructor, copy constructor and operator =() are fine

        /// Access to general purpose tag.
        int tag() const { return m_tag; }
        void tag(int tag) { m_tag = tag; }

        /// Access to 'seaming status' field.
        int seaming() const { return m_seaming; }
        void seaming(int seaming) { m_seaming = seaming; }


        /// Access to "parameterized?" field.
        bool is_parameterized() const { return m_is_parameterized; }
        void is_parameterized(bool is)  { m_is_parameterized = is; }

        /// Access to 'index' field.
        int index() const { return m_index; }
        void index(int i) { m_index = i; } 

        int vindex() const { return m_vindex; }
        void vindex(int i) { m_vindex = i; }

        int on_border() const { return m_on_border; }
        void on_border(int i) { m_on_border = i; }
    };

    /// Additional info attached to vertices.
    class Vertex_info
    {
    private:
        int m_tag;                  ///< general purpose tag
        int m_seaming;              ///< seaming status
        int m_index;                ///< unique index
        int m_on_border;           ///<  0=inner, 2= main border

    public:
        /// Default constructor.
        Vertex_info()
        {
            m_index = -1;           // uninitialized
            m_tag = -1;             // uninitialized
            m_seaming = 0;         // uninitialized
            m_on_border = 0;
        }

        // Default destructor, copy constructor and operator =() are fine.

        /// Access to 'index' field.
        int index() const { return m_index; }
        void index(int i) { m_index = i; }

        /// Access to 'tag' field.
        int tag() const { return m_tag; }
        void tag(int tag) { m_tag = tag; }

        /// Access to 'seaming status' field.
        int seaming() const { return m_seaming; }
        void seaming(int seaming) { m_seaming = seaming; }

        int on_border() const { return m_on_border; }
        void on_border(int i) { m_on_border = i; }
    };

    /// \endcond
    
    //@} // end of Interface specific to Parameterization_polyhedron_adaptor_3

    #ifndef DOXYGEN_RUNNING
    //*******************************************************************
    /// @name ParameterizationMesh_3 INTERFACE
    //*******************************************************************
    //@{

    /// Number type to represent coordinates.
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type PPmap;
  typedef typename boost::property_traits<PPmap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;

    typedef typename Kernel::FT NT;

    /// 2D point that represents (u,v) coordinates computed
    /// by parameterization methods. Must provide X() and Y() methods.
    typedef typename Kernel::Point_2  Point_2;
    /// 3D point that represents vertices coordinates. Must provide X() and Y() methods.
    typedef typename Kernel::Point_3   Point_3;
    /// 2D vector. Must provide X() and Y() methods.
    typedef typename Kernel::Vector_2   Vector_2;
    /// 3D vector. Must provide X() and Y() methods.
    typedef typename Kernel::Vector_3   Vector_3;

    /// Iterator over vertices of the mesh "main border".
    /// Model of the ForwardIterator concept.
    typedef typename std::list<halfedge_descriptor>::const_iterator Border_halfedge_iterator;

   
    //@} // end of ParameterizationMesh_3 INTERFACE
  #endif //DOXYGEN_RUNNING

// Public operations
public:

    //*******************************************************************
    // @name Interface specific to Parameterization_polyhedron_adaptor_3
    //*******************************************************************
    //@{

    /// Create an adaptator for an existing Polyhedron_3 mesh.
    /// The input mesh can be of any genus.
    /// It can have have any number of borders.
    /// border will be the main border
  Parameterization_polyhedron_adaptor_3(Polyhedron& mesh,
                                        typename boost::graph_traits<Polyhedron>::halfedge_descriptor border,
                                        VertexIndexMap vim,
                                        HalfedgeIndexMap him,
                                        HalfedgeUvMap huvm)
        // Store reference to adapted mesh
    : m_polyhedron(mesh), m_vim(vim), m_him(him), m_huvm(huvm)
    {
        typedef typename Halfedge_info_map::value_type Halfedge_info_pair;
        typedef typename Vertex_info_map::value_type Vertex_info_pair;

        m_halfedge_info.resize(num_halfedges(mesh));
        m_vertex_info.resize(num_vertices(mesh));

        // Allocate extra info for each vertex
        m_num_vertices = 0;
        BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
          set_vertex_index(vd, m_num_vertices++);
        }

        BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
          info(hd)->vindex(info(target(hd,mesh))->index());
        }
        // if there is no seam each halfedge stores the index of the target
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(border,mesh)){
          m_main_border.push_back(hd);
          info(hd)->on_border(2);
          info(target(hd,mesh))->on_border(2);
        }
    }



  /// Create an adaptator for an existing Polyhedron_3 mesh.
  /// The input mesh can be of any genus.
  /// It can have have any number of borders.
  /// the seam will be the main border
  
  ///
  Parameterization_polyhedron_adaptor_3(Polyhedron& mesh,
                                        std::vector<edge_descriptor>& seam,
                                        VertexIndexMap vim,
                                        HalfedgeIndexMap him,
                                        HalfedgeUvMap huvm)
        // Store reference to adapted mesh
    : m_polyhedron(mesh), m_vim(vim), m_him(him), m_huvm(huvm)
  {
    typedef typename Halfedge_info_map::value_type Halfedge_info_pair;
    typedef typename Vertex_info_map::value_type Vertex_info_pair;

    m_halfedge_info.resize(num_halfedges(mesh));
    m_vertex_info.resize(num_vertices(mesh));

    // each vertex shall know on how many seaming edges it is
    BOOST_FOREACH(edge_descriptor e, seam){
      info(halfedge(e,mesh))->seaming(1);
      info(opposite(halfedge(e,mesh),mesh))->seaming(1);
      vertex_descriptor vd = target(e,mesh);
      info(vd)->seaming(info(vd)->seaming()+1);
      vd = source(e,mesh);
      info(vd)->seaming(info(vd)->seaming()+1);
    }

    // we only want to use the index stored in info if it is not inside a seam
    m_num_vertices = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      if(info(vd)->seaming() == 0){
        set_vertex_index(vd, m_num_vertices++);
      } else {
        assert(get_vertex_index(vd) == -1);
      } 
    }
    BOOST_FOREACH(edge_descriptor e, seam){
      halfedge_descriptor hd = halfedge(e,mesh);
      info(target(hd,mesh))->on_border(2);
      set_vertex_index(hd,m_num_vertices);
      while(info(next(hd,mesh))->seaming()==0){
        hd = opposite(next(hd,mesh),mesh);
        set_vertex_index(hd,m_num_vertices);
      }
      ++m_num_vertices;
      hd = opposite(halfedge(e,mesh),mesh);
      info(target(hd,mesh))->on_border(2);
      set_vertex_index(hd,m_num_vertices);
      while(info(next(hd,mesh))->seaming()==0){
        hd = opposite(next(hd,mesh),mesh);
        set_vertex_index(hd,m_num_vertices);
      
      }
      ++m_num_vertices;
    }

    // all halfedges where target is not on a seam will get it from the target
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
      if(info(target(hd,mesh))->seaming()==0){
        assert(info(hd)->seaming()==0);
        assert(info(target(hd,mesh))->index() != -1);
        set_vertex_index(hd, info(target(hd,mesh))->index());
      }
      assert(get_vertex_index(hd)!= -1);
    }
    

    halfedge_descriptor hd = halfedge(seam[0],mesh), done(hd);
    do {
      m_main_border.push_back(hd);
      while(info(next(hd,mesh))->seaming()!=1){
        hd = opposite(next(hd,mesh),mesh);
      }
      hd =next(hd,mesh);
    } while(hd != done);

    std::cerr << "main border initialized\n";
#if 0
    // to get started deal with a polyline border 
    for(int i=0; i < seam.size(); i++){
      m_main_border.push_back(halfedge(seam[i],mesh));
    }
    for(int i=seam.size(); i>0 ; i--){
      m_main_border.push_back(opposite(halfedge(seam[i-1],mesh),mesh));
    }
#endif
  }



    // Default destructor, copy constructor and operator =() are fine

    /// Get the adapted mesh.
    Polyhedron&       get_adapted_mesh()       { return m_polyhedron; }
    const Polyhedron& get_adapted_mesh() const { return m_polyhedron; }

    /// Get halfedge from source and target vertices.
    /// Will assert if such a halfedge doesn't exist.
    halfedge_descriptor get_halfedge(
        vertex_descriptor source, vertex_descriptor target) const
    {
      std::pair<halfedge_descriptor,bool> res = halfedge(source,target, m_polyhedron);
      assert(res.second);
      return res.first;
    }

    /// Access to additional info attached to halfedges.
    const Halfedge_info* info(halfedge_descriptor halfedge) const
    {
      return & m_halfedge_info[get(m_him,halfedge)];
    }

    Halfedge_info* info(halfedge_descriptor halfedge)
    {
      return & m_halfedge_info[get(m_him,halfedge)];
    }

    /// Access to additional info attached to vertices.
    const Vertex_info* info(vertex_descriptor vertex) const
    {
      return & m_vertex_info[get(m_vim,vertex)];
    }

    Vertex_info* info(vertex_descriptor vertex)
    {
      return & m_vertex_info[get(m_vim,vertex)];
    }

    //@} // end of Interface specific to Parameterization_polyhedron_adaptor_3

    #ifndef DOXYGEN_RUNNING
    //*******************************************************************
    /// @name ParameterizationMesh_3 INTERFACE
    //*******************************************************************
    //@{

    // MESH INTERFACE

    /// Indicate if the mesh matches the ParameterizationMesh_3 concept.
    bool is_valid() const {
        return m_polyhedron.is_valid();
    }

    /// Count the number of vertices of the mesh.
    int  count_mesh_vertices() const {
      return m_num_vertices;
    }


    const std::list<halfedge_descriptor>& main_border() const { return m_main_border; }

    /// Return the border containing seed_vertex.
    /// Return an empty list if not found.
    std::list<vertex_descriptor> get_border(vertex_descriptor seed_vertex)
    {
        std::list<vertex_descriptor> border;    // returned list

        // if isolated vertex
        if(halfedge(seed_vertex, m_polyhedron) ==  boost::graph_traits<Polyhedron>::null_halfedge()){
          border.push_back(seed_vertex);
          return border;
        }

        // Get seed_vertex' border halfedge
        halfedge_descriptor  seed_halfedge = boost::graph_traits<Polyhedron>::null_halfedge();
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(seed_vertex,m_polyhedron)) {
          if(is_border(hd, m_polyhedron)) {
                seed_halfedge = hd;
                break;
            }
        }

        // if inner vertex
        if (seed_halfedge == boost::graph_traits<Polyhedron>::null_halfedge())
            return border;                  // return empty list

        // Add seed vertex
        border.push_back(seed_vertex);

        // fill border
        int size = 1;
        halfedge_descriptor current_halfedge = seed_halfedge;
        do
        {
            // Stop if end of loop
          halfedge_descriptor next_halfedge = next(current_halfedge, m_polyhedron);
          vertex_descriptor next_vertex = target(next_halfedge, m_polyhedron);
            if(next_vertex == seed_vertex)
                break;

            // Add vertex
            border.push_back(next_vertex);

            current_halfedge = next_halfedge;
            size++;
        }
        while(1);

        return border;
    }

    /// Get iterator over first facet of mesh.
    face_iterator  mesh_facets_begin()const {
      return faces(m_polyhedron).first;
    }

    face_iterator  mesh_facets_end() const {
      return faces(m_polyhedron).second;
    }

    /// Count the number of facets of the mesh.
    int  count_mesh_facets() const {
        face_iterator b,e;
        boost::tie(b,e) = faces(m_polyhedron);

        return std::distance(b,e);
    }

    /// Count the number of halfedges of the mesh.
    int  count_mesh_halfedges() const {
      halfedge_iterator b,e;
      boost::tie(b,e) = halfedges(m_polyhedron);
      return std::distance(b,e);
    }

    // VERTEX INTERFACE

    /// Get/set the 2D position (u/v pair) of a vertex. Default value is undefined.
    /// (stored in halfedges sharing the same vertex).
    Point_2  get_vertex_uv(vertex_descriptor vertex) const {
      return get_corners_uv(vertex, vertex_descriptor(), vertex_descriptor());
    }
    void  set_vertex_uv(vertex_descriptor vertex, const Point_2& uv) {
      set_corners_uv(vertex, vertex_descriptor(), vertex_descriptor(), uv);
    }

    void  set_vertex_uv(halfedge_descriptor halfedge, const Point_2& uv) {
      if(is_border(halfedge,m_polyhedron)){
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(halfedge,m_polyhedron)){
          put(m_huvm, hd, uv);
        }
      }else {
        put(m_huvm, halfedge, uv);
        while(info(next(halfedge,m_polyhedron))->seaming()==0){
          halfedge = opposite(next(halfedge,m_polyhedron),m_polyhedron);
          put(m_huvm, halfedge, uv);
        }
      }
    }

    void  set_vertex_parameterized(halfedge_descriptor halfedge, bool b) {
      if(is_border(halfedge,m_polyhedron)){
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(halfedge,m_polyhedron)){
          info(hd)->is_parameterized(b);
        }
      }else {
        info(halfedge)->is_parameterized(true);
        while(info(next(halfedge,m_polyhedron))->seaming()==0){
          halfedge = opposite(next(halfedge,m_polyhedron),m_polyhedron);
          info(halfedge)->is_parameterized(true);
        }
      }
    }
    Point_2  get_vertex_uv(halfedge_descriptor halfedge) const {
      return get(m_huvm, halfedge);
    }

    /// Get/set "is parameterized" field of vertex. Default value is undefined.
    /// (stored in halfedges sharing the same vertex).
    bool  is_vertex_parameterized(vertex_descriptor vertex) const {
      return are_corners_parameterized(vertex, vertex_descriptor(), vertex_descriptor());
    }
    void  set_vertex_parameterized(vertex_descriptor vertex, bool parameterized) {
      set_corners_parameterized(vertex, vertex_descriptor(), vertex_descriptor(), parameterized);
    }

    /// Get/set vertex index. Default value is undefined.
    /// (stored in Polyhedron vertex for debugging purpose).
    int  get_vertex_index(vertex_descriptor vertex) const {
        return info(vertex)->index();
    }

    void  set_vertex_index(vertex_descriptor vertex, int index)
    {
        info(vertex)->index(index);
    }

    int get_vertex_index(halfedge_descriptor halfedge) const
    {
      return info(halfedge)->vindex();
    }

    void  set_vertex_index(halfedge_descriptor halfedge, int index)
    {
        info(halfedge)->vindex(index);
    }


    /// Get/set vertex' all purpose tag. Default value is undefined.
    /// (stored in halfedges sharing the same vertex).
    int  get_vertex_tag(vertex_descriptor vertex) const {
      return get_corners_tag(vertex, vertex_descriptor(), vertex_descriptor());
    }
    void set_vertex_tag(vertex_descriptor vertex, int tag) {
      set_corners_tag(vertex, vertex_descriptor(), vertex_descriptor(), tag);
    }

    /// Return true if a vertex belongs to ANY mesh's border.
    bool  is_vertex_on_border(vertex_descriptor vertex) const
    {
      return info(vertex)->on_border() > 0;
    }

    /// Return true if a vertex belongs to the UNIQUE mesh's main border,
    /// i.e. the mesh's LONGEST border.
    bool  is_vertex_on_main_border(vertex_descriptor vertex) const {
        return info(vertex)->on_border() == 2;
    }

    /// Get circulator over the vertices incident to 'vertex'.
    /// 'start_position' defines the optional initial position of the circulator.
    vertex_around_target_circulator vertices_around_vertex_begin(
                            vertex_descriptor vertex,
                            vertex_descriptor start_position = vertex_descriptor())
    {
      if (start_position == vertex_descriptor())
          return vertex_around_target_circulator(halfedge(vertex,m_polyhedron));
        else
            return vertex_around_target_circulator(
                                                   get_halfedge(start_position, vertex));
    }
  

    //@} // end of ParameterizationMesh_3 INTERFACE

    //*******************************************************************
    /// @name ParameterizationPatchableMesh_3 INTERFACE
    //*******************************************************************
    //@{

    // VERTEX INTERFACE

    /// Get/set vertex seaming flag. Default value is undefined.
    int  get_vertex_seaming(vertex_descriptor vertex) const {
        return info(vertex)->seaming();
    }
    void set_vertex_seaming(vertex_descriptor vertex, int seaming) {
        info(vertex)->seaming(seaming);
    }

    // EDGE INTERFACE

    /// Get/set oriented edge's seaming flag, i.e. position of the oriented edge
    /// w.r.t. to the UNIQUE main border.
    int  get_halfedge_seaming(vertex_descriptor source, vertex_descriptor target) const {
        return info(get_halfedge(source, target))->seaming();
    }
    void set_halfedge_seaming(vertex_descriptor source, vertex_descriptor target, int seaming) {
        info(get_halfedge(source, target))->seaming(seaming);
    }

    // CORNER INTERFACE

    /// Get/set the 2D position (= (u,v) pair) of corners at the <i>right</i>
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    /// (stored in incident halfedges).
    Point_2 get_corners_uv(vertex_descriptor vertex,
                           vertex_descriptor prev_vertex,
                           vertex_descriptor next_vertex) const
    {
        // if inner vertex
      if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // get (u,v) pair from any incident halfedge
          return get(m_huvm, halfedge(vertex,m_polyhedron));
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // TODO no need for a circulator
            // get (u,v) pair from first inner halfedge (clockwise)
            halfedge_around_target_circulator cir(
                                                        get_halfedge(next_vertex, vertex), m_polyhedron );
            return get(m_huvm, *cir);
        }
    }
    void set_corners_uv(vertex_descriptor vertex,
                        vertex_descriptor prev_vertex,
                        vertex_descriptor next_vertex,
                        const Point_2& uv)
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
          // Loop over all incident halfedges
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vertex,m_polyhedron)){
            put(m_huvm, hd, uv);
          }
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // first inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex), m_polyhedron );

            // past-the-end inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir_end(
                                                      get_halfedge(prev_vertex, vertex), m_polyhedron );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
              put(m_huvm, *cir, uv);
        }
    }

    /// Get/set "is parameterized" field of corners at the <i>right</i>
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    /// (stored in incident halfedges).
    bool are_corners_parameterized(vertex_descriptor vertex,
                                   vertex_descriptor prev_vertex,
                                   vertex_descriptor next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // get "is parameterized" field from any incident halfedge
          return info(halfedge(vertex,m_polyhedron))->is_parameterized();
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // get "is parameterized" field from first inner halfedge (clockwise)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex),m_polyhedron );
            return info(*cir)->is_parameterized();
        }
    }
    void set_corners_parameterized(vertex_descriptor vertex,
                                   vertex_descriptor prev_vertex,
                                   vertex_descriptor next_vertex,
                                   bool parameterized)
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
          // Loop over all incident halfedges
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vertex,m_polyhedron)){
            info(hd)->is_parameterized(parameterized);
          }
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // first inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex),m_polyhedron );

            // past-the-end inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir_end(
                                                      get_halfedge(prev_vertex, vertex),m_polyhedron );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                info(*cir)->is_parameterized(parameterized);
        }
    }

    // AF: only used by Parameterization_mesh_patch_3::vertex_index()
    /// Get/set index of corners at the <i>right</i>
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    /// (stored in incident halfedges).
    int get_corners_index(vertex_descriptor vertex,
                          vertex_descriptor prev_vertex,
                          vertex_descriptor next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // get index from any incident halfedge
          return info(halfedge(vertex, m_polyhedron))->index();
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // get index from first inner halfedge (clockwise)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex), m_polyhedron );
            return info(*cir)->index();
        }
    }
    void set_corners_index(vertex_descriptor vertex,
                           vertex_descriptor prev_vertex,
                           vertex_descriptor next_vertex,
                           int index)
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // Loop over all incident halfedges
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vertex, m_polyhedron)){
                info(cir)->index(index);
            }
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // first inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex),m_polyhedron );

            // past-the-end inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir_end(
                                                      get_halfedge(prev_vertex, vertex),m_polyhedron );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                info(*cir)->index(index);
        }
    }

    /// Get/set all purpose tag of corners at the <i>right</i>
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    /// (stored in incident halfedges).
    int get_corners_tag(vertex_descriptor vertex,
                        vertex_descriptor prev_vertex,
                        vertex_descriptor next_vertex) const
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // get tag from any incident halfedge
          return info(halfedge(vertex, m_polyhedron))->tag();
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // get tag from first inner halfedge (clockwise)
            // TODO no need for a circulator
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex),m_polyhedron );
            return info(*cir)->tag();
        }
    }
    void set_corners_tag(vertex_descriptor vertex,
                         vertex_descriptor prev_vertex,
                         vertex_descriptor next_vertex,
                         int tag)
    {
        // if inner vertex
        if (prev_vertex == boost::graph_traits<Polyhedron>::null_vertex() && next_vertex == boost::graph_traits<Polyhedron>::null_vertex())
        {
            // Loop over all incident halfedges
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vertex,m_polyhedron)){
              info(hd)->tag(tag);
            }
        }
        else // if seam vertex
        {
            CGAL_surface_mesh_parameterization_precondition(prev_vertex != boost::graph_traits<Polyhedron>::null_vertex());
            CGAL_surface_mesh_parameterization_precondition(next_vertex != boost::graph_traits<Polyhedron>::null_vertex());

            // first inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir(
                                                  get_halfedge(next_vertex, vertex),m_polyhedron );

            // past-the-end inner halfedge (for a clockwise rotation)
            halfedge_around_target_circulator cir_end(
                                get_halfedge(prev_vertex, vertex),m_polyhedron  );

            // Loop over incident halfedges at the "right"
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                info(*cir)->tag(tag);
        }
        }

    //@} // end of ParameterizationPatchableMesh_3 INTERFACE
  #endif //DOXYGEN_RUNNING

// Private operations
private:

  struct Fct {
    Adaptor& a;
    int i;

    Fct(Adaptor& a, int i)
      : a(a), i(i)
    {}

    void operator()(vertex_descriptor vd) const
    {
      a.info(vd)->on_border(i);
    }
  };
  
  
    /// Extract mesh's longest border.
    halfedge_descriptor extract_longest_border(Polyhedron& mesh)
    {
        std::list<vertex_descriptor> result;    // returned list

        Fct fct1(*this,1);
        Fct fct2(*this,2);
        halfedge_descriptor hd = longest_border(mesh,fct1);
        if(hd != boost::graph_traits<Polyhedron>::null_halfedge()){
          BOOST_FOREACH(halfedge_descriptor hdf, halfedges_around_face(hd,mesh)){
            fct2(target(hdf,mesh));
          }
        }
        return hd;
    }
// Fields
private:

    /// The adapted mesh (cannot be NULL).
    Polyhedron&                 m_polyhedron;
  
    VertexIndexMap              m_vim;
    HalfedgeIndexMap            m_him;
    HalfedgeUvMap               m_huvm;
    /// Additional info attached to halfedges.
    Halfedge_info_map           m_halfedge_info;
    /// Additional info attached to vertices.
    Vertex_info_map             m_vertex_info;

    /// Main border of a topological disc inside m_polyhedron (may be empty).
    std::list<halfedge_descriptor> m_main_border;
    int m_num_vertices;

    }; // Parameterization_polyhedron_adaptor_3

  

  template <typename PolygonMesh, typename Fct>
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
  longest_border(const PolygonMesh& mesh, Fct fct)
  {
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor; 
    typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(CGAL::vertex_point,mesh);
    double length=0;
    halfedge_descriptor result;
    std::set<halfedge_descriptor> visited;
    // TODO use halfedge_index
    // std::vector<bool> visited(num_halfedges(mesh));
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
      if(is_border(hd,mesh)&& (visited.find(hd)== visited.end())){
        double len=0;
        BOOST_FOREACH(halfedge_descriptor hdf, halfedges_around_face(hd,mesh)){
          fct(target(hdf,mesh));
          visited.insert(hdf);
          len += std::sqrt(squared_distance(get(ppmap,source(hdf,mesh)), get(ppmap,target(hdf,mesh))));
        }
        if(len > length){
          length = len;
          result = hd;
        }
      }
    }
    return result;
  }
  
  template <typename PolygonMesh>
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
  longest_border(const PolygonMesh& mesh)
  {
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor; 
    typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(CGAL::vertex_point,mesh);
    double length=0;
    halfedge_descriptor result;
    std::set<halfedge_descriptor> visited;
    // TODO use halfedge_index
    // std::vector<bool> visited(num_halfedges(mesh));
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
      if(is_border(hd,mesh)&& (visited.find(hd)== visited.end())){
        double len=0;
        BOOST_FOREACH(halfedge_descriptor hdf, halfedges_around_face(hd,mesh)){
          visited.insert(hdf);
          len += std::sqrt(squared_distance(get(ppmap,source(hdf,mesh)), get(ppmap,target(hdf,mesh))));
        }
        if(len > length){
          length = len;
          result = hd;
        }
      }
    }
    return result;
  }
  



} //namespace CGAL

#endif //CGAL_SURFACE_MESH_PARAMETERIZATION_POLYHEDRON_ADAPTOR3_H
