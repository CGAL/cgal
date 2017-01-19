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

#ifndef CGAL_COREFINEMENT_OPERATIONS_H
#define CGAL_COREFINEMENT_OPERATIONS_H

#include <CGAL/license/Polygon_mesh_processing.h>


#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/internal/corefinement/Combinatorial_map_output_builder.h>

namespace CGAL{
/** \cond */
namespace internal{
  
template <class HDS>
class Import_volume_as_polyhedron : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Halfedge_handle Halfedge_handle;
  typedef typename HDS::Vertex_handle   Vertex_handle;
  typedef typename HDS::Face_handle     Face_handle;
  typedef typename HDS::Vertex          Vertex;
  typedef typename HDS::Halfedge        Halfedge;
  typedef typename HDS::Face            Face;
  
  //data members
  Face_handle current_face;
  std::vector< typename Vertex::Point >   points;
  std::size_t nb_edges;
  std::vector<CGAL::cpp11::tuple<unsigned,unsigned,unsigned> >           faces;

public:
  
  //to import each piece individually
  template <class Combinatorial_map_3>
  Import_volume_as_polyhedron(const Combinatorial_map_3& map,typename Combinatorial_map_3::Dart_const_handle dart):nb_edges(0)
  {
    typedef Combinatorial_map_3 CMap;
    typedef std::map<const typename Vertex::Point*,unsigned int> Vertex_map;
    Vertex_map vertex_map;
    unsigned int index=0;
    //recover all the vertices in the current volumeiterator over all the point
    //map the vertex to the index of the point in the vector
    for (typename CMap::template One_dart_per_incident_cell_const_range<0,3>::const_iterator 
         it=map.template one_dart_per_incident_cell<0,3>(dart).begin(),
         itend=map.template one_dart_per_incident_cell<0,3>(dart).end();
       it!=itend; ++it)
    {
      points.push_back(map.template attribute<0>(it)->point());
      vertex_map.insert(std::make_pair(&map.template attribute<0>(it)->point(),index++));
    }
    
    //count the number of edges    
    nb_edges+=map.template one_dart_per_incident_cell<1,3>(dart).size();

    //recover one dart per face
    for (typename CMap::template  One_dart_per_incident_cell_const_range<2,3>::const_iterator 
         it=map.template one_dart_per_incident_cell<2,3>(dart).begin(),
         itend=map.template one_dart_per_incident_cell<2,3>(dart).end();
       it!=itend; ++it)
    {
      //warning: the convention used into a polyhedron is that the normal 
      //         of a triangle indicates the outside of the object; thus 
      //         we need to reverse the orientation of the faces of the
      //         combinatorial map.
      unsigned int i=vertex_map[&map.template attribute<0>(it)->point()];
      unsigned int j=vertex_map[&map.template attribute<0>(map.beta(it, 0))->point()];
      unsigned int k=vertex_map[&map.template attribute<0>(map.beta(it, 1))->point()];
      faces.push_back(CGAL::cpp11::make_tuple(i,j,k));
    }
  }
  
  //for intersection and symetric difference
  template <class Combinatorial_map_3,class Iterator>
  Import_volume_as_polyhedron(const Combinatorial_map_3& map,
                              Iterator dart_begin,
                              Iterator dart_end):nb_edges(0)
  {
    typedef Combinatorial_map_3 CMap;
    typedef std::map<const typename Vertex::Point*,unsigned int> Vertex_map;
    Vertex_map vertex_map;
    unsigned int index=0;
    
    for (Iterator it=dart_begin;it!=dart_end;++it)
    {
      typename Combinatorial_map_3::Dart_const_handle dart=*it;
      //recover all the vertices in the current volumeiterator over all the point
      //map the vertex to the index of the point in the vector
      for (typename CMap::template  One_dart_per_incident_cell_const_range<0,3>::const_iterator 
           it=map.template one_dart_per_incident_cell<0,3>(dart).begin(),
           itend=map.template one_dart_per_incident_cell<0,3>(dart).end();
         it!=itend; ++it)
      {
        if ( vertex_map.insert(std::make_pair(&map.template attribute<0>(it)->point(),index)).second )
        {          
          points.push_back(map.template attribute<0>(it)->point());
          ++index;
        }
      }
      
       //count the number of edges    
      nb_edges+=map.template one_dart_per_incident_cell<1,3>(dart).size();

      //recover one dart per face
      for (typename CMap::template One_dart_per_incident_cell_const_range<2,3>::const_iterator 
           it=map.template one_dart_per_incident_cell<2,3>(dart).begin(),
           itend=map.template one_dart_per_incident_cell<2,3>(dart).end();
         it!=itend; ++it)
      {
        //warning: the convention used into a polyhedron is that the normal 
        //         of a triangle indicates the outside of the object; thus 
        //         we need to reverse the orientation of the faces of the
        //         combinatorial map.        
        unsigned int i=vertex_map[&map.template attribute<0>(it)->point()];
        unsigned int j=vertex_map[&map.template attribute<0>(map.beta(it, 0))->point()];
        unsigned int k=vertex_map[&map.template attribute<0>(map.beta(it, 1))->point()];
        faces.push_back(CGAL::cpp11::make_tuple(i,j,k));
      }
    }
  }
  
  //for union : use the inverse of the complementary
  template <class Combinatorial_map_3,class Iterator>
  Import_volume_as_polyhedron(const Combinatorial_map_3& map,
                              Iterator dart_begin,
                              Iterator dart_end,bool):nb_edges(0)
  {
    typedef Combinatorial_map_3 CMap;
    typedef std::map<const typename Vertex::Point*,unsigned int> Vertex_map;
    Vertex_map vertex_map;
    unsigned int index=0;
    
    for (Iterator it=dart_begin;it!=dart_end;++it)
    {
      typename Combinatorial_map_3::Dart_const_handle dart=*it;
      //recover all the vertices in the current volumeiterator over all the point
      //map the vertex to the index of the point in the vector
      for (typename CMap::template One_dart_per_incident_cell_const_range<0,3>::const_iterator 
           it=map.template one_dart_per_incident_cell<0,3>(dart).begin(),
           itend=map.template one_dart_per_incident_cell<0,3>(dart).end();
         it!=itend; ++it)
      {
        if (vertex_map.insert(std::make_pair(&map.template attribute<0>(it)->point(),index)).second )
        {
          points.push_back(map.template attribute<0>(it)->point());
          ++index;
        }
      }
      
      //count the number of edges  
      nb_edges+=map.template one_dart_per_incident_cell<1,3>(dart).size();

      //recover one dart per face
      for (typename CMap::template One_dart_per_incident_cell_const_range<2,3>::const_iterator 
           it=map.template one_dart_per_incident_cell<2,3>(dart).begin(),
           itend=map.template one_dart_per_incident_cell<2,3>(dart).end();
         it!=itend; ++it)
      {
        //warning: the convention used into a polyhedron is that the normal 
        //         of a triangle indicates the outside of the object; thus 
        //         we need to reverse the orientation of the faces of the
        //         combinatorial map. Since to get the complementary we 
        //         also need to reverse the orientation, we finally do
        //         not change it.
        unsigned int i=vertex_map[&map.template attribute<0>(it)->point()];
        unsigned int j=vertex_map[&map.template attribute<0>(map.beta(it, 1))->point()];
        unsigned int k=vertex_map[&map.template attribute<0>(map.beta(it, 0))->point()];
        faces.push_back(CGAL::cpp11::make_tuple(i,j,k));
      }
    }
  }
    
  
  void operator()( HDS& hds)
  {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface( points.size(), faces.size(),2*nb_edges);
    
    //insert vertices
    for (typename std::vector<typename Vertex::Point>::iterator it=points.begin();it!=points.end();++it)
      B.add_vertex(*it);

    //create faces
    for (std::vector<CGAL::cpp11::tuple<unsigned,unsigned,unsigned> >::iterator it=faces.begin();it!=faces.end();++it)
    {
      B.begin_facet();
      B.add_vertex_to_facet(CGAL::cpp11::get<0>(*it));
      B.add_vertex_to_facet(CGAL::cpp11::get<1>(*it));
      B.add_vertex_to_facet(CGAL::cpp11::get<2>(*it));
      B.end_facet();
    }
    B.end_surface();    
  }
};

}
/** \endcond */

/*! \class Polyhedron_corefinement corefinement_operations.h CGAL/corefinement_operations.h
  * Function object to compute the decomposition of the space induced by two polyhedra.
  * @tparam Polyhedron must be an instantiation of  CGAL::Polyhedron_3<Traits>.
  * @tparam Kernel must be a CGAL Kernel compatible with the underlying kernel of Polyhedron.
  * @tparam Output_polyhedron is a polyhedron type used as output. `Kernel` must be compatible with the underlying kernel of `Output_polyhedron`.
  */
template <class Polyhedron,class Kernel=typename Polyhedron::Traits::Kernel, class Output_polyhedron=Polyhedron>
class Polyhedron_corefinement
{
  typedef internal::Import_volume_as_polyhedron<typename Output_polyhedron::HalfedgeDS> Volume_import_modifier;
  
public:
  /** Enumeration of the different feature tags, listing all kind of space decomposition the functor can compute given two input polyhedra P and Q.*/
  enum Boolean_operation_tag
  {
    Join_tag=1,           /*!< the union of P and Q. */  
    Intersection_tag=2,   /*!< the intersection of P and Q. */  
    P_minus_Q_tag=4,      /*!< P minus Q. */  
    Q_minus_P_tag=8,      /*!< Q minus P. */  
    Parts_of_P_tag=32,    /*!< decomposition of the volume bounded by P induced by Q. */  
    Parts_of_Q_tag=64,    /*!< decomposition of the volume bounded by Q induced by P. */  
    Decomposition_tag=16  /*!< Both decompositions of P and Q. */  
  };

  /**
    * This computes different polyhedra according to the value of features.
    * Each input polyhedron is expected to bound a volume: the volume is bounded by the surface polyhedron with facet
    * normals pointing outside the volume. Each facet is supposed to be counterclockwise oriented, that is its vertex 
    * sequence (induced by halfedges) is seen counterclockwise from the side of the facet pointed by its normal. If one or
    * both polyhedra are not closed, the algorithm will end up correctly if each intersection polyline separates each surface
    * in two components. In that case for each triangle of each surface path boundary, the interior of the volume is considered
    * to be on the side opposite of that pointed by it normals (but the orientation must be consistent on each patch). 
    * the surface the volume bounded by an open surface is considered to be an infinite volume above
    * or below the surface (the side not pointed by normal vectors). 
    * @tparam Polyline_output_iterator must be an output iterator of std::vector<Kernel::Point_3>.
    * @tparam Polyhedron_ptr_and_type_output_iterator an output iterator of std::pair<Polyhedron*,int>.
    * @param P is the first input triangulated polyhedron. Note that a reference is taken and P will be updated to contain the 1D intersection between the two surfaces P and Q.
    * @param Q is the second input triangulated polyhedron. Note that a reference is taken and Q will be updated to contain the 1D intersection between the two surfaces P and Q.
    * @param polyline_output is an output iterator that collects intersection polylines between P and Q.
    * @param poly_output is an output iterator that collects output polyhedra. Note that each polyhedron is allocated within this function (thus must be explicitly deleted when no longer used). The integer is a feature tag corresponding to the volume represented by the polyhedron of the pair.
    * @param features is an integer indicating what polyhedra the function must compute. If several kind of polyhedra are expected, feature tags must be combined by |. For example if features = Polyhedron_corefinement<Polyhedron,Kernel>::Join_tag | Polyhedron_corefinement<Polyhedron,Kernel>::Intersection_tag, then poly_output will collect two polyhedra, the union and the intersection of P and Q. 
    */
  template <class Polyline_output_iterator,class Polyhedron_ptr_and_type_output_iterator>
  void operator()(  Polyhedron& P, Polyhedron& Q,
                    Polyline_output_iterator polyline_output,
                    Polyhedron_ptr_and_type_output_iterator poly_output,
                    int features) const
  {
    typedef CGAL::Corefinement::Combinatorial_map_output_builder<Polyhedron> Output_builder;
    Output_builder output_builder;
    typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron, Output_builder> Split_visitor;
    Split_visitor visitor(output_builder);
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,
      Kernel,
      Split_visitor> polyline_intersections(visitor);

    polyline_intersections(P, Q, polyline_output);

    typedef typename Output_builder::Combinatorial_map_3  Combinatorial_map_3;
    typedef typename Output_builder::Volume_info  Volume_info;
    typedef typename Combinatorial_map_3::Dart_const_handle Dart_const_handle;
    
    const Combinatorial_map_3& final_map=output_builder.combinatorial_map();
   
    typename Combinatorial_map_3::template One_dart_per_cell_const_range<3> cell_range=final_map.template one_dart_per_cell<3>();

    std::list<Dart_const_handle> intersection;
    std::list<Dart_const_handle> union_;
    std::list<Dart_const_handle> P_minus_Q;
    std::list<Dart_const_handle> Q_minus_P;
        
    for (typename Combinatorial_map_3::template One_dart_per_cell_const_range<3>::const_iterator 
      it = cell_range.begin(), it_end= cell_range.end();
      it!= it_end;
      ++it )
    {

      const Volume_info& info=final_map.template attribute<3>(it)->info();
      std::size_t inside_size=info.inside.size();
      std::size_t outside_size=info.outside.size();

      if ( inside_size + outside_size != 2){
        std::cerr << "Error: One volume cannot be represented using a polyhedron. Aborted.\n";
        break;
      }

      switch (outside_size)
      {
        case 2:
          if (features & Join_tag) union_.push_back(it);
          break;
        case 0:
          if (features & Intersection_tag) intersection.push_back(it);
          break;
        default:
          if ( *info.inside.begin() == &P )
          { if (features & P_minus_Q_tag) P_minus_Q.push_back(it); }
          else
          { if (features & Q_minus_P_tag) Q_minus_P.push_back(it); }
      }
      
      if ( features&Decomposition_tag )
      {
        Volume_import_modifier modifier(final_map,it);
        Output_polyhedron* new_poly=new Output_polyhedron();
        new_poly->delegate(modifier);
        *poly_output++ = std::make_pair( new_poly,static_cast<int>(Decomposition_tag) );
      }
      
      if ( features&Parts_of_P_tag && info.inside.find(&P)!=info.inside.end() )
      {
        Volume_import_modifier modifier(final_map,it);
        Output_polyhedron* new_poly=new Output_polyhedron();
        new_poly->delegate(modifier);
        *poly_output++ = std::make_pair( new_poly,static_cast<int>(Parts_of_P_tag) );
      }

      if ( features&Parts_of_Q_tag && info.inside.find(&Q)!=info.inside.end() )
      {
        Volume_import_modifier modifier(final_map,it);
        Output_polyhedron* new_poly=new Output_polyhedron();
        new_poly->delegate(modifier);
        *poly_output++ = std::make_pair( new_poly,static_cast<int>(Parts_of_Q_tag) );
      }
    }
    
    if (!intersection.empty())
    {
      Volume_import_modifier modifier(final_map,intersection.begin(),intersection.end());
      Output_polyhedron* new_poly=new Output_polyhedron();
      new_poly->delegate(modifier);
      *poly_output++=std::make_pair( new_poly,static_cast<int>(Intersection_tag) );
    }
    
    if (!P_minus_Q.empty())
    {
      Volume_import_modifier modifier(final_map,P_minus_Q.begin(),P_minus_Q.end());
      Output_polyhedron* new_poly=new Output_polyhedron();
      new_poly->delegate(modifier);
      *poly_output++=std::make_pair( new_poly,static_cast<int>(P_minus_Q_tag) );
    }

    if (!Q_minus_P.empty())
    {
      Volume_import_modifier modifier(final_map,Q_minus_P.begin(),Q_minus_P.end());
      Output_polyhedron* new_poly=new Output_polyhedron();
      new_poly->delegate(modifier);
      *poly_output++=std::make_pair( new_poly,static_cast<int>(Q_minus_P_tag) );
    }
    
    if (!union_.empty())
    {
      Volume_import_modifier modifier(final_map,union_.begin(),union_.end(),true);
      Output_polyhedron* new_poly=new Output_polyhedron();
      new_poly->delegate(modifier);
      *poly_output++=std::make_pair( new_poly,static_cast<int>(Join_tag) );
    }
  }
  /**
    * This updates P according to the value of features.
    * Each input polyhedron is expected to bound a volume: the volume is bounded by the surface polyhedron with facet
    * normals pointing outside the volume. Each facet is supposed to be counterclockwise oriented, that is its vertex 
    * sequence (induced by halfedges) is seen counterclockwise from the side of the facet pointed by its normal. If one or
    * both polyhedra are not closed, the algorithm will end up correctly if each intersection polyline separates each surface
    * in two components. In that case for each triangle of each surface path boundary, the interior of the volume is considered
    * to be on the side opposite of that pointed by it normals (but the orientation must be consistent on each patch). 
    * the surface the volume bounded by an open surface is considered to be an infinite volume above
    * or below the surface (the side not pointed by normal vectors). 
    * @tparam Polyline_output_iterator must be an output iterator of std::vector<Kernel::Point_3>.
    * @param P is the first input triangulated polyhedron. Note that a reference is taken and P will be updated to contain the 1D intersection between the two surfaces P and Q.
    * @param Q is the first input triangulated polyhedron. Note that a reference is taken and Q will be updated to contain the 1D intersection between the two surfaces P and Q.
    * @param polyline_output is an output iterator that collects intersection polylines between P and Q.
    * @param features is either Join_tag, Intersection_tag, P_minus_Q_tag or Q_minus_P_tag
    */
  template <class Polyline_output_iterator>
  void operator()(  Polyhedron& P, Polyhedron& Q,
                    Polyline_output_iterator polyline_output,
                    Boolean_operation_tag features) const
  {
    typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron> Split_visitor;
    Split_visitor visitor;
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,
      Kernel,
      Split_visitor> polyline_intersections(visitor);

    polyline_intersections(P, Q, polyline_output);

    typedef typename Split_visitor::Combinatorial_map_3  Combinatorial_map_3;
    typedef typename Split_visitor::Volume_info  Volume_info;
    typedef typename Combinatorial_map_3::Dart_const_handle Dart_const_handle;
    
    const typename Split_visitor::Combinatorial_map_3& final_map=visitor.combinatorial_map();
   
    typename Combinatorial_map_3::template One_dart_per_cell_const_range<3> cell_range=final_map.template one_dart_per_cell<3>();

    std::list<Dart_const_handle> darts;
        
    for (typename Combinatorial_map_3::template One_dart_per_cell_const_range<3>::const_iterator 
      it = cell_range.begin(), it_end= cell_range.end();
      it!= it_end;
      ++it )
    {

      const Volume_info& info=it->template attribute<3>()->info();
      std::size_t inside_size=info.inside.size();
      std::size_t outside_size=info.outside.size();

      if ( inside_size + outside_size != 2){
        std::cerr << "Error: One volume cannot be represented using a polyhedron. Aborted.\n";
        break;
      }

      switch (outside_size)
      {
        case 2:
          if ( features==Join_tag ) darts.push_back(it);
          break;
        case 0:
          if ( features==Intersection_tag ) darts.push_back(it);
          break;
        default:
          if ( *info.inside.begin() == &P )
          { if (features == P_minus_Q_tag) darts.push_back(it); }
          else
          { if (features == Q_minus_P_tag) darts.push_back(it); }
      }
    }
    
    P.clear();
    
    if (!darts.empty())
    {
      Volume_import_modifier modifier=
      features==Join_tag?
          Volume_import_modifier(final_map,darts.begin(),darts.end(),true)
        : Volume_import_modifier(final_map,darts.begin(),darts.end());
      P.delegate(modifier);
    }
  }
  
  /** \cond */
  static std::string get_type_str(const std::string& Pname,const std::string& Qname,int i)
  {
    switch (i)
    {
      case Join_tag: 
        return Pname+std::string("_union_")+Qname;
      case P_minus_Q_tag: 
        return Pname+std::string("_minus_")+Qname;
      case Q_minus_P_tag:
        return Qname+std::string("_minus_")+Pname;
      case Intersection_tag:
        return Pname+std::string("_inter_")+Qname;
      case Decomposition_tag:
        return std::string("Decomposition");
    }
    return std::string("Unknow");
  }
  /** \endcond */ 
};




}

#endif // CGAL_COREFINEMENT_OPERATIONS_H
