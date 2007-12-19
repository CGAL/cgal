// Copyright (c) 2005, 2006  INRIA Sophia-Antipolis (France).
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
// author(s)     : Raphaelle Chaine (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//                 Jerome Piovano
//                 Andreas Fabri (GeometryFactory)


#ifndef CGAL_RECONSTRUCTION_FROM_SLICES_3_H
#define CGAL_RECONSTRUCTION_FROM_SLICES_3_H

#include <iostream>
#include <vector>
#include <CGAL/IO/VRML_2_ostream.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_from_slices_3.h>
#include <CGAL/TFS_polyline_vertex_base_3.h>
#include <CGAL/TFS_cell_base_3.h>

#ifdef GEOMVIEW_DUMP //TO RM
#include <CGAL/IO/Geomview_stream.h>
#endif // GEOMVIEW_DUMP TO RM

CGAL_BEGIN_NAMESPACE

typedef Exact_predicates_exact_constructions_kernel default_Gt;
typedef Triangulation_from_slices_3<default_Gt, 
  Triangulation_data_structure_3< TFS_polyline_vertex_base_3< Triangulation_vertex_base_3<default_Gt> >, 
  TFS_cell_base_3<Triangulation_cell_base_3<default_Gt> > > > default_Tr;

// template < class Tr >   = Triangulation_from_slices_3< Gt, 
//           Triangulation_data_structure_3< TFS_polyline_vertex_base_3< Triangulation_vertex_base_3<Gt> >, 
//                                           TFS_cell_base_3<Triangulation_cell_base_3<Gt> > > 
// Gt (default) = Exact_predicates_exact_constructions_kernel


/*===========================================================================*/
/*         FUNCTIONS LIST                                                    */
/*===========================================================================*/
template < class Tr >
bool check_polylines_respect(const Tr & tr);

template < class Tr >
void check_and_conform_polylines(Tr & tr);
// Guaranteed only for non intersecting polylines without small angles

template < class Tr >
void tag_inner_outer(const Tr & tr);
// Combinatorial tag of tr cells, extends 2D inner/outer property with regards to polylines 
//                               into 3D inner/outer property with regards to an induced volume (cf. Nuage Boissonnat/Geiger)
// tr is a triangulation constrained to slices and the vertices must be organised in polylines 
// Precondition : tr must not have been tagged previously, check_polylines_respect(tr)==true                    

template <class Tr>
void non_solid_connections_removal(Tr & tr, bool heavy_removal=false);
// Removes non manifoldiness of the surface enclosed between inner and outer cells
//      -- if heavy_removal=false : removes cluster of facets between two slides
//         that are not solidly attached ot the upper and the lower slice
//      -- if heavy_removal=true then stable_non_solid_connections_heavy_removal 
//         should be call afterwards
// Precondition : tag_inner_outer(tr) previously called

template <class Tr>
void stable_non_solid_connections_heavy_removal(Tr & tr);
// Removes further local non manifoldiness of the surface enclosed between inner and outer cells
//   (those located at the level of single edges)
// Precondition : tag_inner_outer(tr) previously called

template <class Tr>
void save_in_off_file(const Tr & tr, const char * fname_head_off, const char * fname_body_off);
//precondition : tag_inner_outer(tr) (and possibly non_solid_connections_removal(tr)) previously called
//postcondition : cat of the two files remain to be done ...

template < class Tr >
typename Tr::Facet get_facet_towards_external_infinite_cell(const Tr & tr);
// Returns a facet oriented towards an external infinite cell, 
// the facet has its 3 vertices located on 2 different slices (here the first and the second one)

/*===========================================================================*/
/*         INTERNAL FUNCTIONS LIST                                           */
/*===========================================================================*/

template < class Tr, class OutputIterator>
static typename Tr::Facet get_facet_towards_external_infinite_cell(const Tr & tr,typename Tr::Cell_handle cell, OutputIterator cells);
// Internal recursive function, used to implement get_facet_towards_external_infinite_cell(tr)

template < class Tr > 
static void tag_inner_outer(const Tr & tr,typename Tr::Facet f_start,bool intern_tag, int slice_down, typename Tr::Facet & for_next_slice);
// Internal recursive function, used to implement tag_inner_outer(tr)
//
// Precondition : f_start is oriented either towards an infinite cell 
//                or towards a cell having 3 vertices on slice number slice_down+1 and 0 on slice number slice_down
//                intern_tag true if f_start is enclosed inside a polyline, false else
// Postcondition : the process can be relaunched on slice slice_down+1, from facet for_next_slice

template <class Tr>
static void determine_solidity(Tr & tr, typename Tr::Cell_handle c, bool heavy_removal=true);
// Internal function, used to implement non_solid_connections_removal(Tr & tr, bool heavy_removal=true)
// 
// Determines the solidity of the cluster of internal cells containing c, with regards to its upper slice first 
//       and with regards to its lower slice then
// (Cluster through facet adjacency between cells having vertices on the upper slice first, 
//       and with regards to its lower slice then)
// Precondition : c is internal

template <class Tr, class Iterator >
static bool test_cluster_solidity(Tr & tr, typename Tr::Cell_handle ch_start, int slice_index, bool & solid, Iterator it, 
				  int tet, int li, int lj, int other_slice_index, bool heavy_removal=true);
// Internal recursive function, used to implement determine_solidity(Tr & tr, typename Tr::Cell_handle c, bool heavy_removal=true)
// 
// Recursively determines a cluster of internal cells having vertices in slice_index slice and including ch_start,
// The cluster is pushed in Iterator it
// tet, li, lj, other_slice_index : additional information on ch_start (tetrahedra type)
//
//  Preconditions :                 ch_start is Internal
//                                  tr.slices[slice_index].test_tetrahedra_type(ch_start,tet,li,lj,other_slice_index)
//                                  The solidity of ch_start with regard to slice_index slice has not been previously determined
// Return the boolean value of solid :
//         true if the cluster contains tetrahedra having 3 vertices in slice_index slice
//         false else

template <class Tr>
static bool determine_heavy_solidity(Tr & tr, typename Tr::Cell_handle c, 
				     int slice1, int li1, int lj1,
				     int slice2, int li2, int lj2);

template <class Tr, class IteratorTet2, class IteratorTet3>
static bool test_heavy_solidity_around_edge(Tr & tr, typename Tr::Cell_handle ch_start, int slice_index, bool & heavy_solid, 
					    IteratorTet2 itTet2, IteratorTet3 itTet3,
					    int tet, int li, int lj, int lk, int ll, int other_slice_index);
// Recursively determines a cluster of internal cells having a common edge ch_start->vertex(li), ch_start->vertex(lj) and including ch_start.
// The cluster is pushed in IteratorTet2 itTet2 or IteratorTet3 itTet3 depending on the type of the tetrahedra
// 
// tet, other_slice_index : additional information on ch_start (tetrahedra type)
// lk and lh : indexes of the two remaining vertices
//
// Preconditions :                 ch_start is Internal
//                                 ch_start->vertex(li)->slice()==ch_start->vertex(lj)->slice()==slice_index
//                                 ch_start->vertex(lk)->slice() ou ch_start->vertex(ll)->slice() == other_slice_index 
//                                 The heavy solidity of ch_start with regard to slice_index slice has not been previously determined
//
// Return the boolean value of heavy_solid :
//         true if the cluster contains tetrahedra having 3 vertices in slice_index slice
//         false else
//
// Poscondition : if heavy_solid==false, the set of cells pushed in itTet2 must be set non_solid (none where pushed in itTet3)
//                if heavy_solid==true, the set of cells pushed in itTet3 must be unset with regards to heavy solid determination.

/*===========================================================================*/
/*         IN/OUT TAGGING                                                    */
/*===========================================================================*/
template < class Tr >
bool check_polylines_respect(const Tr & tr)
{
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;  
  for(typename Tr::Finite_vertices_iterator it=tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
    {
      Cell_handle ctemp;
      int i,j;
      Vertex_handle vh1(it),vh2(it->next());
      if (!tr.is_edge(vh1,vh2,ctemp,i,j)){
	std::cerr << "Error: Cannot find constraint edge " << vh1->point() << " [" << vh1->slice() << "]   " << vh2->point() << " [" << vh1->slice() << "]" << std::endl;
	return false;
      }
    }
  return true;
}

#ifdef CGAL_RFS_CONFORM // If no access to mesh_2 conforming (for example when slices are not horizontal)
template < class Tr >
void check_and_conform_polylines(Tr & tr)
{  
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle; 
  
  for (unsigned int i=0 ; i<tr.represented_slices_number ; i++)
    {
      bool conform;
      do
	{
	  conform=true;
	  for (typename std::vector<Vertex_handle>::iterator it =tr.slices[i].set_of_vertices.begin();
	       it != tr.slices[i].set_of_vertices.end(); ++it)
	    {
	      Cell_handle ctemp;
	      int ci,cj; 
	      Vertex_handle vh1(*it),vh2(vh1->next());
	      if (!tr.is_edge(vh1,vh2,ctemp,ci,cj))
		{ //The polyline going through vh1 and vh2 is not conform
		  conform=false;
#ifdef CGAL_DUMP
		  std::cout<<"Slice " << i << "is not conformal"<< std::endl;
		  std::cout<<"Do you want to conform it?"<< std::endl;		  
		  std::getchar();
#endif
		  // Conform the entire edge  vh1, vh2
		  Vertex_handle vstep1=vh1;
		  Vertex_handle vstep2;
		  do //circulate over the polyline between vh1 and vh2
		    {
		      vstep2=vstep1->next();
		      if (!tr.is_edge(vstep1,vstep2,ctemp,ci,cj))
			{
#ifdef CGAL_DUMP		  
			  std::cout << "Conform an edge by midpoint insertion" << std::endl;
#endif
			  typename Tr::Point P;
			  P=midpoint(vstep1->point(),vstep2->point());
			  Vertex_handle v=tr.insert(P,vstep1->slice());
			  vstep1->set_next(v);
			  v->set_next(vstep2);
			}
		      else
			vstep1=vstep2;
		    }while(vstep1!=vh2);
		  // the polyline between vh1 and vh2 is conformed
		  break;//for each vertex) 
		}
	    }//for each vertex
	}while(!conform); 
#ifdef CGAL_DUMP
      std::cout << "Slice is conformed" << std::endl;      
#endif
    }
}
#endif // CGAL_RFS_CONFORM

template < class Tr > 
void tag_inner_outer(const Tr & tr)
{
  CGAL_triangulation_precondition((tr.dimension()>=3)); 
  //--------------------------------------------------------------------------
#ifdef GEOMVIEW_DUMP //TO RM
  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Segment Segment;

  gv.set_edge_color (CGAL::Color(0,0,0));
  //gv.set_face_color (CGAL::Color(255,255,255));
  gv.set_bg_color (CGAL::Color(180,180,180));

  for (Finite_edges_iterator it=tr.finite_edges_begin(); it != tr.finite_edges_end(); ++it)
    {
      if (((*it).first->vertex((*it).second)->next()==(*it).first->vertex((*it).third))
	  ||((*it).first->vertex((*it).third)->next()==(*it).first->vertex((*it).second)))
	gv.set_edge_color (CGAL::Color(0,0,180));
      else
	gv.set_edge_color (CGAL::Color(0,0,0));
      gv << Segment((*it).first->vertex((*it).second)->point(),(*it).first->vertex((*it).third)->point());
    }
  gv.look_recenter();
  std::getchar();
#endif // GEOMVIEW_DUMP TO RM
  //--------------------------------------------------------------------------
  typename Tr::Facet f_infinite=get_facet_towards_external_infinite_cell(tr); 
#ifdef CGAL_DUMP 
  std::cout << "first facet_towards_external_infinite_cell found" << std::endl;
#endif
  typename Tr::Facet f_infinite_next(NULL,0);
  int in=0;
  while (f_infinite.first!=NULL)
    {
#ifdef GEOMVIEW_DUMP //To RM
      std::cout << "IN/OUT Tagging slice " << in << std::endl;
      gv.clear();
      gv.set_edge_color(CGAL::Color(180,0,0));    
      typedef typename Tr::Segment Segment;
      for (int i=0;i<4;i++)
	{ 
	  for (int j=i+1;j<4;j++)
	    gv << Segment( f_infinite.first->vertex(i)->point(), f_infinite.first->vertex(j)->point() );
	}      
      std::getchar(); 
      std::cout << " ...from an infinite exterior cell" << std::endl;
#endif // GEOMVIEW_DUMP TO RM
      tag_inner_outer(tr,f_infinite,false,in,f_infinite_next);
#ifdef CGAL_DUMP 
      std::cout << "IN/OUT Tagging slice" << in << " over" << std::endl;
#endif // CGAL_DUMP
      f_infinite=f_infinite_next;
      f_infinite_next.first=NULL;
      in++;
    }      
}

template < class Tr > 
static void tag_inner_outer(const Tr & tr, typename Tr::Facet f_start, bool intern_tag, int slice_down, typename Tr::Facet & for_next_slice)
{
  CGAL_triangulation_precondition(((tr.is_infinite(f_start.first->vertex(f_start.second)))&&(tr.slices[slice_down].tetrahedra_type(f_start.first)>0)) 
				  ||(tr.slices[slice_down+1].tetrahedra_type(f_start.first)==3));
  CGAL_triangulation_precondition((!intern_tag)
			       ||(tr.slices[slice_down].tetrahedra_type(f_start.first)==3)
			       ||(tr.slices[slice_down+1].tetrahedra_type(f_start.first)==3));
  CGAL_triangulation_precondition(!f_start.first->is_explicit_in_out_tagged());
#ifdef GEOMVIEW_DUMP //TO RM
  gv.set_edge_color(CGAL::Color(0,180,0));     
  typedef typename Tr::Segment Segment;
  for (int i=0;i<3;i++) //for each edge of the facet
    { 
      typename Tr::Edge edge=tr.get_opposite_edge_in_facet(f_start,i);
      gv << Segment( edge.first->vertex(edge.second)->point(), edge.first->vertex(edge.third)->point());
    }
  std::cout << "Tag (IN/OUT) to be transmitted " << intern_tag << std::endl;
  std::getchar(); 
#endif // GEOMVIEW_DUMP TO RM
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Cell_circulator Cell_circulator;

  if(tr.is_infinite(f_start.first))
    {
      f_start.first->set_explicit_external();
      if((tr.slices[slice_down].tetrahedra_type(f_start.first)==3) // on first slice
	 ||(tr.slices[slice_down+1].tetrahedra_type(f_start.first)==3)) // on last slice
	{
	  if(intern_tag)
	    f_start.first->neighbor(f_start.second)->set_explicit_internal();
	  else
	    f_start.first->neighbor(f_start.second)->set_explicit_external();	    	
	}
    }
  else
    {
      CGAL_triangulation_assertion(!tr.is_infinite(f_start.first->neighbor(f_start.second)));
      CGAL_triangulation_assertion((tr.slices[slice_down].tetrahedra_type(f_start.first)==3)
				   ||(tr.slices[slice_down+1].tetrahedra_type(f_start.first)==3));
      if(intern_tag)
	{
	  f_start.first->set_explicit_internal();
	  f_start.first->neighbor(f_start.second)->set_explicit_internal();
	}
      else
	{
	  f_start.first->set_explicit_external();
	  f_start.first->neighbor(f_start.second)->set_explicit_external();
	}
    } 
  
  Facet frec[3];
  frec[0]=frec[1]=frec[2]=f_start;
  for (int i=0;i<3;i++) //for each edge of the facet
    { 
      Edge edge=tr.get_opposite_edge_in_facet(f_start,i);
      bool boundary=false;
      if (edge.first->vertex(edge.second)->slice()==edge.first->vertex(edge.third)->slice())
	{
	  boundary=((edge.first->vertex(edge.second)->next()==edge.first->vertex(edge.third))
		    ||(edge.first->vertex(edge.third)->next()==edge.first->vertex(edge.second)));
	  int sindex=edge.first->vertex(edge.second)->slice();
	  Cell_circulator circ = tr.incident_cells(edge,f_start.first),endcirc(circ);
	  ++circ;
	  while(circ!=endcirc) 
	    {
	      Cell_handle c(circ);
	      int li, lj, other_slice;
	      int tet=tr.slices[sindex].tetrahedra_type(c,li,lj,other_slice);
	      if (tet==3)
		{
		  if(((sindex==slice_down)&&(slice_down+1!=other_slice))
		     ||((sindex==slice_down+1)&&(slice_down!=other_slice)))
		    /* 		  if(((sindex==slice_down)&&(tr.slices[slice_down+1].tetrahedra_type(c)==0)) */
		    /* 		     ||((sindex==slice_down+1)&&(tr.slices[slice_down].tetrahedra_type(c)==0))) */
		    frec[i]=Facet(c,li);
		}
	      else 
		{ //tet==2
		  CGAL_triangulation_precondition(tet==2);
		  if(tr.is_infinite(c))
		    {
		      if((for_next_slice.first==NULL)&&(sindex==slice_down+1)&&(slice_down!=other_slice)) //&&(tr.slices[slice_down].tetrahedra_type(c)==0))
			{
			  for_next_slice=Facet(c,c->index(tr.infinite_vertex()));
#ifdef CGAL_DUMP 
			  std::cout << "A for_next_slice candidate has been found"<< std::endl;
#endif //CGAL_DUMP
			}
		    }
		  else if (intern_tag||boundary)
		    {
		      if(((sindex==slice_down)&&(other_slice==slice_down+1)) //&&(tr.slices[slice_down+1].tetrahedra_type(c)==2))
			 ||((sindex==slice_down+1)&&(other_slice!=slice_down))) //&&(tr.slices[slice_down].tetrahedra_type(c)==0)))
			{
			  if (c->is_T2_2_down_internal())
			    break;
			  else
			    c->set_T2_2_down_internal();
			}
		      else
			{
			  if (c->is_T2_2_up_internal())
			    break;
			  else
			    c->set_T2_2_up_internal();
			}		    
		    }
		} //tet==2
	      ++circ;
	    }//while
/* 	  if (boundary) */
/* 	    intern_tag=!intern_tag; */
	}
      else //(e.first->vertex(e.second)->slice()!=e.first->vertex(e.third)->slice())
	{
	  CGAL_triangulation_assertion(tr.is_infinite((f_start.first)->vertex(f_start.second)));
	  Cell_handle c_i=f_start.first->neighbor(tr.vertex_triple_index(f_start.second,i));
	  frec[i]=Facet(c_i,c_i->index(tr.infinite_vertex()));
	  CGAL_triangulation_assertion(tr.is_infinite(frec[i].first->vertex(frec[i].second)));	
	  CGAL_triangulation_assertion(!intern_tag); 
	}	  
      if(!frec[i].first->is_explicit_in_out_tagged())
	tag_inner_outer(tr,frec[i],(boundary?!intern_tag:intern_tag),slice_down,for_next_slice);
      //CGAL_triangulation_assertion(frec[i]!=f_start);   NO   
    }
}

template < class Tr >
typename Tr::Facet get_facet_towards_external_infinite_cell(const Tr & tr)
{
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;
  Cell_handle cell=NULL;
  std::vector<Cell_handle> temp_cells;
  tr.incident_cells(tr.slices[0].get_vertex_start(),std::back_inserter(temp_cells));

  for (typename std::vector<Cell_handle>::iterator cit = temp_cells.begin();cit != temp_cells.end(); ++cit)
    {
      if(tr.is_infinite(*cit))
	{
	  cell=*cit;
	  break;
	}
    }

  CGAL_triangulation_assertion(cell!=NULL);
  std::vector<Cell_handle> visited_cells;
  Facet res=get_facet_towards_external_infinite_cell(tr,Facet(cell,cell->index(tr.infinite_vertex())),std::back_inserter(visited_cells));
  for(typename std::vector<Cell_handle>::iterator cit = visited_cells.begin();cit != visited_cells.end(); ++cit) 
    (*cit)->set_in_conflict_flag(0);

  if(res.first==NULL)
    { 
      //TO CHANGE LATER : ....
      std::cout << "TO DO The research of an external infinite cell will be improved" << std::endl;
      std::vector<Cell_handle> infinite_cells;
      tr.incident_cells(tr.infinite_vertex(),std::back_inserter(infinite_cells));
      for (typename std::vector<Cell_handle>::iterator cit = infinite_cells.begin(); cit != infinite_cells.end(); ++cit)
	{
	  CGAL_triangulation_assertion(tr.is_infinite(*cit));
	  int ind=(*cit)->index(tr.infinite_vertex());
	  
	  int i1=(*cit)->vertex((ind+1)%4)->slice();
	  int i2=(*cit)->vertex((ind+2)%4)->slice();
	  int i3=(*cit)->vertex((ind+3)%4)->slice();
	  
	  if ((i1==i2)!=(i1==i3))
	    {
	      if (((i1==0)&&((i2==1)||(i3==1)))
		  ||((i1==1)&&((i2==0)||(i3==0))))
		return Facet(*cit,ind);
	    }	  
	}
    }      
  CGAL_triangulation_assertion((res.first!=NULL)&&(tr.is_infinite(res.first))); 
  return res;
}

template < class Tr, class OutputIterator>
static typename Tr::Facet get_facet_towards_external_infinite_cell(const Tr & tr,typename Tr::Facet f, OutputIterator cellsit)
{ 
  CGAL_triangulation_precondition(tr.is_infinite(f.first->vertex(f.second)));

  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Facet Facet;

  Facet res(NULL,0);

  Cell_handle cell=f.first;
  int ind = f.second; 

  // if the 3 non infinite vertices are not on the same slice ...
  if(!( (cell->vertex((ind+1)%4)->slice() == cell->vertex((ind+2)%4)->slice()) && 
	(cell->vertex((ind+1)%4)->slice() == cell->vertex((ind+3)%4)->slice()) ) )
    return f;
  else
    {  
      cell->set_in_conflict_flag(1);
      *(cellsit++)=cell;
      for(int i=1;i<4;i++)
	{
	  Cell_handle next=cell->neighbor((ind+i)%4);
	  if (next->get_in_conflict_flag() != 0)
	    {
	      res=get_facet_towards_external_infinite_cell(tr,Facet(next,next->index(tr.infinite_vertex())),cellsit);
	      if (res.first!=NULL)
		return res;
	    }
	}
    }
  return res;
}

/*===========================================================================*/
/*         NON SOLID CELLS REMOVAL                                           */
/*===========================================================================*/

template <class Tr, class Iterator>
static bool test_cluster_solidity(Tr & tr, typename Tr::Cell_handle ch_start, int slice_index, bool & solid, Iterator it, 
				  int tet, int li, int lj, int other_slice_index, bool heavy_removal)
{
  CGAL_triangulation_precondition(ch_start->is_internal()&&(tr.slices[slice_index].test_tetrahedra_type(ch_start,tet,li,lj,other_slice_index))
				  &&(((other_slice_index==slice_index+1)&&(!ch_start->is_solid_tested_down()))
				     ||((other_slice_index==slice_index-1)&&(!ch_start->is_solid_tested_up()))));
  if(tet==3)
    solid = true; //The cluster of tetrahedra pushed in Iterator "it" is solid
  
  if(other_slice_index==slice_index+1)
    ch_start->set_solid_tested_down();   
  else // other_slice_index==slice_index-1
    ch_start->set_solid_tested_up();
  
  *(it++)=ch_start;

  typedef typename Tr::Cell_handle Cell_handle;
  Cell_handle next;
  
  for(int i=0; i<4; i++)
    {
      next = ch_start->neighbor(i);  
      if(next->is_internal()&&((tet=tr.slices[slice_index].tetrahedra_type(next,li,lj,other_slice_index))!=0)) //!tr.is_infinite(next)&&
	if(((other_slice_index==slice_index+1)&&(!next->is_solid_tested_down()))
	   ||((other_slice_index==slice_index-1)&&(!next->is_solid_tested_up())))
	  {
	    if ((heavy_removal)&&(tet==2))
	      {		
		std::pair<int,int> opp_indexes=get_opposite_index_pair(li,lj);
		determine_heavy_solidity(tr,next,slice_index,li,lj,other_slice_index,opp_indexes.first,opp_indexes.second);
	      }
	    if(!next->is_non_solid())
	      test_cluster_solidity(tr,next,slice_index,solid,it,
				    tet,li,lj,other_slice_index,heavy_removal);
	  }
    }
  return solid;
}

template <class Tr>
static bool determine_heavy_solidity(Tr & tr, typename Tr::Cell_handle c, 
				     int slice1, int li1, int lj1,
				     int slice2, int li2, int lj2)
{
  CGAL_triangulation_precondition(c->is_internal());
  CGAL_triangulation_precondition((c->vertex(li1)->slice()==slice1)&&(c->vertex(lj1)->slice()==slice1));
  CGAL_triangulation_precondition((c->vertex(li2)->slice()==slice2)&&(c->vertex(lj2)->slice()==slice2));
  typedef typename Tr::Cell_handle Cell_handle;
  
  int tempindex=slice1, templi=li1, templj=lj1, templk=li2, templl=lj2, tempother=slice2;
  for(int cmpt=0;cmpt<2;cmpt++)
    {
      if(!c->is_non_solid())	      
	if(((tempother==tempindex+1)&&(!c->is_heavy_solid_tested_down()))
	   ||((tempother==tempindex-1)&&!c->is_heavy_solid_tested_up()))
	  // c not heavy solid determined with regard to tempindex slice
	  {
	    std::vector<Cell_handle> hCellGroupTet2;	  
	    typename std::vector<Cell_handle>::iterator hitTet2;
	    std::vector<Cell_handle> hCellGroupTet3;	  
	    typename std::vector<Cell_handle>::iterator hitTet3;
	    bool heavy_solid=false;	
	    if(!test_heavy_solidity_around_edge(tr,c,tempindex,heavy_solid,
						std::back_inserter(hCellGroupTet2),std::back_inserter(hCellGroupTet3),
						2,templi,templj,templk,templl,tempother))
	      {
		CGAL_triangulation_assertion(hCellGroupTet3.begin()== hCellGroupTet3.end());
		for(hitTet2 = hCellGroupTet2.begin(); hitTet2 != hCellGroupTet2.end(); hitTet2++) 
		  (*hitTet2)->set_non_solid();
	      }
	    else //hitTet3 iterator over a non empty set		  
	      for(hitTet3 = hCellGroupTet3.begin(); hitTet3 != hCellGroupTet3.end(); hitTet3++) 
		{
		  (*hitTet3)->unset_heavy_solid_tested_down();
		  (*hitTet3)->unset_heavy_solid_tested_up();			  
		}
	    for(hitTet2 = hCellGroupTet2.begin(); hitTet2 != hCellGroupTet2.end(); hitTet2++)
	      { 
		(*hitTet2)->unset_heavy_solid_tested_down();
		(*hitTet2)->unset_heavy_solid_tested_up();
		/* 		if((*hitTet2)->is_internal()) */
		/* 		  { */
		/* 		    (*hitTet2)->unset_heavy_solid_tested_down(); */
		/* 		    (*hitTet2)->unset_heavy_solid_tested_up(); */
		/* 		  } */
	      }
	  }
      tempindex=slice2, templi=li2, templj=lj2, templk=li1, templl=lj1, tempother=slice1;
    }
  return (!c->is_non_solid());
}

template <class Tr>
static void determine_solidity(Tr & tr, typename Tr::Cell_handle c, bool heavy_removal)
// Determines the solidity of the cluster of internal cells containing c, with regards to its up and down slices.
{
  typedef typename Tr::Cell_handle Cell_handle;
  CGAL_triangulation_precondition(c->is_internal());//!tr.is_infinite(c)

  if((!c->is_solid_tested_up())||(!c->is_solid_tested_down())) //&&(!c->is_non_solid())
    {
      int index_slice=c->vertex(0)->slice(), other_slice_index;	  
      int tet1, li1, lj1, tet2, li2, lj2, temp;
      int tempindex, temptet, templi, templj, tempother;
      
      tet1=tr.slices[index_slice].tetrahedra_type(c,li1,lj1,other_slice_index);
      tet2=tr.slices[other_slice_index].tetrahedra_type(c,li2,lj2,temp);
      CGAL_triangulation_assertion((temp==index_slice)&&(tet1+tet2==4));
      
      if ((heavy_removal)&&(tet1==2))
	{
	  determine_heavy_solidity(tr,c,index_slice,li1,lj1,other_slice_index,li2,lj2);
	}
      
      tempindex=index_slice, temptet=tet1, templi=li1, templj=lj1, tempother=other_slice_index;      
      for(int cmpt=0;cmpt<2;cmpt++)
	{
	  if(!c->is_non_solid())	         
	    if(((tempother==tempindex+1)&&!c->is_solid_tested_down())
	       ||((tempother==tempindex-1)&&!c->is_solid_tested_up()))
	      //c not solid determined with regard to tempindex slice
	      {	
		std::vector<Cell_handle> cellGroup;
		typename std::vector<Cell_handle>::iterator it;
		bool solid = false;  
		if (!test_cluster_solidity(tr,c,tempindex,solid,std::back_inserter(cellGroup), 
					   temptet,templi,templj,tempother,heavy_removal))
		  for(it = cellGroup.begin(); it != cellGroup.end(); it++) 
		    {
		      (*it)->set_non_solid();
		      CGAL_triangulation_assertion(tr.slices[tempindex].tetrahedra_type(*it)!=3);//TO RM
		    }
	      }
	  tempindex=other_slice_index, temptet=tet2, templi=li2, templj=lj2, tempother=index_slice;      
	}
    }
}

template <class Tr, class IteratorTet2, class IteratorTet3>
static bool test_heavy_solidity_around_edge(Tr & tr, typename Tr::Cell_handle ch_start, int slice_index, bool & heavy_solid, 
					    IteratorTet2 itTet2, IteratorTet3 itTet3,
					    int tet, int li, int lj, int lk, int ll, int other_slice_index)
// Recursively determines a cluster of internal cells having a common edge ch_start->vertex(li), ch_start->vertex(lj) and including ch_start.
{
  CGAL_triangulation_precondition(ch_start->is_internal()&&(tr.slices[slice_index].test_tetrahedra_type(ch_start,tet)));
  CGAL_triangulation_precondition((ch_start->vertex(li)->slice()==slice_index)&&(ch_start->vertex(lj)->slice()==slice_index)); 
  CGAL_triangulation_precondition((ch_start->vertex(lk)->slice()==other_slice_index)||(ch_start->vertex(ll)->slice()==other_slice_index)); 
  CGAL_triangulation_precondition(((other_slice_index==slice_index+1)&&(!ch_start->is_heavy_solid_tested_down()))
				  ||((other_slice_index==slice_index-1)&&(!ch_start->is_heavy_solid_tested_up())));
      
  if(other_slice_index==slice_index+1)
    ch_start->set_heavy_solid_tested_down();   
  else // other_slice_index==slice_index-1
    ch_start->set_heavy_solid_tested_up();  
  
  if(tet==3)
    {
      heavy_solid = true; //The cluster of tetrahedra pushed in Iterator "it" is heavy solid
      *(itTet3++) = ch_start;
    }
  else//(tet==2)
    *(itTet2++)=ch_start;
  
  typename Tr::Cell_handle next;  
  int lnext=lk, lcommon=ll;
  
  for(int cmpt=0;cmpt<2;cmpt++)
    { 
      next=ch_start->neighbor(lnext); 
      if (next->is_internal())
	{
	  int ntet=2, nli=-1, nlj=-1, nlk=-1, nll=-1, nother=-1;
	  for(int i=0; i<4; i++)
	    {
	      if (next->vertex(i)==ch_start->vertex(li))
		nli=i;
	      else if (next->vertex(i)==ch_start->vertex(lj))
		nlj=i;
	      else 
		{ 
		  if (next->vertex(i)==ch_start->vertex(lcommon))
		    nlk=i;
		  else
		    nll=i;
		  if (next->vertex(i)->slice()==slice_index)
		    ntet++;
		  else 
		    nother=next->vertex(i)->slice();
		}
	    }
	  if(((nother==slice_index+1)&&(!next->is_heavy_solid_tested_down()))
	     ||((nother==slice_index-1)&&(!next->is_heavy_solid_tested_up())))
	    test_heavy_solidity_around_edge(tr,next,slice_index,heavy_solid,itTet2,itTet3,
					    ntet,nli,nlj,nlk,nll,nother);
	}
      lnext=ll;
      lcommon=lk;
    }
  return heavy_solid;
}

template <class Tr>
void non_solid_connections_removal(Tr & tr, bool heavy_removal)
{
  typename Tr::Finite_cells_iterator it;
  for (it  = tr.finite_cells_begin();
       it != tr.finite_cells_end();
       it++)
    {
      typename Tr::Cell_handle c(it);
      if(c->is_internal())
	determine_solidity(tr,c,heavy_removal);
    }
}

/* template <class Tr> */
/* void non_solid_connections_heavy_removal(Tr & tr) */
/* { */
/*   typename Tr::Finite_cells_iterator it; */
/*   for (it  = tr.finite_cells_begin(); */
/*        it != tr.finite_cells_end(); */
/*        it++) */
/*     { */
/*       typename Tr::Cell_handle c(it); */
/*       if(c->is_internal()) */
/* 	{ */
/* 	  int index_slice=c->vertex(0)->slice(), other_slice_index;	   */
/* 	  int tet1, li1, lj1, tet2, li2, lj2, temp; */
	  
/* 	  tet1=tr.slices[index_slice].tetrahedra_type(c,li1,lj1,other_slice_index); */
/* 	  tet2=tr.slices[other_slice_index].tetrahedra_type(c,li2,lj2,temp); */
/* 	  CGAL_triangulation_assertion((temp==index_slice)&&(tet1+tet2==4)); */
	  
/* 	  if (tet1==2) */
/* 	    determine_heavy_solidity(tr,c,index_slice,li1,lj1,other_slice_index,li2,lj2); */
/* 	} */
/*     } */
/* } */

template <class Tr>
void stable_non_solid_connections_heavy_removal(Tr & tr)
{
  static int cmpt=0;
  bool again;
  do
    {
      again=false;
      typename Tr::Finite_cells_iterator it;
      for (it  = tr.finite_cells_begin();
	   it != tr.finite_cells_end();
	   it++)
	{
	  typename Tr::Cell_handle c(it);
	  if(c->is_internal())
	    {
	      int index_slice=c->vertex(0)->slice(), other_slice_index;	  
	      int tet1, li1, lj1, tet2, li2, lj2, temp;
	      
	      tet1=tr.slices[index_slice].tetrahedra_type(c,li1,lj1,other_slice_index);
	      tet2=tr.slices[other_slice_index].tetrahedra_type(c,li2,lj2,temp);
	      CGAL_triangulation_assertion((temp==index_slice)&&(tet1+tet2==4));
	      
	      if (tet1==2)
		if (!determine_heavy_solidity(tr,c,index_slice,li1,lj1,other_slice_index,li2,lj2))
		  {
		    again=true;
		    CGAL_triangulation_assertion(!c->is_internal());
		  }
	    }
	}
      if(again)
	{
	  for (it  = tr.finite_cells_begin();
	       it != tr.finite_cells_end();
	       it++)
	    {
	      it->unset_heavy_solid_tested_down();
	      it->unset_heavy_solid_tested_up();
	      it->unset_solid_tested_down();
	      it->unset_solid_tested_up();
	    }
	  non_solid_connections_removal(tr,true);
	  //non_solid_connections_removal(tr,true); //TO DO determine which call is faster
	  std::cerr << "Info: Heavy solid removal " << ++cmpt << std::endl;
	}
    }while(again);
}
/*===========================================================================*/
/*         OFF_FILE WRITING                                                  */
/*===========================================================================*/

template <class Tr>
void save_in_off_file(const Tr & tr, const char * fname_off)
{  
  int l_pref=std::strlen(fname_off); 
  char *fname_head=new char[l_pref+10];
  char *fname_body=new char[l_pref+10];

  std::strcat(std::strcpy(fname_head,fname_off),".head");
  std::strcat(std::strcpy(fname_body,fname_off),".body");
  
  // TO DO : make the concatenation of head and body 
  // without system command    
  char *system_command=new char[3*l_pref+40];
  std::strcat(std::strcat(std::strcat(std::strcpy(system_command,"cat "),fname_head)," "),fname_body);
  std::strcat(std::strcat(system_command," > "),fname_off);
  std::ofstream headFile(fname_head, std::ios::out);
  std::ofstream oFile(fname_body, std::ios::out);
  headFile << "OFF\n";
  // write nvertex, nface, nedge(=0)
  headFile << tr.number_of_vertices();
  // write vertex coordinates  
  int n=0;
  typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin ();//Vertex_iterator
  while(vit!=tr.finite_vertices_end())
    {
      vit->set_index(n++);
      oFile << *vit << std::endl;
      ++vit;
    }

  int face_number=0;
  // write oriented faces (counterclockwise from the outside)
  for (typename Tr::Cell_iterator it=tr.cells_begin(); it != tr.cells_end(); ++it)
    {
      typename Tr::Cell_handle c(it);
      if(c->is_external())
	for(int i=0;i<4;i++)
	  if(c->neighbor(i)->is_internal())
	    {  
	      face_number++;   
	      oFile << "3";
	      //for (int j=1;j<4;j++)
		//		oFile << " " << c->vertex((i+j)%4)->index();
	      for (int j=0;j<3;j++)
		oFile <<  " " << c->vertex(Tr::vertex_triple_index(i,j))->index();
	      oFile << std::endl;
	    }
    }
  oFile.close();
  headFile << " " << face_number
	   << " " << 0 << std::endl;
  headFile.close();
  // TO DO : make the concatenation of head and body 
  // without system command  
  std::system(system_command);
  std::strcat(std::strcat(std::strcat(std::strcpy(system_command,"rm "),fname_head)," "),fname_body);
  std::system(system_command);
  delete fname_head;
  delete fname_body;
  delete system_command;
}


template <class Tr>
void save_in_wrl_file(const Tr & tr, const char * fname_wrl, bool outputContours)
{
  typedef typename Tr::Vertex_handle Vertex_handle;
  std::ofstream oFile(fname_wrl,std::ios::out);
  VRML_2_ostream vos(oFile);

  
  std::string s("                Shape {\n"\
		"                    appearance Appearance { material USE Material }\n"\
		"                    geometry IndexedFaceSet {\n"\
		"                        convex FALSE\n"\
		"                        solid  FALSE\n"\
		"                        coord DEF TheCoordinates Coordinate {\n"\
		"                            point [\n");
  vos << s.c_str() ;  
  int n=0;
  typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin ();//Vertex_iterator
  while(vit!=tr.finite_vertices_end())
    {
      vit->set_index(n++);
      vos << to_double(vit->point().x()) << " "  << to_double(vit->point().y()) << " " << 10*to_double(vit->point().z()) << ",\n";
      ++vit;
    }

  s = std::string("] #point\n"
		  "   } #coord Coordinate\n"
		  " coordIndex  [\n");
  vos << s.c_str() ;  
  int face_number=0;
  // write oriented faces (counterclockwise from the outside)
  for (typename Tr::Cell_iterator it=tr.cells_begin(); it != tr.cells_end(); ++it)
    {
      typename Tr::Cell_handle c(it);
      if(c->is_external())
	for(int i=0;i<4;i++)
	  if(c->neighbor(i)->is_internal())
	    {  
	      face_number++;   
	      for (int j=0;j<3;j++){
		vos <<  " " << c->vertex(Tr::vertex_triple_index(i,j))->index() << ", ";
	      }
	      vos << "-1,\n";
	    }
    }

  s = std::string("        ] #coordIndex\n"
		  "      } #geometry\n"
		  "    } #Shape\n");
  vos << s.c_str() ;

  if(outputContours){
    vit = tr.finite_vertices_begin ();//Vertex_iterator
    while(vit!=tr.finite_vertices_end()){
      Vertex_handle vh = vit;
      if(vh->index() != -1){
	
	
	s = std::string("     Shape {\n"
			"       geometry IndexedLineSet {\n"
			"         coord USE TheCoordinates\n"
			"         coordIndex [\n");
	vos << s.c_str() ;
	
	int start = vh->index();
	do {
	  vos << vh->index() << " ";
	  vh->set_index(-1);
	  vh = vh->next();
	}while(vh->index() != -1);
	vos << start ;
	
	s = std::string("         ]\n"
			"       }\n"
			"     } #Shape\n");
	vos << s.c_str() ;
      }
      vit++;
    }
  }
  
}

CGAL_END_NAMESPACE

#endif //CGAL_RECONSTRUCTION_FROM_SLICES_3_H
