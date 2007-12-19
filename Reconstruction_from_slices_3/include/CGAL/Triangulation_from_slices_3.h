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
// Author(s)     : Raphaelle Chaine (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//                 Andreas Fabri (GeometryFactory)
//   
//===============================================================================================

/*! \file Triangulation_from_slices_3.h
    \brief A triangulation class that
*/         
#ifndef CGAL_TRIANGULATION_FROM_SLICES_3_H
#define CGAL_TRIANGULATION_FROM_SLICES_3_H


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// Equivalent to Filtered_kernel<Cartesian<Lazy_exact_nt<quotient<mp_float > > > >

#include <iostream>
#include <vector>
#include <utility> 

#include <CGAL/utility.h>
#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/TFS_vertex_base_3.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Parser_SLC.h>

#ifdef GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#endif // GEOMVIEW

#define CGAL_RFS_CONFORM //if no access to mesh_2

CGAL_BEGIN_NAMESPACE


std::pair<int,int> get_opposite_index_pair(int i,int j);//i et j in O..3

/*!
 * In this triangulation the points are in layers and cells only have
 * vertices in neighboring layers.
 * 
*/

template < class Gt = Exact_predicates_exact_constructions_kernel, 
           class Tds = Triangulation_data_structure_3 < TFS_vertex_base_3 <Triangulation_vertex_base_3<Gt> >,
                                                        Triangulation_cell_base_3<Gt> > >
class Triangulation_from_slices_3 : public Triangulation_3<Gt,Tds>
  {  
    
    /*==================================*/
    /*         Type Definitions         */
    /*==================================*/
    
    typedef Triangulation_from_slices_3<Gt,Tds> Self;
    typedef Triangulation_3<Gt,Tds> Tr_Base;
    
    public:
    typedef Tds Triangulation_data_structure;
    typedef Gt  Geom_traits;
    
    typedef typename Gt::Point_3       Point;
    typedef typename Gt::Segment_3     Segment;
    typedef typename Gt::Triangle_3    Triangle;
    typedef typename Gt::Tetrahedron_3 Tetrahedron;
    
    typedef typename Gt::Vector_3 Vector;
    
    typedef typename Gt::Plane_3 Plane;
    
    // to do : types for dual

    typedef typename Tr_Base::Cell   Cell;
    typedef typename Tr_Base::Vertex Vertex;
    typedef typename Tr_Base::Facet  Facet;
    typedef typename Tr_Base::Edge   Edge;

    typedef typename Tr_Base::Cell_handle   Cell_handle;
    typedef typename Tr_Base::Vertex_handle Vertex_handle;

    typedef typename Tr_Base::Cell_circulator Cell_circulator;
    typedef typename Tr_Base::Cell_iterator   Cell_iterator;
    typedef typename Tr_Base::Facet_iterator  Facet_iterator;
    typedef typename Tr_Base::Edge_iterator   Edge_iterator;
    typedef typename Tr_Base::Vertex_iterator Vertex_iterator;

    typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename Tr_Base::Finite_cells_iterator    Finite_cells_iterator;
    typedef typename Tr_Base::Finite_facets_iterator   Finite_facets_iterator;
    typedef typename Tr_Base::Finite_edges_iterator    Finite_edges_iterator;
    typedef typename Tr_Base::Facet_circulator         Facet_circulator;

    typedef typename Tr_Base::All_cells_iterator       All_cells_iterator;

    typedef typename Tr_Base::Locate_type Locate_type;

    typedef Triple<Vertex_handle,Vertex_handle,Vertex_handle> Vertex_triple; 

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
    using Tr_Base::cw;
    using Tr_Base::ccw;
    using Tr_Base::geom_traits;
    using Tr_Base::number_of_vertices;
    using Tr_Base::dimension;
    using Tr_Base::finite_facets_begin;
    using Tr_Base::finite_facets_end;
    using Tr_Base::finite_vertices_begin;
    using Tr_Base::finite_vertices_end;
    using Tr_Base::finite_cells_begin;
    using Tr_Base::finite_cells_end;
    using Tr_Base::finite_edges_begin;
    using Tr_Base::finite_edges_end;
    using Tr_Base::tds;
    using Tr_Base::infinite_vertex;
    using Tr_Base::next_around_edge;
    using Tr_Base::vertex_triple_index;
#endif

    /*==================================*/
    /*      Nested Class        */
    /*==================================*/
    
    class Slice;

    class Conflict_tester_3;
    class Conflict_tester_2;

    class Perturbation_order; // Copy paste from Delaunay_Triangulation_3

    
    friend class Slice;
    friend class Perturbation_order;
    friend class Conflict_tester_3;
    friend class Conflict_tester_2;  

    Edge get_opposite_edge_in_facet(Facet f,int i) const; //i in O..2
     
    private :
    /*==================================*/
    /*          Member Datas            */
    /*==================================*/
    public : //TO RM and replace by private (after access to slices for reconstruction from slices)

    std::vector<Slice> slices;
    unsigned int represented_slices_number; // May be suppressed ?
    //int nb_slices;

      /*==================================*/
      /*     Public Members Functions     */
      /*==================================*/
  
    public:
      //CONSTRUCTOR     
    Triangulation_from_slices_3(const Gt& gt = Gt())
    : Tr_Base(gt), represented_slices_number(0)
    {}
  
    // copy constructor duplicates vertices and cells
    Triangulation_from_slices_3(const Triangulation_from_slices_3 & tr)
    : Tr_Base(tr), represented_slices_number(represented_slices_number)
    { 
      CGAL_triangulation_postcondition( is_valid() );
    }

    template < typename InputIterator >
    Triangulation_from_slices_3(InputIterator first, InputIterator last,
					const Gt& gt = Gt())
    : Tr_Base(gt), represented_slices_number(0)
    {
      insert(first, last);
      CGAL_triangulation_postcondition( is_valid() );
    }

    Triangulation_from_slices_3(const char* slc_file_name,
					const Gt& gt = Gt())
    : Tr_Base(gt), represented_slices_number(0)
    {      
      parser_SLC<Self> slcp(slc_file_name,this);
      slcp.parse();
      CGAL_triangulation_postcondition( is_valid() );
    }

    void clear()
    {
      Tr_Base::clear();
      slices.clear();
    }


    template <class Parser>
    void load(const char* slc_file_name);
    //Precondition : the triangulation must contain no points */

    template <class OutputIteratorBoundaryFacets,
              class OutputIteratorCells,
              class OutputIteratorInternalFacets>
    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           OutputIteratorInternalFacets>
    find_conflicts(const Point &p,int slice, Cell_handle c,
		   OutputIteratorBoundaryFacets bfit,
		   OutputIteratorCells cit,
		   OutputIteratorInternalFacets ifit) const;

    // FIND CONFLICT
    template <class Conflict_test,
              class OutputIteratorBoundaryFacets,
              class OutputIteratorCells,
              class OutputIteratorInternalFacets>
    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           OutputIteratorInternalFacets>
    economic_find_conflicts_3(Cell_handle c, const Conflict_test &tester,
			      Triple<OutputIteratorBoundaryFacets,
			      OutputIteratorCells,
			      OutputIteratorInternalFacets> it) const;
                             //IN PROGRESS : DO NOT USE NOW

    template <class Conflict_test,
              class OutputIteratorBoundaryFacets,
              class OutputIteratorCells,
              class OutputIteratorInternalFacets>
    Triple<OutputIteratorBoundaryFacets,
           OutputIteratorCells,
           OutputIteratorInternalFacets>
    find_conflicts_3(Cell_handle c, const Conflict_test &tester,
	             Triple<OutputIteratorBoundaryFacets,
		     OutputIteratorCells,
		     OutputIteratorInternalFacets> it) const;
    
    /*   Vertex_handle */ //Offered in Delaunay_triangulation_3
    /*   nearest_vertex_in_cell(const Point& p, const Cell_handle& c) const; */

    
    /*     private: */ //Offered in Delaunay_triangulation_3
    /*     Vertex_handle */
    /*     nearest_vertex(const Point &p, Vertex_handle v, Vertex_handle w) const; */


/**
 * Starts a new slice in the triangulation
 * @param slice is the index of the slice
 * @param plane is the plane of the slice.
 * Currently the slices must be inserted in ascending order and all points
 * of a slice must be inserted before the next slice is inserted.
 */
    void insert_slice(int slice, const Plane & plane);

    template < class InputIterator >
    // Iterator on std::pair<Point,int>
    int
    insert(InputIterator first, InputIterator last);
/**
 * Inserts a point in a specific slice of the triangulation
 * @param p the point to insert
 * @param slice is where p gets inserted.
 * Currently the slices must be inserted in ascending order and all points
 * of a slice must be inserted before the next slice is inserted.
 */
    Vertex_handle insert(const Point & p, int slice, Cell_handle cstart = NULL);

    Vertex_handle insert(const std::pair<Point,int> & pt_s, Cell_handle cstart = NULL);
    Vertex_handle insert(const Point & p, int slice, Locate_type lt,
			 Cell_handle c, int li, int);
    Vertex_handle insert(const std::pair<Point,int> & pt_s, Locate_type lt,
			 Cell_handle c, int li, int);

    // VALIDITY
    bool
    is_valid(bool verbose = false, int level = 0) const;
    
    bool 
    is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

    //--------------------------------------------------------------------
    protected:
  
    Oriented_side //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    side_of_oriented_sphere(const Point &p0, const Point &p1, const Point &p2,
			    const Point &p3, const Point &t, bool perturb = false) const;


    Bounded_side //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    coplanar_side_of_bounded_circle(const Point &p, const Point &q,
				    const Point &r, const Point &s, bool perturb = false) const;
 
    bool //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    less_distance(const Point &p, const Point &q, const Point &r) const
    {
      return geom_traits().compare_distance_3_object()(p, q, r) == SMALLER;
    }

    public :
/*   void */
/*   make_canonical(Vertex_triple& t) const; */
  
/*   Vertex_triple */
/*   make_vertex_triple(const Facet& f) const; */
    private :
    Bounded_side//COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    side_of_sphere(const Vertex_handle& v0, const Vertex_handle& v1,
		   const Vertex_handle& v2, const Vertex_handle& v3,
		   const Point &p, bool perturb) const;
    public:
    //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    Bounded_side
    side_of_sphere(const Cell_handle& c, const Point & p,
		   bool perturb = false) const
    {
      return side_of_sphere(c->vertex(0), c->vertex(1),
                            c->vertex(2), c->vertex(3), p, perturb);
    }
    //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    Bounded_side
    side_of_circle( const Facet & f, const Point & p, bool perturb = false) const
    {
      return side_of_circle(f.first, f.second, p, perturb);
    }
    //COPY_PASTE from Delaunay_triangulation_3 (How to access it without heritage ...)
    Bounded_side
    side_of_circle( const Cell_handle& c, int i, const Point & p,
		    bool perturb = false) const;
    
    //-------------------------------------------------------------------

  }; //Triangulation_from_slices_3<Gt,Tds> 

template < class Gt, class Tds >
class Triangulation_from_slices_3<Gt,Tds>::Slice
{ 
 private :      
  int slice_index; // May be suppressed ?
  Plane slice_equation; 
  int vertex_number;
  const Self *t;  
  Vertex_handle vert_start; // To fasten locate operations (Easier to maintain than a Cell_handle)

#ifdef CGAL_RFS_CONFORM
 public : 
  std::vector<Vertex_handle> set_of_vertices; //Is just here to fasten conform operations
                                            //TO DO This container could be a template parameter 
                                            //      that could be set to EmptySet by default
#endif //CGAL_RFS_CONFORM

  // Construction------------------------------------------------------------------------------
 public :    
  Slice(const Self *tr=NULL, int num=-1) 
    : slice_index(num), vertex_number(0), t(tr), vert_start(NULL){}	
    void set_start(Vertex_handle vh) 
    {
      CGAL_triangulation_precondition(vh->slice()==slice_index);
      vert_start=vh;
    }       
    void set_index(int num) {slice_index=num;}
    int get_index() const {return slice_index;}
    void set_plane(const Plane & p) {slice_equation=p;}
    const Plane & get_plane() const {return slice_equation;}
    int get_vertex_number() const {return vertex_number;}
    void set_triangulation(const Self *tr) {t=tr;}
    Vertex_handle get_vertex_start() const {return vert_start;}
    // Locate a point --------------------------------------------------------------------------
    Cell_handle locate_in_slice(const Point & p, Locate_type & lt, int & li, int & lj,  
				Cell_handle s=NULL) const;
    // Returns a finite cell whose restriction to *this slice contains p,
    // or an infinite cell
    // Starts at cell "start"
    // if lt == OUTSIDE_CONVEX_HULL, li is the
    // index of a facet separating p from the rest of the triangulation
    // in dimension 2 :
    // returns a facet (Cell_handle,li) if lt == FACET
    // returns an edge (Cell_handle,li,lj) if lt == EDGE
    // returns a vertex (Cell_handle,li) if lt == VERTEX
    // if lt == OUTSIDE_CONVEX_HULL, li, lj give the edge of c
    // separating p from the rest of the triangulation
    
    // Type of cells with regards to a slice -----------------------------------------------------
    int tetrahedra_type(Cell_handle c, int & li, int & lj, int & index_other_slice) const 
    // Returns the number of c vertices on *this slice
    // case 1 : li index of the vertex that lies on *this slice
    // case 2 : li and lj indexes of vertices belonging to *this slice
    // case 3 : li index of the vertex that does not belong to *this slice
    // index_other_slice represents the index of the other slice where c vertices lie on
    // Postcondition : return value between 0 and 3
    {
      li=lj=index_other_slice=-1;
      int indexes[3]; int opposite=-1;
      int cmpt=0;
      for (int i=0;i<4;i++)
	{
	  if (c->vertex(i)->slice()==slice_index)
	    indexes[cmpt++]=i;
	  else
	    {
	      opposite=i;
	      if (c->vertex(i)->slice()!=-1) //infinite vertex
		index_other_slice=c->vertex(i)->slice();
	    }
	}
      CGAL_triangulation_assertion(cmpt<4);
      li=lj=-1;
      switch(cmpt)
	{
	case 0: break;
	case 1:
	  li=indexes[0]; break;	    
	case 2:
	  li=indexes[0]; lj=indexes[1]; break;
	case 3:
	  li=opposite; break;	    
	}
      return cmpt;
    }
    
    int tetrahedra_type(Cell_handle c) const    
    { int li,lj,other_slice;
      return tetrahedra_type(c,li,lj,other_slice);
    }
    
    bool test_tetrahedra_type(Cell_handle c,int nbvertex,int li=-1,int lj=-1, int other_slice_index=-1) const
    {
      int i,j,other_slice;
      if (nbvertex==tetrahedra_type(c,i,j,other_slice))
	{
	  switch (nbvertex)
	    {
	    case 3:
	    case 1:
	      return ((li==-1)
		      ||((li==i)&&((other_slice_index==-1)
				   ||(other_slice_index==other_slice))));
	    case 2:
	      return ((li==-1)
		      ||((li==i)&&(lj==j)&&((other_slice_index==-1)
					    ||(other_slice_index==other_slice))));
	    case 0:
	      return true;
	    }
	}
      return false;
    } 
    
    // Restriction of the 3D Triangulation to *this slice ------------------------------------------------
    Facet mirror_facet(Facet f) const
    {
      CGAL_triangulation_precondition(test_tetrahedra_type(f.first,3,f.second));
      // c is a T1next or T2prev of *this slice
      return std::make_pair(f.first->neighbor(f.second),f.first->mirror_index(f.second));
    }
    
    Facet neighbor_facet(Facet f, int i) const
    {
      CGAL_triangulation_precondition(test_tetrahedra_type(f.first,3,f.second));	
      // c is a T1next or T2prev of *this slice
      Vertex_handle vopp=f.first->vertex(f.second);
      int ii=vertex_triple_index(f.second,i); //Index of the ith vertex  of the facet opposite to vertex f.second
      
      if(vopp->slice==-1) //ie. if (is_infinite(vopp))
	{
	  CGAL_triangulation_assertion(is_infinite(vopp)); //TO RM 
	  return std::make_pair(f.first->neighbor(ii),
				f.first->neighbor(ii)->index(vopp));
	}
/* 	else */
/* 	  { */
/* 	    Facet mf=mirror_facet(f); */
/* 	    if(mf.first->vertex(mf.second)->slice==-1) */
/* 	      return mirror_facet( TO DO ); */
/* 	  } */
      else
	{
	  int jj=vertex_triple_index(f.second,(i+1)%3);
	  int kk=vertex_triple_index(f.second,(i+2)%3);
	  Facet_circulator fc=incident_facets(f,jj,kk),finit(fc); //TO DO : verify the initialization of fc
	  CGAL_triangulation_assertion(*fc == f); //TO RM 
	  int oppslice_init=(*finit).first->vertex((*finit).second)->slice;
	  int oppslice=oppslice_init;
	  do
	    {
	      fc++;
	      oppslice=(*fc).first->vertex((*fc).second)->slice;
	    }
	  while(oppslice==oppslice_init);
	  fc++;
	  Facet fres((*fc).first->neighbor((*fc).second),(*fc).first->mirror_index((*fc).second));
	  CGAL_triangulation_assertion(((fres.first->vertex(fres.second)->slice==oppslice_init)
					&&(test_tetrahedra_type(fres.first,3,fres.second) ))
				       ||((test_tetrahedra_type(fres.first,2))
					  &&is_infinite(fres.first)));
	  // returned tetrahedra_type 3 or(2 and infinite) //TO RM 
	  return fres;
	}	
    }
    bool has_full_dimension()
    {
      std::vector<Cell_handle> temp;
      t->incident_cells(this->vert_start,std::back_inserter(temp));
      for (typename std::vector<Cell_handle>::iterator cit = temp.begin();
	   cit != temp.end(); ++cit) 
	{
	  if (this->tetrahedra_type(*cit)==3)
	    return true;
	}
      return false;
    }
    friend class Triangulation_from_slices_3<Gt,Tds>;
}; //class  Triangulation_from_slices_3<Gt,Tds>::slice   

template < class Gt, class Tds >
class Triangulation_from_slices_3<Gt,Tds>::Conflict_tester_3
{
 public : //TO RM
  const Point &p;
  int slice_index;
  const Self *t;
  
 public:
  
  Conflict_tester_3(const Point &pt, int num, const Self *tr)
    : p(pt), slice_index(num), t(tr) 
    {
      CGAL_triangulation_precondition((t->slices[num].get_index()==num) 
				      //ie. the slice has previously been constructed
				      &&(t->slices[num].get_plane().has_on(pt)));
    }
    
    bool operator()(const Cell_handle c) const
    {	    
      if (t->slices[slice_index].get_vertex_number()==0)
	if (t->is_infinite(c)&&((t->slices[slice_index-1].tetrahedra_type(c)==3)
				||(t->slices[slice_index+1].tetrahedra_type(c)==3))) //TO COMMENT
	  return true;
	else
	  return false;
      else //(t->slices[slice_index].get_vertex_number()!=0)	      
	switch(t->slices[slice_index].tetrahedra_type(c))
	  {
	  case 3 : //test 2D 
	    // return t->side_of_circle(c, ..., p, true) == ON_BOUNDED_SIDE;
	  case 1 : //test 3D
	  case 2 : //test 3D
	    return t->side_of_sphere(c,p,true) == ON_BOUNDED_SIDE;
	  default :
	    return false;
	  }
    }
    
    bool operator()(const Cell_handle c,int tet, int li, int lj, int index_other_slice) const
    {
      CGAL_triangulation_precondition(t->slices[slice_index].tetrahedra_type(c,tet,li,lj,index_other_slice));
      switch(tet)
	{
	case 3 : //test 2D 
	  // return t->side_of_circle(c, ..., p, true) == ON_BOUNDED_SIDE;
	case 1 :
	case 2 : //test 3D
	  return t->side_of_sphere(c,p,true) == ON_BOUNDED_SIDE;
	default :
	  if (t->slices[slice_index].get_vertex_number()!=0)
	    return false;
	  else if (t->is_infinite(c)&&((t->slices[slice_index-1].tetrahedra_type(c)==3)
				       ||(t->slices[slice_index+1].tetrahedra_type(c)==3))) //TO COMMENT
	    return true;
	  else
	    return false;
	}
    }
}; //class Triangulation_from_slices_3<Gt,Tds>::Conflict_tester_3

template < class Gt, class Tds >
class Triangulation_from_slices_3<Gt,Tds>::Conflict_tester_2
{
  const Point &p;      
  int slice_index;
  const Self *t;
  
 public:

  Conflict_tester_2(const Point &pt,int num, const Self *tr)
    : p(pt), slice_index(num), t(tr) 
    {
      CGAL_triangulation_precondition((t->slices[num].get_index()==num)
				      &&(t->slices[num].get_plane().has_on(pt)));
    }
    
    bool operator()(const Cell_handle c) const
    {
      return t->side_of_circle(c, 3, p, true) == ON_BOUNDED_SIDE;
    }
};//class Triangulation_from_slices_3<Gt,Tds>::Conflict_tester_2

template < class Gt, class Tds >
class Triangulation_from_slices_3<Gt,Tds>::Perturbation_order
{
  const Self *t;
  
 public:
  Perturbation_order(const Self *tr)
    : t(tr) {}
    
    bool operator()(const Point *p, const Point *q) const 
    {
      return t->compare_xyz(*p, *q) == SMALLER;
    }
};//class Triangulation_from_slices_3<Gt,Tds>::Perturbation_order


std::pair<int,int> get_opposite_index_pair(int i,int j)
{  
  CGAL_triangulation_precondition((i>=0)&&(i<4));
  CGAL_triangulation_precondition((j>=0)&&(j<4));
  CGAL_triangulation_precondition(i!=j);
 
  int dk, dl;
  for (dk=1;dk<4;dk++)
    if((i+dk)%4!=j)
      break;
  for (dl=dk+1;dl<4;dl++)
    if((i+dl)%4!=j)
      break;
  return std::make_pair((i+dk)%4,(i+dl)%4);  
}

template < class Gt, class Tds >
template <class Parser>
void Triangulation_from_slices_3<Gt,Tds>::load(const char* slc_file_name)
{  
  Parser slcp(slc_file_name,this);
  slcp.parse();
}

template < class Gt, class Tds >
template <class OutputIteratorBoundaryFacets,
          class OutputIteratorCells,
          class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets,
       OutputIteratorCells,
       OutputIteratorInternalFacets>
  Triangulation_from_slices_3<Gt,Tds>::find_conflicts(const Point &p,int slice, Cell_handle c,
							    OutputIteratorBoundaryFacets bfit,
							    OutputIteratorCells cit,
							    OutputIteratorInternalFacets ifit) const
{
  CGAL_triangulation_precondition(dimension() >= 2);
  
  std::vector<Cell_handle> cells;
  cells.reserve(32);
  std::vector<Facet> facets;
  facets.reserve(64);
  
  if (dimension() == 2) 
    {// call find_conflicts_2 (Triangulation_3.h)
      Conflict_tester_2 tester(p,slice,this);
      ifit = find_conflicts_2(c, tester,
			      make_triple(std::back_inserter(facets),
					  std::back_inserter(cells),
					  ifit)).third;
    }
  else 
    {
      Conflict_tester_3 tester(p,slice,this);
      //call redefined find_conflicts_3
      ifit = find_conflicts_3(c, tester,
			      make_triple(std::back_inserter(facets),
					  std::back_inserter(cells),
					  ifit)).third;
    }	    	 	 
  
  // Reset the conflict flag on the boundary.
  for(typename std::vector<Facet>::iterator fit=facets.begin();
      fit != facets.end(); ++fit) 
    {
      fit->first->neighbor(fit->second)->set_in_conflict_flag(0);
      *bfit++ = *fit;
    }
  
  // Reset the conflict flag in the conflict cells.
  for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
      ccit != cells.end(); ++ccit) 
    {
      (*ccit)->set_in_conflict_flag(0);
      *cit++ = *ccit;
    }
  return make_triple(bfit, cit, ifit);
}

template < class Gt, class Tds >
template <class Conflict_test,
          class OutputIteratorBoundaryFacets,
          class OutputIteratorCells,
          class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  Triangulation_from_slices_3<Gt,Tds>::find_conflicts_3(Cell_handle c, const Conflict_test &tester,
							      Triple<OutputIteratorBoundaryFacets,
							      OutputIteratorCells,
							      OutputIteratorInternalFacets> it) const
{
  CGAL_triangulation_precondition( dimension()==3 );
  CGAL_triangulation_precondition( tester(c) );
  
  it = Tr_Base::find_conflicts_3(c,tester,it);
  //it=economic_find_conflicts_3(c,tester,it);
  
  return it;      
}


template < class Gt, class Tds >
void Triangulation_from_slices_3<Gt,Tds>::insert_slice(int slice_num,const Plane & equ)
{  
  slices.push_back(Slice(this,slice_num));
  slices[slice_num].set_plane(equ);   
}

template < class Gt, class Tds >
template < class InputIterator >
// Iterator on std::pair<Point,int>
int 
Triangulation_from_slices_3<Gt,Tds>::insert(InputIterator first, InputIterator last)
{
  int n = number_of_vertices();
  while(first != last)
    {
      insert(*first);
      ++first;
    }
  return number_of_vertices() - n;
}

template < class Gt, class Tds >
  typename Triangulation_from_slices_3<Gt,Tds>::Edge
  Triangulation_from_slices_3<Gt,Tds>::get_opposite_edge_in_facet(typename Triangulation_from_slices_3<Gt,Tds>::Facet f,int j) const
{
  CGAL_triangulation_precondition((j>=0)&&(j<3));
  return Edge(f.first,this->vertex_triple_index(f.second,(j+1)%3),this->vertex_triple_index(f.second,(j+2)%3));  
}

#ifdef GEOMVIEW
template < class Gt, class Tds> 
  CGAL::Geomview_stream & 
  operator<< ( CGAL::Geomview_stream &gv, const Triangulation_from_slices_3<Gt,Tds> & tr);
#endif //GEOMVIEW

template < class Gt, class Tds >
typename Triangulation_from_slices_3<Gt,Tds>::Vertex_handle
Triangulation_from_slices_3<Gt,Tds>::
insert(const Point & p, int slice, Cell_handle start)
{
  CGAL_triangulation_precondition(static_cast<unsigned int>(slice)<=represented_slices_number);
  CGAL_triangulation_precondition(((dimension()<2)&&(slice==0))
				  ||((dimension()==2)&&((slice==0)||(slice==1)))
				  ||((dimension()==3)&&((static_cast<unsigned int>(slice)<represented_slices_number)
							||((static_cast<unsigned int>(slice)==represented_slices_number)
							   &&(slices.size()>represented_slices_number)
							   &&(slices[represented_slices_number-1].has_full_dimension())
							   ))));
  Locate_type lt;
  int li, lj;
  if(dimension()>2)
    if(start==NULL)
      if(slices[slice].vert_start!=NULL)
	start=slices[slice].vert_start->cell();
      else if(slices[slice-1].vert_start!=NULL)
	start=slices[slice-1].vert_start->cell();
  Cell_handle c = slices[slice].locate_in_slice(p,lt, li, lj, start);
  return insert(p, slice, lt, c, li, lj);
}

template < class Gt, class Tds >
typename Triangulation_from_slices_3<Gt,Tds>::Vertex_handle
Triangulation_from_slices_3<Gt,Tds>::
insert(const Point & p, int slice, Locate_type lt, Cell_handle c, int li, int)
{ 
  CGAL_triangulation_precondition(static_cast<unsigned int>(slice)<=represented_slices_number);
  CGAL_triangulation_precondition(((dimension()<2)&&(slice==0))
				  ||((dimension()==2)&&((slice==0)||(slice==1)))
				  ||((dimension()==3)&&((static_cast<unsigned int>(slice)<represented_slices_number)
							||((static_cast<unsigned int>(slice)==represented_slices_number)
							   &&(slices.size()>represented_slices_number)
							   &&(slices[represented_slices_number-1].has_full_dimension())))));
  Vertex_handle v;
  switch (dimension()) 
    {
    case 3:
      {
	if ( lt == Tr_Base::VERTEX )
	  return c->vertex(li);
	// TEST call insert_conflicts_3 (Triangulation_3.h) with redefined find_conflicts_3
	Conflict_tester_3 tester(p,slice,this);
	v = insert_conflict(c, tester);//v = insert_conflict_3(c, tester);
#ifdef CGAL_RFS_CONFORM
	slices[slice].set_of_vertices.push_back(v);
#endif
	v->set_point(p);
	break;// dim 3 
      }
    case 2:
      switch (lt) 
	{
	case Tr_Base::OUTSIDE_CONVEX_HULL:
	case Tr_Base::CELL:
	case Tr_Base::FACET:
	case Tr_Base::EDGE:
	  {
	    Conflict_tester_2 tester(p,slice,this);
	    v = insert_conflict(c, tester);//v = insert_conflict_2(c, tester);
#ifdef CGAL_RFS_CONFORM
	    slices[slice].set_of_vertices.push_back(v);
#endif
	    v->set_point(p); 
	    break;
	  }
	case Tr_Base::VERTEX:
	  return c->vertex(li);
	case Tr_Base::OUTSIDE_AFFINE_HULL:
	  // if the 2d triangulation is Delaunay, the 3d
	  // triangulation will be Delaunay 
	  v=Tr_Base::insert_outside_affine_hull(p);
	  break;
	}
      break;
      //dim 2
    default :
      // dimension <= 1
      v=Tr_Base::insert(p, c);
      break;
    }  
  v->set_slice(slice);
  if(slices[slice].vertex_number==0)
    represented_slices_number++;
  slices[slice].vertex_number++;
  slices[slice].set_start(v);
  return v;
}

template < class Gt, class Tds >
  bool
  Triangulation_from_slices_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if ( ! tds().is_valid(verbose,level) ) 
    {
      if (verbose)
	std::cerr << "invalid data structure" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }
  
  if ( infinite_vertex() == NULL ) 
    {
      if (verbose)
	std::cerr << "no infinite vertex" << std::endl;
      CGAL_triangulation_assertion(false);
      return false;
    }

  switch ( dimension() ) 
    {
    case 3:
      {
	Finite_cells_iterator it;
	for ( it = finite_cells_begin(); it != finite_cells_end(); ++it ) 
	  {
	    is_valid_finite(it);
	    int li,lj,other_slice,num,ty,w=-1;
	    num=it->vertex(0)->slice();
	    ty=slices[num].tetrahedra_type(*it,li,lj,other_slice);
	    CGAL_triangulation_assertion(ty!=0);
	    if((ty==1)||(ty==3))
	      w=li;	    
	    for(int i=0; i<4; i++ ) 
	      {
		if((i!=w)&&(!is_infinite(it->neighbor(i)->vertex(it->neighbor(i)->index(it))))) 
		  {
		    if(side_of_sphere(it, 
				      it->neighbor(i)->vertex(it->neighbor(i)->index(it))->point())
		       == ON_BOUNDED_SIDE ) 
		      {//Non locally Delaunay
			if (verbose)
			  std::cerr << "non-empty sphere " << std::endl;
			CGAL_triangulation_assertion(false);
			return false;
		      }
		  }
	      }
	  }
	break;
      }
  case 2:
    {
      Finite_facets_iterator it;
      for ( it = finite_facets_begin(); it != finite_facets_end(); ++it ) 
	{
	  is_valid_finite((*it).first);
	  for (int i=0; i<3; i++ ) 
	    {
	      if( !is_infinite
		  ((*it).first->neighbor(i)->vertex( (((*it).first)->neighbor(i))
						     ->index((*it).first))) ) 
		{
		  if ( side_of_circle ( (*it).first, 3,
					(*it).first->neighbor(i)->
					vertex( (((*it).first)->neighbor(i))
						->index((*it).first) )->point() )
		       == ON_BOUNDED_SIDE ) 
		    {
		      if (verbose)
			std::cerr << "non-empty circle " << std::endl;
		      CGAL_triangulation_assertion(false);
		      return false;
		    }
		}
	    }
	}
      break;
    }
  case 1:
    {
      Finite_edges_iterator it;
      for ( it = finite_edges_begin(); it != finite_edges_end(); ++it )
	is_valid_finite((*it).first);
      break;
    }
  }
  if (verbose)
    std::cerr << "Delaunay valid triangulation" << std::endl;
  return true;
}

template < class Gt, class Tds >
  bool
  Triangulation_from_slices_3<Gt,Tds>::
is_valid(Cell_handle c, bool verbose, int level) const
{
  if ( ! c->is_valid(dimension(),verbose,level) ) 
    {
      if (verbose) { 
	std::cerr << "combinatorically invalid cell" ;
	for (int i=0; i <= dimension(); i++ )
	  std::cerr << c->vertex(i)->point() << ", " ;
	std::cerr << std::endl;
      }
      CGAL_triangulation_assertion(false);
      return false;
    }
  switch ( dimension() ) 
    {    
    case 3:
      {
	if ( ! is_infinite(c) ) 
	  {
	    is_valid_finite(c,verbose,level);
	    int li,lj,other_slice,num,ty,w=-1;
	    num=c->vertex(0)->slice;
	    ty=slices[num]->tetrahedra_type(c,li,lj,other_slice);
	    CGAL_triangulation_assertion(ty!=0);
	    if((ty==1)||(ty==3))
	      w=li;	    	    
	    for (int i=0; i<4; i++ ) 
	      {
		if ((i!=w)&&(side_of_sphere(c, c->vertex((c->neighbor(i))->index(c))->point()))
		    == ON_BOUNDED_SIDE )
		  {
		    if (verbose)
		      std::cerr << "non-empty sphere " << std::endl;
		    CGAL_triangulation_assertion(false);
		    return false;
		  }
	      }
	  }
	break;
      }
    case 2:
      {
	if ( ! is_infinite(c,3) ) 
	  {
	    for (int i=0; i<2; i++ ) 
	      {
		if (side_of_circle(c, 3, c->vertex(c->neighbor(i)->index(c))->point())
		    == ON_BOUNDED_SIDE ) 
		  {
		    if (verbose)
		      std::cerr << "non-empty circle " << std::endl;
		    CGAL_triangulation_assertion(false);
		    return false;
		  }
	      }
	  }
	break;
      }
    }
  if (verbose)
    std::cerr << "Delaunay valid cell" << std::endl;
  return true;
}

template < class Gt, class Tds >
typename Triangulation_from_slices_3<Gt,Tds>::Cell_handle Triangulation_from_slices_3<Gt,Tds>::Slice::locate_in_slice(const Point & p, Locate_type & lt, int & li, int & lj, Cell_handle s) const      
    {
      CGAL_triangulation_precondition(slice_equation.has_on(p));
      Cell_handle cres;
      if((t->dimension()>2)&&(vertex_number==0))
	{
	  std::vector<Cell_handle> temp;
	  t->incident_cells(t->infinite_vertex(),std::back_inserter(temp));
	  for (typename std::vector<Cell_handle>::iterator cit = temp.begin();cit != temp.end(); ++cit) 
	    {
	      int temp=-1;
	      if (t->slices[slice_index-1].tetrahedra_type(*cit,li,lj,temp)==3)
		{
		  cres=*cit;
		  lt=Tr_Base::OUTSIDE_CONVEX_HULL;
		  //li=cres->index(t->infinite_vertex()); 
#ifdef CGAL_DUMP
		  std::cout << "Tr_Base::OUTSIDE_CONVEX_HULL (modifie)" << std::endl;	
#endif	    
		  break;
		}
	    }
	}
      else
	{
	  if ((s==NULL)&&(vert_start!=NULL))
	    s=vert_start->cell();
	  cres=t->locate(p, lt, li, lj, s);
	  if((t->dimension()==3)&&(lt==Tr_Base::OUTSIDE_CONVEX_HULL)&&(tetrahedra_type(cres)==0))
	    {
	      std::cerr << "Error: Case that should not happen thanks to vert_start but to be handled by the author " << std::endl;
	      //TO DO : Visit infinite cells having a vertex in *this to find a good candidate		
	    }
#ifdef CGAL_DUMP
	  if (lt==Tr_Base::OUTSIDE_CONVEX_HULL)
	    std::cout << "Tr_Base::OUTSIDE_CONVEX_HULL" << std::endl;	
	  else if (lt==Tr_Base::VERTEX)
	    std::cout << "Tr_Base::VERTEX" << std::endl;
	  else if (lt==Tr_Base::EDGE)	
	    std::cout << "Tr_Base::EDGE" << std::endl;
	  else if (lt==Tr_Base::FACET)
	    std::cout << "Tr_Base::FACET" << std::endl;
	  else if (lt==Tr_Base::OUTSIDE_AFFINE_HULL)
	    std::cout << "Tr_Base::OUTSIDE_AFFINE_HULL" << std::endl;
	  else
	    std::cout << "Non satisfied precondition on slices" << std::endl;
#endif
	}
      CGAL_triangulation_assertion((t->dimension()<3)
				   ||(vertex_number==0)
				   ||((lt==Tr_Base::VERTEX)&&(cres->vertex(li)->slice()==slice_index))
				   ||((lt==Tr_Base::OUTSIDE_CONVEX_HULL)&&(tetrahedra_type(cres)>0))
				   ||((lt==Tr_Base::EDGE)&&(tetrahedra_type(cres)>1)
				      &&(cres->vertex(li)->slice()==slice_index)
				      &&(cres->vertex(lj)->slice()==slice_index))
				   ||((lt==Tr_Base::FACET)&&(tetrahedra_type(cres)>1)
				      &&((cres->vertex(li)->slice()==slice_index+1)
					 ||(cres->vertex(li)->slice()==slice_index-1)
					 ||(cres->vertex(li)->slice()==-1)/*infinite vertex*/)));
      return cres;
    } 

/* template < class Gt, class Tds > */
/* void */
/* Delaunay_triangulation_3<Gt,Tds>:: */
/* make_canonical(Vertex_triple& t) const */
/* { */
/*   int i = (&*(t.first) < &*(t.second))? 0 : 1; */
/*   if(i==0) { */
/*     i = (&*(t.first) < &*(t.third))? 0 : 2; */
/*   } else { */
/*     i = (&*(t.second) < &*(t.third))? 1 : 2; */
/*   } */
/*   Vertex_handle tmp;  */
/*   switch(i){ */
/*   case 0: return; */
/*   case 1: */
/*     tmp = t.first; */
/*     t.first = t.second; */
/*     t.second = t.third; */
/*     t.third = tmp; */
/*     return; */
/*   default: */
/*     tmp = t.first; */
/*     t.first = t.third; */
/*     t.third = t.second; */
/*     t.second = tmp; */
/*   } */
/* } */




/* template < class Gt, class Tds > */
/* typename Delaunay_triangulation_3<Gt,Tds>::Vertex_triple */
/* Delaunay_triangulation_3<Gt,Tds>:: */
/* make_vertex_triple(const Facet& f) const */
/* { */
/*   static const int vertex_triple_index[4][3] = { {1, 3, 2}, {0, 2, 3}, */
/*                                                  {0, 3, 1}, {0, 1, 2} }; */
/*   Cell_handle ch = f.first; */
/*   int i = f.second; */
  
/*   return Vertex_triple(ch->vertex(vertex_triple_index[i][0]), */
/* 		       ch->vertex(vertex_triple_index[i][1]), */
/* 		       ch->vertex(vertex_triple_index[i][2]));  */
/* } */


//------------------------------------------------------------------------------
// Copy paste from Delaunay_Triangulation_3 and replace Delaunay_triangulation_3 
// with Triangulation_from_slices_3

template < class Gt, class Tds >
  Bounded_side
  Triangulation_from_slices_3<Gt,Tds>::
coplanar_side_of_bounded_circle(const Point &p0, const Point &p1,
				const Point &p2, const Point &p, bool perturb) const
{
  // In dim==2, we should even be able to assert orient == POSITIVE.
  CGAL_triangulation_precondition( coplanar_orientation(p0, p1, p2)
				   != COLLINEAR );

  Bounded_side bs =
    geom_traits().coplanar_side_of_bounded_circle_3_object()(p0, p1, p2, p);

  if (bs != ON_BOUNDARY || !perturb)
    return bs;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[4] = {&p0, &p1, &p2, &p};
  std::sort(points, points+4, Perturbation_order(this) );

  Orientation local = coplanar_orientation(p0, p1, p2);

  // we successively look whether the leading monomial, then 2nd monimial,
  // then 3rd monomial, of the determinant which has non null coefficient
  // [syl] : TODO : Probably it can be stopped earlier like the 3D version
  for (int i=3; i>0; --i) {
    if (points[i] == &p)
      return Bounded_side(NEGATIVE); // since p0 p1 p2 are non collinear
                                     // but not necessarily positively oriented
    Orientation o;
    if (points[i] == &p2
	&& (o = coplanar_orientation(p0,p1,p)) != COLLINEAR )
      // [syl] : TODO : I'm not sure of the signs here (nor the rest :)
      return Bounded_side(o*local);
    if (points[i] == &p1
	&& (o = coplanar_orientation(p0,p,p2)) != COLLINEAR )
      return Bounded_side(o*local);
    if (points[i] == &p0
	&& (o = coplanar_orientation(p,p1,p2)) != COLLINEAR )
      return Bounded_side(o*local);
  }

  // case when the first non null coefficient is the coefficient of 
  // the 4th monomial
  // moreover, the tests (points[] == &p) were false up to here, so the
  // monomial corresponding to p is the only monomial with non-zero
  // coefficient, it is equal to coplanar_orient(p0,p1,p2) == positive 
  // so, no further test is required
  return Bounded_side(-local); //ON_UNBOUNDED_SIDE;
}

template < class Gt, class Tds >
  Bounded_side
  Triangulation_from_slices_3<Gt,Tds>::
side_of_circle(const Cell_handle& c, int i,
	       const Point & p, bool perturb) const
// precondition : dimension >=2
// in dimension 3, - for a finite facet
// returns ON_BOUNDARY if the point lies on the circle,
// ON_UNBOUNDED_SIDE when exterior, ON_BOUNDED_SIDE
// interior
// for an infinite facet, considers the plane defined by the
// adjacent finite facet of the same cell, and does the same as in 
// dimension 2 in this plane
// in dimension 2, for an infinite facet
// in this case, returns ON_BOUNDARY if the point lies on the 
// finite edge (endpoints included) 
// ON_BOUNDED_SIDE for a point in the open half-plane
// ON_UNBOUNDED_SIDE elsewhere
{
  CGAL_triangulation_precondition( dimension() >= 2 );
  int i3 = 5;

  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    // the triangulation is supposed to be valid, ie the facet
    // with vertices 0 1 2 in this order is positively oriented
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
      return coplanar_side_of_bounded_circle( c->vertex(0)->point(),
					      c->vertex(1)->point(),
					      c->vertex(2)->point(),
					      p, perturb);
    // else infinite facet
    // v1, v2 finite vertices of the facet such that v1,v2,infinite
    // is positively oriented
    Vertex_handle v1 = c->vertex( ccw(i3) ),
      v2 = c->vertex( cw(i3) );
    CGAL_triangulation_assertion(coplanar_orientation(v1->point(), v2->point(),
						      (c->mirror_vertex(i3))->point()) == NEGATIVE);
    Orientation o = coplanar_orientation(v1->point(), v2->point(), p);
    if ( o != COLLINEAR )
      return Bounded_side( o );
    // because p is in f iff
    // it does not lie on the same side of v1v2 as vn
    int i_e;
    Locate_type lt;
    // case when p collinear with v1v2
    return side_of_segment( p,
			    v1->point(), v2->point(),
			    lt, i_e );
  }

  // else dimension == 3
  CGAL_triangulation_precondition( i >= 0 && i < 4 );
  if ( ( ! c->has_vertex(infinite_vertex(),i3) ) || ( i3 != i ) ) {
    // finite facet
    // initialization of i0 i1 i2, vertices of the facet positively 
    // oriented (if the triangulation is valid)
    int i0 = (i>0) ? 0 : 1;
    int i1 = (i>1) ? 1 : 2;
    int i2 = (i>2) ? 2 : 3;
    CGAL_triangulation_precondition( coplanar( c->vertex(i0)->point(),
				               c->vertex(i1)->point(),
				               c->vertex(i2)->point(),
					       p ) );
    return coplanar_side_of_bounded_circle( c->vertex(i0)->point(),
					    c->vertex(i1)->point(),
					    c->vertex(i2)->point(),
					    p, perturb);
  }

  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle v1 = c->vertex( next_around_edge(i3,i) ),
    v2 = c->vertex( next_around_edge(i,i3) );
  Orientation o = (Orientation)
    (coplanar_orientation( v1->point(), v2->point(),
			   c->vertex(i)->point()) *
     coplanar_orientation( v1->point(), v2->point(), p ));
  // then the code is duplicated from 2d case
  if ( o != COLLINEAR )
    return Bounded_side( -o );
  // because p is in f iff 
  // it is not on the same side of v1v2 as c->vertex(i)
  int i_e;
  Locate_type lt;
  // case when p collinear with v1v2
  return side_of_segment( p,
			  v1->point(), v2->point(),
			  lt, i_e );
}

template < class Gt, class Tds >
  Bounded_side
  Triangulation_from_slices_3<Gt,Tds>::
side_of_sphere(const Vertex_handle& v0, const Vertex_handle& v1,
	       const Vertex_handle& v2, const Vertex_handle& v3,
	       const Point &p, bool perturb) const
{
  CGAL_triangulation_precondition( dimension() == 3 );

  // TODO :
  // - avoid accessing points of infinite vertex
  // - share the 4 codes below (see old version)
  const Point &p0 = v0->point();
  const Point &p1 = v1->point();
  const Point &p2 = v2->point();
  const Point &p3 = v3->point();

  if (is_infinite(v0)) {
    Orientation o = orientation(p2, p1, p3, p);
    if (o != COPLANAR)
      return Bounded_side(o);
    return coplanar_side_of_bounded_circle(p2, p1, p3, p, perturb);
  }

  if (is_infinite(v1)) {
    Orientation o = orientation(p2, p3, p0, p);
    if (o != COPLANAR)
      return Bounded_side(o);
    return coplanar_side_of_bounded_circle(p2, p3, p0, p, perturb);
  }

  if (is_infinite(v2)) {
    Orientation o = orientation(p1, p0, p3, p);
    if (o != COPLANAR)
      return Bounded_side(o);
    return coplanar_side_of_bounded_circle(p1, p0, p3, p, perturb);
  }

  if (is_infinite(v3)) {
    Orientation o = orientation(p0, p1, p2, p);
    if (o != COPLANAR)
      return Bounded_side(o);
    return coplanar_side_of_bounded_circle(p0, p1, p2, p, perturb);
  }

  return (Bounded_side) side_of_oriented_sphere(p0, p1, p2, p3, p, perturb);
}

template < class Gt, class Tds >
  Oriented_side
  Triangulation_from_slices_3<Gt,Tds>::
side_of_oriented_sphere(const Point &p0, const Point &p1, const Point &p2,
	                const Point &p3, const Point &p, bool perturb) const
{
  CGAL_triangulation_precondition( orientation(p0, p1, p2, p3) == POSITIVE );
  
  Oriented_side os =
    geom_traits().side_of_oriented_sphere_3_object()(p0, p1, p2, p3, p);

  if (os != ON_ORIENTED_BOUNDARY || !perturb)
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[5] = {&p0, &p1, &p2, &p3, &p};
  std::sort(points, points+5, Perturbation_order(this) );

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for (int i=4; i>2; --i) {
      if (points[i] == &p)
	return ON_NEGATIVE_SIDE; // since p0 p1 p2 p3 are non coplanar
                                 // and positively oriented
      Orientation o;
      if (points[i] == &p3 && (o = orientation(p0,p1,p2,p)) != COPLANAR )
	return Oriented_side(o);
      if (points[i] == &p2 && (o = orientation(p0,p1,p,p3)) != COPLANAR )
	return Oriented_side(o);
      if (points[i] == &p1 && (o = orientation(p0,p,p2,p3)) != COPLANAR )
	return Oriented_side(o);
      if (points[i] == &p0 && (o = orientation(p,p1,p2,p3)) != COPLANAR )
	return Oriented_side(o);
  }

  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE; 
}

// End Copy paste
//------------------------------------------------------------------------------

/**
   Operator << overload : usefull to flush a triangulation in a geomview stream
   ??? Should be located in Triangulation_3 or in Triangulation_data_structure_3
*/
#ifdef GEOMVIEW
template < class Gt, class Tds> 
CGAL::Geomview_stream & 
operator<<( CGAL::Geomview_stream &gv, const Triangulation_from_slices_3<Gt,Tds> & tr)
{  
  typedef Triangulation_from_slices_3<Gt,Tds> Triangulation;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Triangulation::Segment Segment;

  gv.set_edge_color (CGAL::Color(0,0,0));
  //gv.set_face_color (CGAL::Color(255,255,255));
  gv.set_bg_color (CGAL::Color(180,180,180));

  for (Finite_edges_iterator it=tr.finite_edges_begin(); it != tr.finite_edges_end(); ++it)
    gv << Segment((*it).first->vertex((*it).second)->point(),(*it).first->vertex((*it).third)->point());
  gv.look_recenter();
  
  return gv;
}
#endif //GEOMVIEW

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_FROM_SLICES_3_H
