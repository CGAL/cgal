// ============================================================================
//
// $Id$
//
// vertex of a combinatotial triangulation of any dimension <=3
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_VERTEX_H
#define CGAL_TRIANGULATION_DS_VERTEX_H

#include <CGAL/Triangulation_short_names.h>
#include <CGAL/Triangulation_ds_iterators.h>

template <class Vb, class Cb >
class  CGAL_Triangulation_ds_vertex 
  : public Vb
{

public:
  typedef typename Vb::Point Point;
  typedef CGAL_Triangulation_ds_vertex<Vb,Cb> Vertex;
  typedef CGAL_Triangulation_ds_cell<Vb,Cb> Cell;
//   typedef CGAL_Triangulation_ds_cell_circulator_2<Vertex,Cell> Cell_circulator;
//   typedef CGAL_Triangulation_ds_vertex_circulator_2<Vertex,Cell> Vertex_circulator;
//   typedef CGAL_Triangulation_ds_edge_circulator_2<Vertex,Cell> Edge_circulator;

  // CONSTRUCTORS

  inline 
  CGAL_Triangulation_ds_vertex()
    : Vb()
  {}
    
  inline 
  CGAL_Triangulation_ds_vertex(const Point & p)
    :  Vb(p)
  {}
    
  inline 
  CGAL_Triangulation_ds_vertex(const Point & p, Cell * c)
    :  Vb(p, c )
  {}

  // ACCESS

  // point()
  //inherited from Vb

  inline 
  Cell* cell() const
  {
    return ( (Cell *) (Vb::cell()) );
  }
    
  // SETTING

  // set_point()
  //inherited from Vb

  inline 
  void set_cell(Cell* c)
  {
    Vb::set_cell(c);
  }

  // USEFUL CONSTANT TIME FUNCTIONS

  // USEFUL NON CONSTANT TIME FUNCTIONS

//   int degree() const
//   {
// TBD with edge circulator
//   }
  
//   inline Vertex_circulator incident_vertices()
//   {
//     return Vertex_circulator(this, cell());
//   }
    
//   inline Cell_circulator incident_cells()
//   {
//     return Cell_circulator(this, cell());
//   }
    
//   inline Cell_circulator incident_cells(const Cell* c)
//   {
//     return Cell_circulator(this, c);
//   }
    
//   inline Edge_circulator incident_edges()
//   {
//     return Edge_circulator(this, cell());
//   }
       
//   inline Edge_circulator incident_edges(const Cell* c)
//   {
//     return Edge_circulator(this, c);
//   }

    
  // CHECKING

   bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Vb::is_valid();
    result = result && cell()->has_vertex(this);
    return result;
  }
};

#endif CGAL_TRIANGULATION_DS_VERTEX_H
