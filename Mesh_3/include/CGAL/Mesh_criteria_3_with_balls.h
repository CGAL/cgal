//
// C++ Interface: Mesh_criteria_3_with_balls
//
// Description: 
//
//
// Author:  <>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef CGAL_MESH_CRITERIA_3_WITH_BALLS_H
#define CGAL_MESH_CRITERIA_3_WITH_BALLS_H

#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>
#include <CGAL/Mesh_edge_criteria_3.h>
#include <CGAL/Mesh_3/Facet_criteria_visitor_with_balls.h>
#include <CGAL/Mesh_3/Cell_criteria_visitor_with_balls.h>
#include <CGAL/Mesh_3/Facet_on_same_surface_criterion.h>

namespace CGAL {

// Class Mesh_criteria_3
// Provides default meshing criteria to drive Mesh_3 process
template <class Tr>
class Mesh_criteria_3_with_balls
{
public:
  typedef Mesh_facet_criteria_3<Tr,Mesh_3::Facet_criterion_visitor_with_balls<Tr>  >  Facet_criteria;
  typedef Mesh_cell_criteria_3<Tr,Mesh_3::Cell_criteria_visitor_with_balls<Tr>  >     Cell_criteria;
  typedef Mesh_edge_criteria_3<Tr>                                                    Edge_criteria;

  // Constructor (backwards compatibility with Polyhedron_demo)
  // TODO: remove it !
  Mesh_criteria_3_with_balls(const Facet_criteria& facet_criteria,
                             const Cell_criteria& cell_criteria,
                             bool activate_same_surface_criteria = true)
    : edge_criteria_(0)
    , facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria)
  {
    typedef Mesh_3::Facet_criterion_visitor_with_balls<Tr> Facet_visitor;
    typedef Mesh_3::Facet_on_same_surface_criterion<Tr,Facet_visitor>
      Same_surface_criterion;
    
    if(activate_same_surface_criteria) {
      facet_criteria_.add(new Same_surface_criterion());
    }
  }
  
  /// Constructor
  Mesh_criteria_3_with_balls(const Edge_criteria& edge_criteria,
                             const Facet_criteria& facet_criteria,
                             const Cell_criteria& cell_criteria,
                             bool activate_same_surface_criteria = true)
    : edge_criteria_(edge_criteria)
    , facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria)
  {
    typedef Mesh_3::Facet_criterion_visitor_with_balls<Tr> Facet_visitor;
    typedef Mesh_3::Facet_on_same_surface_criterion<Tr,Facet_visitor>
    Same_surface_criterion;
    
    if(activate_same_surface_criteria) {
      facet_criteria_.add(new Same_surface_criterion());
    }
  }

  // Destructor
  ~Mesh_criteria_3_with_balls() { };

  const Edge_criteria& edge_criteria() const { return edge_criteria_; }
  const Facet_criteria& facet_criteria() const { return facet_criteria_; };
  const Cell_criteria& cell_criteria() const { return cell_criteria_; };
  
private:
  Edge_criteria edge_criteria_;
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;

};  // end class Mesh_criteria_3


}  // end namespace CGAL


#endif // CGAL_MESH_CRITERIA_3_WITH_BALLS_H
