#ifndef WHICH_DIAGRAM_H
#define WHICH_DIAGRAM_H

#include <CGAL/basic.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>

CGAL_BEGIN_NAMESPACE

template<class Matching_class> struct Which_diagram;

template<class Gt, class DS, class LTag>
struct Which_diagram< Segment_Voronoi_diagram_2<Gt,DS,LTag> >
{
  typedef Tag_false Is_hierarchy;
};

template<class Gt, class STag, class DS, class LTag>
struct Which_diagram< Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
{
  typedef Tag_true  Is_hierarchy;
};

CGAL_END_NAMESPACE


#endif // WHICH_DIAGRAM_H
