#ifndef CEP_SUPPORT_FUNCTIONS_H
#define CEP_SUPPORT_FUNCTIONS_H

// some simple support functions ...

#include <CGAL/Origin.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

inline CGAL::Orientation reverse_orientation(CGAL::Orientation ori)
{
  switch(ori){
    case CGAL::LEFTTURN:  { return CGAL::RIGHTTURN; }
    case CGAL::RIGHTTURN: { return LEFTTURN; }
    case CGAL::COLLINEAR: { return CGAL::COLLINEAR; }
    //case CLOCKWISE: return COUNTERCLOCKWISE; (like RIGHTTURN)
    //case COUNTERCLOCKWISE: return CLOCKWISE; (like LEFTTURN)
  }

  return ori;
}

CGAL_END_NAMESPACE

#endif
