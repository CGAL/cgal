#ifndef CGAL_POLYGONITEM_H
#define CGAL_POLYGONITEM_H

#include "Typedefs.h"

class PolygonItem
{
public:
  PolygonItem(PolygonWithHoles* polygon) {
    m_polygon = polygon;
    m_graphics = new PolygonGraphicsItem(m_polygon);
  }
  ~PolygonItem() { };

  PolygonWithHoles* polygon() { return m_polygon; };
  PolygonGraphicsItem* graphics() { return m_graphics; };

private:
  PolygonWithHoles* m_polygon;
  PolygonGraphicsItem* m_graphics;
};

#endif // CGAL_POLYGOITEM_H


