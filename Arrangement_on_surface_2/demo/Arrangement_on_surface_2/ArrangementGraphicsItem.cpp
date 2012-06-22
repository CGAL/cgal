#include "ArrangementGraphicsItem.h"

namespace CGAL {
namespace Qt {

ArrangementGraphicsItemBase::
ArrangementGraphicsItemBase( ):
    bb( 0, 0, 0, 0 ),
    bb_initialized( false ),
    visible_edges( true ),
    visible_vertices( true ),
    verticesPen( QPen( ::Qt::black, 3. ) ),
    edgesPen( QPen( ::Qt::black, 1. ) )
{ }

const QPen&
ArrangementGraphicsItemBase::
getVerticesPen() const
{
    return this->verticesPen;
}

const QPen& 
ArrangementGraphicsItemBase::
getEdgesPen() const
{
    return this->edgesPen;
}

void 
ArrangementGraphicsItemBase::
setVerticesPen(const QPen& pen)
{
    this->verticesPen = pen;
}

void 
ArrangementGraphicsItemBase::
setEdgesPen(const QPen& pen)
{
    this->edgesPen = pen;
}

bool 
ArrangementGraphicsItemBase::
visibleVertices() const
{
    return this->visible_vertices;
}

void 
ArrangementGraphicsItemBase::
setVisibleVertices(const bool b)
{
    this->visible_vertices = b;
    this->update();
}

bool 
ArrangementGraphicsItemBase::
visibleEdges() const
{
    return this->visible_edges;
}

void 
ArrangementGraphicsItemBase::
setVisibleEdges(const bool b)
{
    this->visible_edges = b;
    this->update();
}

} // namespace Qt
} // namespace CGAL
