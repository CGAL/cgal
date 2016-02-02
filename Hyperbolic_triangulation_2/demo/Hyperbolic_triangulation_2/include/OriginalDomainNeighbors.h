
#ifndef CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS
#define CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS

#include <CGAL/circulator.h>

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>

#include <PointGraphicsItem.h>

namespace CGAL {
namespace Qt {

template <typename DT, typename Translation>
class OriginalDomainNeighbors : public GraphicsViewInput
{
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;
  
  OriginalDomainNeighbors(QGraphicsScene& scene_, DT  * dt_, QObject* parent, Point p_ = Point(0, 0), int color = 0);
  
  template <typename InputIterator>
  void assign(InputIterator begin, InputIterator end)
  {
    insert(begin, end);
  }
  
protected:
  template <typename InputIterator>
  void insert(InputIterator begin, InputIterator end);
  
  Point do_point_translation(const Translation& translation) const;
  virtual Vertex_handle insert_point(const Point& p);

private:
  DT * dt;
  
  // center of original domain
  Point p;
  int color;
};


template <typename DT, typename Translation>
OriginalDomainNeighbors<DT, Translation>::OriginalDomainNeighbors(QGraphicsScene& ,
                                                                  DT * dt_,
                                                                  QObject* parent,
                                                                  Point p_,
                                                                  int color_)
  :  GraphicsViewInput(parent), dt(dt_), p(p_), color(color_)
{
}

  
template <typename DT, typename Translation>
template <typename InputIterator>
void
OriginalDomainNeighbors<DT, Translation>::insert(InputIterator begin, InputIterator end)
{
  insert_point(p);
  
  for(InputIterator it = begin; it != end; ++it) {
    Translation t = it->g;
    Point neighbor = do_point_translation(it->g);
      
    insert_point(neighbor);
  }  
}


template <typename DT, typename Translation>
typename OriginalDomainNeighbors<DT, Translation>::Point 
OriginalDomainNeighbors<DT, Translation>::do_point_translation(const Translation& translation) const
{
  return translation.DoAction(p);
}

  
template <typename DT, typename Translation>
typename OriginalDomainNeighbors<DT, Translation>::Vertex_handle 
OriginalDomainNeighbors<DT, Translation>::insert_point(const Point& p)
{
  Vertex_handle v = dt->insert(p);
  v->info().setColor(color);
  
  emit(modelChanged());
  
  return v;
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS
