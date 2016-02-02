
#ifndef CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS_COMMON
#define CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS_COMMON

#include <CGAL/circulator.h>

#include <CGAL/Qt/GraphicsViewInput.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>

#include <PointGraphicsItem.h>

namespace CGAL {
namespace Qt {

template <typename DT, typename Translation>
class OriginalDomainNeighborsCommon : public GraphicsViewInput
{
public:
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Point Point;
  
  OriginalDomainNeighborsCommon(QGraphicsScene& scene_, DT  * dt_, QObject* parent, Point p_, int color);
  
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
  
  // extra - not common
  template<typename Info>
  void add_info(Vertex_handle v, Info new_info);
  void add_color(Vertex_handle v);
  
private:
  DT * dt;
  
  // center of original domain
  Point p;
  int color;
};


template <typename DT, typename Translation>
OriginalDomainNeighborsCommon<DT, Translation>::OriginalDomainNeighborsCommon(QGraphicsScene& ,
                                                                  DT * dt_,
                                                                  QObject* parent,
                                                                  Point p_ = Point(0, 0),
                                                                  int color_ = 0)
  :  GraphicsViewInput(parent), dt(dt_), p(p_), color(color_)
{
}

  
template <typename DT, typename Translation>
template <typename InputIterator>
void
OriginalDomainNeighborsCommon<DT, Translation>::insert(InputIterator begin, InputIterator end)
{
  insert_point(p);
  
  for(InputIterator it = begin; it != end; ++it) {
    Translation t = *it;
    Point neighbor = do_point_translation(*it);
      
    Vertex_handle v = insert_point(neighbor);
    
    add_info(v, it->info());
    add_color(v);
  }  
}


template <typename DT, typename Translation>
typename OriginalDomainNeighborsCommon<DT, Translation>::Point 
OriginalDomainNeighborsCommon<DT, Translation>::do_point_translation(const Translation& translation) const
{
  return translation.DoAction(p);
}

  
template <typename DT, typename Translation>
typename OriginalDomainNeighborsCommon<DT, Translation>::Vertex_handle 
OriginalDomainNeighborsCommon<DT, Translation>::insert_point(const Point& p)
{
  Vertex_handle v = dt->insert(p);  
  emit(modelChanged());

  return v;
}

  
template <typename DT, typename Translation>
template <typename Info>
void 
OriginalDomainNeighborsCommon<DT, Translation>::add_info(Vertex_handle v, Info new_info)
{
  v->info() = new_info;
}


template <typename DT, typename Translation>
void 
OriginalDomainNeighborsCommon<DT, Translation>::add_color(Vertex_handle v)
{
  v->info().setColor(color);
}
  

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ORIGINAL_DOMAIN_NEIGHBORS_COMMON
