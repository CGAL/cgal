#ifndef CGAL_QT_GRAPHICS_VIEW_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_INPUT_H

#include <CGAL/Object.h>
#include <QObject>

namespace CGAL {
namespace Qt {
class GraphicsViewInput  : public QObject
{
  Q_OBJECT

public:
  GraphicsViewInput(QObject* parent) 
    : QObject(parent)
  {}

signals:
  void generate(CGAL::Object o);
  void modelChanged();
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_INPUT_H
