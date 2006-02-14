#ifndef QT_WIDGET_STYLE_EDITOR_H
#define QT_WIDGET_STYLE_EDITOR_H

#include <qframe.h>
#include "Qt_widget_styled_layer.h"

namespace CGAL {

class Qt_widget_style_editor : public QFrame {
  Q_OBJECT
public:

  typedef Qt_widget_styled_layer::Style Style;

  Qt_widget_style_editor(Style* style,
			 QWidget *parent = 0 , const char *name = 0);

  virtual ~Qt_widget_style_editor() {}

signals:
  void styleChanged();

private slots:
  void map(QColor);
  void map(int);
  void map(bool);

private:
  Style* style;
  QMap<const QObject*, QString> mapper;
}; // end of class Qt_widget_style_editor

} // end namespace CGAL

#endif // QT_WIDGET_STYLE_EDITOR_H
