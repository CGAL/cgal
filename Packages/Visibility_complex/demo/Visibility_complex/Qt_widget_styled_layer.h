#ifndef QT_WIDGET_STYLED_LAYER_H
#define QT_WIDGET_STYLED_LAYER_H

#include <qvariant.h>
#include <qstring.h>
#include <qmap.h>

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/function_objects.h>

namespace CGAL {

class Qt_widget_style : public QObject {
  Q_OBJECT
public:
  Qt_widget_style()
    : map() {};

public slots:
  void setBool(QString name, bool b);
  void setInt(QString name, int i);
  void setColor(QString name, QColor c);
  void setPointStyle(QString name, PointStyle s);

public:
  bool getBool(QString name);
  int getInt(QString name);
  QColor getColor(QString name);
  PointStyle getPointStyle(QString name);

private:
  QMap<QString,QVariant> map;
};

class Qt_widget_styled_layer : public Qt_widget_layer {
  Q_OBJECT
public:
  typedef Qt_widget_style Style;

  Qt_widget_styled_layer(Style* style = 0,
			 QObject * parent=0, const char * name=0);

  ~Qt_widget_styled_layer();

  void setStyle(Style* style);
  Style * style()  { return style_; }
private:
  Style* style_;
  bool private_style;
};

} // namespace CGAL

#endif // QT_WIDGET_STYLED_LAYER_H
