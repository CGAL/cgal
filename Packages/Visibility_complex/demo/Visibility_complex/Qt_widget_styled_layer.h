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
<<<<<<< Qt_widget_styled_layer.h

  typedef QMap<QString,QVariant> Map;
  typedef Map::iterator iterator;
=======

  typedef QMap<QString,QVariant> Map;
>>>>>>> 1.2
public:

  typedef Map::const_iterator const_iterator;
  typedef Map::size_type size_type;

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

  const_iterator begin() const
  {
    return map.begin();
  }

  const_iterator end() const
  {
    return map.end();
  }

  size_type size() const
  {
    return map.size();
  }

private:
  Map map;
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
