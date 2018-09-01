#include "PolygonListItem.h"

PolygonListItem::PolygonListItem(QString name) {
  m_name = name;
  m_color = *new QColor();
  m_visibility = true;
  m_polygonItem = NULL;
}
PolygonListItem::PolygonListItem(QString name, QColor color, bool visibility, 
                                 PolygonItem* polygonItem) {
  m_name = name;
  m_color = color;
  m_visibility = visibility;
  m_polygonItem = polygonItem;
}

QString PolygonListItem::getName() {
  return m_name;
}
QColor PolygonListItem::getColor() { 
  return m_color;
}
bool PolygonListItem::isVisible() {
  return m_visibility;
}
PolygonItem* PolygonListItem::getPolygonItem() {
  return m_polygonItem;
}

void PolygonListItem::setName(QString name) {
  m_name = name;
}
void PolygonListItem::setColor(QColor color) {
  m_color = color; 
}
void PolygonListItem::setVisible(bool visibility) {
  m_visibility = visibility;
}
void PolygonListItem::setPolygonItem(PolygonItem* polygonItem) {
  m_polygonItem = polygonItem;
}
