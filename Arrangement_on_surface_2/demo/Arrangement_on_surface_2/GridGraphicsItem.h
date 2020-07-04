#ifndef ARRANGEMENT_DEMO_GRID_GRAPHICS_ITEM_H
#define ARRANGEMENT_DEMO_GRID_GRAPHICS_ITEM_H

#include <QGraphicsItem>
#include <QColor>

class QWidget;
class QPainter;
class QStyleOptionGraphicsItem;

class GridGraphicsItem : public QGraphicsItem
{
public:
  GridGraphicsItem();
  int getXPower5();
  int getXPower2();
  int getYPower5();
  int getYPower2();
  float getXUnit();
  float getYUnit();

  QColor getGridColor();
  QColor getAxesColor();
  QColor getLabelsColor();
  int getSpacing();

  void setGridColor(const QColor&);
  void setAxesColor(const QColor&);
  void setLabelsColor(const QColor&);
  void setSpacing(int);

protected:
  void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;
  QRectF viewportRect(QPainter*);
  QRectF boundingRect() const override;

private:
  QColor gridColor;
  QColor axesColor;
  QColor labelsColor;

  int spacing;
  int x_power5;
  int x_power2;
  int y_power5;
  int y_power2;
  float x_unit;
  float y_unit;
};

#endif
