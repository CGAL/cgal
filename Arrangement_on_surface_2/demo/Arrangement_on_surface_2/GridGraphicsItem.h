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
  int getXPower5() const;
  int getXPower2() const;
  int getYPower5() const;
  int getYPower2() const;
  float getXUnit() const;
  float getYUnit() const;

  QColor getGridColor() const;
  QColor getAxesColor() const;
  QColor getLabelsColor() const;
  int getSpacing() const;

  void setGridColor(const QColor&);
  void setAxesColor(const QColor&);
  void setLabelsColor(const QColor&);
  void setSpacing(int);

protected:
  void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;
  QRectF viewportRect(QPainter*) const;
  QRectF viewportRect() const;
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
