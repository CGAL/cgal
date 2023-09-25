#include "Scene_aff_transformed_item.h"

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>

#include <QApplication>

using namespace CGAL::Three;

Scene_aff_transformed_item::
Scene_aff_transformed_item(const CGAL::qglviewer::Vec& pos)
{
  d = new Scene_aff_transformed_item_priv(pos);
}

Scene_aff_transformed_item::
~Scene_aff_transformed_item()
{
  delete d;
}

void
Scene_aff_transformed_item::
itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);
  if(d && item == this)
  {
    if(d->frame)
    {
      delete d->frame;
      d->frame = nullptr;
    }
  }
}

void
Scene_aff_transformed_item::
setFMatrix(double matrix[16])
{
  for(int i=0; i<16; ++i)
    d->f_matrix.data()[i] = float(matrix[i]);
}

bool
Scene_aff_transformed_item::
keyPressEvent(QKeyEvent* e)
{
  if(e->key() == Qt::Key_S)
  {
    Q_EMIT applyTransformation();
    return true;
  }
  return false;
}
