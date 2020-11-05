#include <CGAL/Three/TextRenderer.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_print_item_interface.h>
#include "Scene_polyhedron_selection_item.h"
void TextRenderer::draw(CGAL::Three::Viewer_interface *viewer, const QVector3D& scaler)
{
    QPainter *painter = viewer->getPainter();
    if (!painter->isActive())
      painter->begin(viewer);
    QRect rect;
    CGAL::qglviewer::Camera* camera = viewer->camera();
    //Display the items textItems
    Q_FOREACH(TextListItem* list, textItems)
    {
      CGAL::Three::Scene_print_item_interface* item =
      qobject_cast<CGAL::Three::Scene_print_item_interface*>(scene->item(scene->mainSelectionIndex()));
      if( item &&
          item->shouldDisplayIds(list->item())
         )
        Q_FOREACH(TextItem* item, list->textList())
        {
          CGAL::qglviewer::Vec src(item->position().x(), item->position().y(),item->position().z());
          if(viewer->testDisplayId(src.x, src.y, src.z))
          {
            if(item->is_3D())
            {
              src.x *= scaler.x();
              src.y *= scaler.y();
              src.z *= scaler.z();
              rect = QRect(int(camera->projectedCoordinatesOf(src).x -item->width()/2),
                           int(camera->projectedCoordinatesOf(src).y -item->height()/2),
                           int(item->width()),
                           int(item->height()));
            }
            else
              rect = QRect(int(src.x-item->width()/2),
                           int(src.y-item->height()/2),
                           int(item->width()),
                           int(item->height()));

            painter->setFont(item->font());
            QColor c = item->color().toHsv();
            c.setHsv((c.hsvHue()+180)%360, 255,255,100);
            painter->setBrush(QBrush(c));
            painter->setPen(QPen(QColor(0,0,0,0)));
            painter->drawRect(rect);
            painter->setPen(QPen(item->color()));
            painter->drawText(rect, item->text());
          }
        }
    }

    //Display the local TextItems
    Q_FOREACH(TextItem* item, local_textItems)
    {
      CGAL::qglviewer::Vec src(item->position().x(), item->position().y(),item->position().z());
      if(item->is_3D())
      {
        if(item->is_always_visible() || viewer->testDisplayId(src.x, src.y, src.z))
        {
            rect = QRect(int(camera->projectedCoordinatesOf(src).x-item->width()/2),
                         int(camera->projectedCoordinatesOf(src).y-item->height()/2),
                         int(item->width()),
                         int(item->height()));
        }
      }
      else
      {
          rect = QRect(int(src.x-item->width()/2),
                       int(src.y-item->height()/2),
                       int(item->width()),
                       int(item->height()));
      }
      painter->setFont(item->font());
      QColor c = item->color().toHsv();
      c.setHsv((c.hsvHue()+180)%360, 255,255,100);
      painter->setBrush(QBrush(c));
      painter->setPen(QPen(QColor(0,0,0,0)));
      painter->drawRect(rect);
      painter->setPen(QPen(item->color()));
      painter->drawText(rect, item->text());
    }
}

 void TextRenderer::addTextList(TextListItem *tl)
 {
   if(tl->textList().size() > max_textItems)
   {
     Q_EMIT sendMessage("There are too many textItems to display.",5000);
     return;
   }
     textItems.append(tl);
 }

 void TextRenderer::addText(TextItem *ti)
 {
     local_textItems.append(ti);
 }

 void TextRenderer::addText(float p_x, float p_y, float p_z, QString p_text, bool p_3D, QFont p_font , QColor p_color )
 {
     local_textItems.append(new TextItem(p_x, p_y, p_z, p_text, p_3D, p_font, p_color));
 }

 void TextRenderer::removeText(TextItem *item)
 {
     local_textItems.removeAll(item);
 }

 void TextRenderer::removeTextList(TextListItem *p_list)
 {
     Q_FOREACH(TextListItem *list, textItems)
         if(list == p_list)
             textItems.removeAll(list);
 }
