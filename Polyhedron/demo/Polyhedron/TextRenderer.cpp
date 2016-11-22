#include <CGAL/Three/TextRenderer.h>
#include <CGAL/Three/Scene_item.h>
void TextRenderer::draw(CGAL::Three::Viewer_interface *viewer)
{
    QPainter *painter = viewer->getPainter();
    if (!painter->isActive())
      painter->begin(viewer);
    QRect rect;
    qglviewer::Camera* camera = viewer->camera();
    //Display the items textItems
    Q_FOREACH(TextListItem* list, textItems)
      if(list->item() == scene->item(scene->mainSelectionIndex()))
        Q_FOREACH(TextItem* item, list->textList())
        {
          qglviewer::Vec src(item->position().x(), item->position().y(),item->position().z());
          if(viewer->testDisplayId(src.x, src.y, src.z))
          {
            if(item->is_3D())
              rect = QRect(camera->projectedCoordinatesOf(src).x-item->width()/2,
                           camera->projectedCoordinatesOf(src).y-item->height()/2,
                           item->width(),
                           item->height());
            else
              rect = QRect(src.x-item->width()/2,
                           src.y-item->height()/2,
                           item->width(),
                           item->height());

            painter->setFont(item->font());
            painter->setPen(QPen(item->color()));
            painter->drawText(rect, item->text());
          }
        }

    //Display the local TextItems
    Q_FOREACH(TextItem* item, local_textItems)
    {
      qglviewer::Vec src(item->position().x(), item->position().y(),item->position().z());
      if(item->is_3D())
      {
        if(item->is_always_visible() || viewer->testDisplayId(src.x, src.y, src.z))
        {
            rect = QRect(camera->projectedCoordinatesOf(src).x-item->width()/2,
                       camera->projectedCoordinatesOf(src).y-item->height()/2,
                       item->width(),
                       item->height());
        }
      }
      else
      {
          rect = QRect(src.x-item->width()/2,
                     src.y-item->height()/2,
                     item->width(),
                     item->height());
      }
      painter->setFont(item->font());
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
