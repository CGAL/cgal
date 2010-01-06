#include "Scene_segmented_image_item.h"
#include "Image_type.h"

Scene_segmented_image_item::Bbox
Scene_segmented_image_item::bbox() const
{
  if(!m_image) return Bbox();
  return Bbox(0, 0, 0,
              m_image->xdim() * m_image->vx(),
              m_image->ydim() * m_image->vy(),
              m_image->zdim() * m_image->vz());
}

void Scene_segmented_image_item::draw() const
{
  if(m_image) {
    m_image->gl_draw_bbox(3.0f,0,0,0);
  }
}

QString Scene_segmented_image_item::toolTip() const
{
  return tr("<p>Image <b>%1</b></p>"
            "<p>Word type: %2</p>"
            "<p>Dimensions: %3 x %4 x %5</p>"
            "<p>Spacings: ( %6 , %7 , %8 )</p>")
    .arg(this->name())
    .arg("...")
    .arg(m_image->xdim()) 
    .arg(m_image->ydim())
    .arg(m_image->zdim())
    .arg(m_image->vx())
    .arg(m_image->vy())
    .arg(m_image->vz());
}

#include "Scene_segmented_image_item.moc"

