#ifndef SCENE_SEGMENTED_IMAGE_ITEM_H
#define SCENE_SEGMENTED_IMAGE_ITEM_H

#include "Scene_item.h"
#include "Scene_interface.h"
#include "Image_type_fwd.h"
#include "Scene_segmented_image_item_config.h"

typedef CGAL::Image_3 Image;

class SCENE_SEGMENTED_IMAGE_ITEM_EXPORT Scene_segmented_image_item 
  : public Scene_item
{
  Q_OBJECT
public:

  Scene_segmented_image_item(Image* im)
    : m_image(im)
  {
  }

  ~Scene_segmented_image_item() {}

  bool isFinite() const { return true; }
  bool isEmpty() const { return false; }
  Bbox bbox() const;

  Scene_segmented_image_item* clone() const { return NULL; }

  bool supportsRenderingMode(RenderingMode) const { return false; }

  void draw() const;
  QString toolTip() const;

  const Image* image() const {
    return m_image;
  }
public:
  Image* m_image;
};

#endif // SCENE_SEGMENTED_IMAGE_ITEM_H
