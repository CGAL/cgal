#ifndef SCENE_LCC_ITEM_H
#define SCENE_LCC_ITEM_H

#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include <QMenu>

#ifdef scene_lcc_item_EXPORTS
#  define SCENE_LCC_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_LCC_ITEM_EXPORT Q_DECL_IMPORT
#endif

struct lcc_priv;
namespace CGAL{
namespace Three{
class Viewer_interface;
}
}

class SCENE_LCC_ITEM_EXPORT Scene_lcc_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public:
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
  Scene_lcc_item(const LCC &lcc);
  ~Scene_lcc_item();

  bool isEmpty()const override;
  bool isFinite() const override { return true; }

  Scene_lcc_item* clone() const override ;
  bool supportsRenderingMode(RenderingMode m) const override ;
  QString toolTip() const override ;
  void compute_bbox()const override;


  void draw(CGAL::Three::Viewer_interface*) const override ;
  void drawEdges(CGAL::Three::Viewer_interface*) const override ;
  void drawPoints(CGAL::Three::Viewer_interface*) const override ;

  void computeElements() const override;
  void initializeBuffers(CGAL::Three::Viewer_interface *) const override;

  QMenu* contextMenu() override ;

public Q_SLOTS:
    void invalidateOpenGLBuffers() override;
    void randomFaceColors();
    void randomVolumeColors();
    void resetColors();
private:
  friend struct lcc_priv;
  lcc_priv* d;
};

#endif // SCENE_LCC_ITEM_H

