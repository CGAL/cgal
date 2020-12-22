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

  bool isEmpty()const Q_DECL_OVERRIDE;
  bool isFinite() const Q_DECL_OVERRIDE { return true; }

  Scene_lcc_item* clone() const Q_DECL_OVERRIDE ;
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE ;
  QString toolTip() const Q_DECL_OVERRIDE ;
  void compute_bbox()const Q_DECL_OVERRIDE;


  void draw(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE ;
  void drawEdges(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE ;
  void drawPoints(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE ;

  void computeElements() const Q_DECL_OVERRIDE;
  void initializeBuffers(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;

  QMenu* contextMenu() Q_DECL_OVERRIDE ;

public Q_SLOTS:
    void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;
    void randomFaceColors();
    void randomVolumeColors();
    void resetColors();
private:
  friend struct lcc_priv;
  lcc_priv* d;
};

#endif // SCENE_LCC_ITEM_H

