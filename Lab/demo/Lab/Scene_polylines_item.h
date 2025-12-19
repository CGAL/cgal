#ifndef SCENE_POLYLINES_ITEM_H
#define SCENE_POLYLINES_ITEM_H
#include "Scene_polylines_item_config.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Three/Scene_group_item.h>

#include <QString>
#include <QMenu>

#include <list>
#include <vector>

struct Scene_polylines_item_private;
class Scene_spheres_item;

class SCENE_POLYLINES_ITEM_EXPORT Scene_polylines_item
    : public CGAL::Three::Scene_group_item
{
    Q_OBJECT
public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef std::vector<Point_3> Polyline;
    typedef std::list<Polyline> Polylines_container;

    Scene_polylines_item();
    virtual ~Scene_polylines_item();
    enum STATS {
      NB_VERTICES = 0,
      NB_EDGES,
      MIN_LENGTH,
      MAX_LENGTH,
      MEAN_LENGTH
    };
    bool has_stats()const override {return true;}
    QString computeStats(int type)override ;
    CGAL::Three::Scene_item::Header_data header() const override ;

    bool isFinite() const override { return true; }
    bool isEmpty() const override ;
    void compute_bbox() const override ;
    Bbox bbox() const override ;

    Scene_polylines_item* clone() const override ;

    QString toolTip() const override ;

    // Indicate if rendering mode is supported
    bool supportsRenderingMode(RenderingMode m) const override ;

    QMenu* contextMenu() override ;

    void draw(CGAL::Three::Viewer_interface*) const override ;
    void drawEdges(CGAL::Three::Viewer_interface*) const override ;
    void drawPoints(CGAL::Three::Viewer_interface*) const override ;


    void smooth(std::vector<Point_3>& polyline);
    //When selecting a polylineitem, we don't want to select its children, so
    //we can still apply Operations to it
    QList<Scene_interface::Item_id> getChildrenForSelection() const override {
      return QList<Scene_interface::Item_id>();
    }
    void setWidth(int i);
    void computeElements() const override;
    void initializeBuffers(Viewer_interface *) const override;

public Q_SLOTS:
    void invalidateOpenGLBuffers() override;
    void change_corner_radii(double);
    void change_corner_radii();
    void split_at_sharp_angles();
    void reset_spheres();

    void merge(Scene_polylines_item*);

    void smooth();
    void point_set_from_polyline();
public:
    Polylines_container polylines;
protected:
    // https://en.wikipedia.org/wiki/D-pointer
    friend struct Scene_polylines_item_private;
    Scene_polylines_item_private* d;

}; // end class Scene_polylines_item

#endif
