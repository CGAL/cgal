#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include <CGAL/Three/Scene_print_item_interface.h>
#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/TextRenderer.h>
#include "Polyhedron_type_fwd.h"
#include "Polyhedron_type.h"
#include <iostream>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <set>
#include <vector>

#include <QColor>
#include <CGAL/Three/Scene_zoomable_item_interface.h>

class QMenu;
struct Scene_polyhedron_item_priv;

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_ITEM_EXPORT Scene_polyhedron_item
        : public CGAL::Three::Scene_item,
          public CGAL::Three::Scene_zoomable_item_interface,
          public CGAL::Three::Scene_print_item_interface{
    Q_INTERFACES(CGAL::Three::Scene_print_item_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PrintInterface/1.0")
    Q_OBJECT
    Q_INTERFACES(CGAL::Three::Scene_zoomable_item_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.ZoomInterface/1.0")
public:
    typedef Polyhedron Face_graph;
    typedef boost::property_map<Face_graph, boost::vertex_index_t>::type Vertex_selection_map;
    typedef boost::property_map<Face_graph, boost::face_index_t>::type Face_selection_map;

    enum STATS {
      NB_VERTICES = 0,
      NB_CONNECTED_COMPOS,
      NB_BORDER_EDGES,
      IS_PURE_TRIANGLE,
      NB_DEGENERATED_FACES,
      HOLES,
      AREA,
      VOLUME,
      SELFINTER,
      NB_FACETS,
      MIN_AREA,
      MAX_AREA,
      MED_AREA,
      MEAN_AREA,
      MIN_ALTITUDE,
      MIN_ASPECT_RATIO,
      MAX_ASPECT_RATIO,
      MEAN_ASPECT_RATIO,
      GENUS,
      NB_EDGES,
      MIN_LENGTH,
      MAX_LENGTH,
      MID_LENGTH,
      MEAN_LENGTH,
      NB_NULL_LENGTH,
      MIN_ANGLE,
      MAX_ANGLE,
      MEAN_ANGLE
    };

    bool has_stats()const Q_DECL_OVERRIDE{return true;}
    QString computeStats(int type)Q_DECL_OVERRIDE;
    CGAL::Three::Scene_item::Header_data header() const Q_DECL_OVERRIDE;
    TextListItem* textItems;
    Scene_polyhedron_item();
    //   Scene_polyhedron_item(const Scene_polyhedron_item&);
    Scene_polyhedron_item(const Polyhedron& p);
    Scene_polyhedron_item(Polyhedron* const p);
    ~Scene_polyhedron_item();

    Scene_polyhedron_item* clone() const Q_DECL_OVERRIDE;

    // IO
    bool load(std::istream& in);
    bool load_obj(std::istream& in);
    bool save(std::ostream& out) const;
    bool save_obj(std::ostream& out) const;

    // Function for displaying meta-data of the item
    virtual QString toolTip() const Q_DECL_OVERRIDE;

    // Function to override the context menu
    QMenu* contextMenu() Q_DECL_OVERRIDE;

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE;
    // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
    void draw() const Q_DECL_OVERRIDE{}
    virtual void draw(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    virtual void drawEdges() const Q_DECL_OVERRIDE{}
    virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
    virtual void drawPoints(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;

    // Get wrapped polyhedron
    Polyhedron*       polyhedron();
    const Polyhedron* polyhedron() const;

    Face_graph*       face_graph() { return polyhedron(); }
    const Face_graph* face_graph() const { return polyhedron(); }

    // Get dimensions
    bool isFinite() const Q_DECL_OVERRIDE { return true; }
    bool isEmpty() const Q_DECL_OVERRIDE;
    void compute_bbox() const Q_DECL_OVERRIDE;
    
    Vertex_selection_map vertex_selection_map();
    Face_selection_map face_selection_map();
 
    std::vector<QColor>& color_vector();
    void set_color_vector_read_only(bool on_off);
    bool is_color_vector_read_only();

    int getNumberOfNullLengthEdges();
    int getNumberOfDegeneratedFaces();
    bool triangulated();
    bool self_intersected();
    //! If b is true, the item will use buffers to render the color.
    //! If b is false, it will use a uniform value. For example, when
    //! using the mesh segmentation plugin, the item must be multicolor.
    void setItemIsMulticolor(bool b);
    //! @returns `true` if the item has multiple colors at the same time.
    bool isItemMulticolor();

    void printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface*viewer)Q_DECL_OVERRIDE;
    void printPrimitiveIds(CGAL::Three::Viewer_interface*viewer) const Q_DECL_OVERRIDE;
    bool testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface*)const Q_DECL_OVERRIDE;


    //! @returns `true` if `f` is the first facet intersected by a raytracing
    bool intersect_face(double orig_x,
                        double orig_y,
                        double orig_z,
                        double dir_x,
                        double dir_y,
                        double dir_z,
                        Polyhedron::Facet_handle f);

public Q_SLOTS:
    virtual void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;
    virtual void selection_changed(bool) Q_DECL_OVERRIDE;
    virtual void setColor(QColor c) Q_DECL_OVERRIDE;
    virtual void show_feature_edges(bool);
    void show_only_feature_edges(bool);
    void enable_facets_picking(bool);
    void set_erase_next_picked_facet(bool);
    void set_flat_disabled(bool b);

    void select(double orig_x,
                double orig_y,
                double orig_z,
                double dir_x,
                double dir_y,
                double dir_z) Q_DECL_OVERRIDE;
    void update_vertex_indices();
    void update_facet_indices();
    void update_halfedge_indices();
    void invalidate_aabb_tree();
    void itemAboutToBeDestroyed(Scene_item *) Q_DECL_OVERRIDE;
    void resetColors();

Q_SIGNALS:
    void selection_done();
    void selected_vertex(void*);
    void selected_facet(void*);
    void selected_edge(void*);
    void selected_halfedge(void*);
    void item_is_about_to_be_changed(); // emitted in invalidateOpenGLBuffers()
public:
    typedef Scene_item Base;
    typedef Polyhedron::Facet_iterator Facet_iterator;
protected:
    friend struct Scene_polyhedron_item_priv;
    Scene_polyhedron_item_priv* d;

   public:
    void zoomToPosition(const QPoint &point, CGAL::Three::Viewer_interface *)const Q_DECL_OVERRIDE;

}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
