#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include "Scene_item.h" //<- modif ?
#include "Polyhedron_type_fwd.h"
#include "Polyhedron_type.h"
#include "Viewer.h"
#include <iostream>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>

#include <set>
#include <vector>

#include <QColor>

class QMenu;

// This class represents a polyhedron in the OpenGL scene
class SCENE_POLYHEDRON_ITEM_EXPORT Scene_polyhedron_item 
        : public Scene_item{
    Q_OBJECT
public:  
    Scene_polyhedron_item();
    //   Scene_polyhedron_item(const Scene_polyhedron_item&);
    Scene_polyhedron_item(const Polyhedron& p);
    Scene_polyhedron_item(Polyhedron* const p);
    ~Scene_polyhedron_item();

    Scene_polyhedron_item* clone() const;

    // IO
    bool load(std::istream& in);
    bool save(std::ostream& out) const;
    mutable bool is_Triangle;

    // Function for displaying meta-data of the item
    virtual QString toolTip() const;

    // Function to override the context menu
    QMenu* contextMenu();

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals && m!=Splatting); }
    // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
    void draw() const {}
    virtual void draw(Viewer_interface*) const;
    virtual void draw_edges() const {}
    virtual void draw_edges(Viewer_interface* viewer) const;
    virtual void draw_points(Viewer_interface*) const;

    // Get wrapped polyhedron
    Polyhedron*       polyhedron();
    const Polyhedron* polyhedron() const;

    // Get dimensions
    bool isFinite() const { return true; }
    bool isEmpty() const;
    Bbox bbox() const;
    std::vector<QColor>& color_vector() {return colors_;}
    void set_color_vector_read_only(bool on_off) {plugin_has_set_color_vector_m=on_off;}

public Q_SLOTS:
    virtual void invalidate_buffers();
    virtual void contextual_changed();
    virtual void selection_changed(bool);
    virtual void setColor(QColor c);
    void show_only_feature_edges(bool);
    void enable_facets_picking(bool);
    void set_erase_next_picked_facet(bool);

    void select(double orig_x,
                double orig_y,
                double orig_z,
                double dir_x,
                double dir_y,
                double dir_z);

    void update_vertex_indices();
    void update_facet_indices();
    void update_halfedge_indices();
    void invalidate_aabb_tree();

Q_SIGNALS:
    void selected_vertex(void*);
    void selected_facet(void*);
    void selected_edge(void*);
    void selected_halfedge(void*);
    void item_is_about_to_be_changed(); // emitted in invalidate_buffers()

private:
    // Initialization
    void init();

private:
    Polyhedron* poly;

private:
    typedef Scene_item Base;
    typedef std::vector<QColor> Color_vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;

    Color_vector colors_;

    bool show_only_feature_edges_m;
    bool facet_picking_m;
    bool erase_next_picked_facet_m;
    //the following variable is used to indicate if the color vector must not be automatically updated.
    bool plugin_has_set_color_vector_m;


    mutable std::vector<float> positions_lines;
    mutable std::vector<float> positions_facets;
    mutable std::vector<float> normals_flat;
    mutable std::vector<float> normals_gouraud;
    mutable std::vector<float> color_lines;
    mutable std::vector<float> color_facets;
    mutable std::vector<float> color_lines_selected;
    mutable std::vector<float> color_facets_selected;
    mutable std::size_t nb_facets;
    mutable std::size_t nb_lines;
    mutable QOpenGLShaderProgram *program;

    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer = 0) const;
    void compute_normals_and_vertices(void) const;
    void compute_colors() const;
    void triangulate_facet(Facet_iterator ) const;
    void triangulate_facet_color(Facet_iterator ) const;
    void is_Triangulated() const;
    double volume, area;

}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
