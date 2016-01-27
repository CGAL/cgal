#ifndef SCENE_POLYHEDRON_ITEM_H
#define SCENE_POLYHEDRON_ITEM_H

#include "Scene_polyhedron_item_config.h"
#include  <CGAL/Three/Scene_item.h> //<- modif ?
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
        : public CGAL::Three::Scene_item{
    Q_OBJECT
public:
    enum STATS {
      NB_VERTICES = 0,
      NB_FACETS,
      NB_CONNECTED_COMPOS,
      NB_BORDER_EDGES,
      NB_DEGENERATED_FACES,
      HOLES,
      AREA,
      VOLUME,
      SELFINTER,
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
    QString compute_stats(int type);
    CGAL::Three::Scene_item::Header_data header() const;
    Scene_polyhedron_item();
    //   Scene_polyhedron_item(const Scene_polyhedron_item&);
    Scene_polyhedron_item(const Polyhedron& p);
    Scene_polyhedron_item(Polyhedron* const p);
    ~Scene_polyhedron_item();

    Scene_polyhedron_item* clone() const;

    // IO
    bool load(std::istream& in);
    bool load_obj(std::istream& in);
    bool save(std::ostream& out) const;

    // Function for displaying meta-data of the item
    virtual QString toolTip() const;

    // Function to override the context menu
    QMenu* contextMenu();

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const { return (m!=PointsPlusNormals && m!=Splatting); }
    // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
    void draw() const {}
    virtual void draw(CGAL::Three::Viewer_interface*) const;
    virtual void draw_edges() const {}
    virtual void draw_edges(CGAL::Three::Viewer_interface* viewer) const;
    virtual void draw_points(CGAL::Three::Viewer_interface*) const;

    // Get wrapped polyhedron
    Polyhedron*       polyhedron();
    const Polyhedron* polyhedron() const;

    // Get dimensions
    bool isFinite() const { return true; }
    bool isEmpty() const;
    void compute_bbox() const;
    std::vector<QColor>& color_vector() {return colors_;}
    void set_color_vector_read_only(bool on_off) {plugin_has_set_color_vector_m=on_off;}
    int getNumberOfNullLengthEdges(){return number_of_null_length_edges;}
    int getNumberOfDegeneratedFaces(){return number_of_degenerated_faces;}
    bool triangulated(){return poly->is_pure_triangle();}
    bool self_intersected(){return !self_intersect;}

public Q_SLOTS:
    virtual void invalidateOpenGLBuffers();
    virtual void selection_changed(bool);
    virtual void setColor(QColor c);
	virtual void show_feature_edges(bool);
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
    void item_is_about_to_be_changed(); // emitted in invalidateOpenGLBuffers()

private:
    // Initialization
    void init();
    void invalidate_stats();

private:
    Polyhedron* poly;

private:
    typedef Scene_item Base;
    typedef std::vector<QColor> Color_vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;


    Color_vector colors_;
    bool show_only_feature_edges_m;
    bool show_feature_edges_m;
    bool facet_picking_m;
    bool erase_next_picked_facet_m;
    //the following variable is used to indicate if the color vector must not be automatically updated.
    bool plugin_has_set_color_vector_m;

    enum VAOs {
        Facets=0,
        Edges,
        Feature_edges,
        Gouraud_Facets,
        NbOfVaos = Gouraud_Facets+1
    };
    enum VBOs {
        Facets_vertices = 0,
        Facets_normals_flat,
        Facets_color,
        Edges_vertices,
        Feature_edges_vertices,
        Edges_color,
        Facets_normals_gouraud,
        NbOfVbos = Facets_normals_gouraud+1
    };

    mutable std::vector<float> positions_lines;
    mutable std::vector<float> positions_feature_lines;
    mutable std::vector<float> positions_facets;
    mutable std::vector<float> normals_flat;
    mutable std::vector<float> normals_gouraud;
    mutable std::vector<float> color_lines;
    mutable std::vector<float> color_facets;
    mutable std::size_t nb_facets;
    mutable std::size_t nb_lines;
    mutable std::size_t nb_f_lines;
    mutable QOpenGLShaderProgram *program;
    mutable unsigned int number_of_null_length_edges;
    mutable unsigned int number_of_degenerated_faces;
    mutable bool self_intersect;

    using CGAL::Three::Scene_item::initialize_buffers;
    void initialize_buffers(CGAL::Three::Viewer_interface *viewer = 0) const;
    void compute_normals_and_vertices(void) const;
    void compute_colors() const;
    void triangulate_facet(Facet_iterator ) const;
    void triangulate_facet_color(Facet_iterator ) const;
    double volume, area;

}; // end class Scene_polyhedron_item

#endif // SCENE_POLYHEDRON_ITEM_H
