#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include <CGAL/AABB_intersections.h> 
#include "types.h"
#include "Color_ramp.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <QtCore/qglobal.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLFunctions_2_1>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>

class Texture{
private:
     int Width;
     int Height;
     int size;
    GLubyte *data;
public:
    Texture(int w, int h)
    {
        Width = w;
        Height = h;
        size = 3*Height*Width;
        data = new GLubyte[Height*Width*3];
    }
    int getWidth() const {return Width;}
    int getHeight() const {return Height;}
    int getSize() const {return size;}
    void setData(int i, int j, int r, int g, int b){
        data[3*(Width*j+i) + 0] = r;
        data[3*(Width*j+i) + 1] = g;
        data[3*(Width*j+i) + 2] = b;
    }

    GLubyte* getData(){return data; }

};
class Viewer;
class Scene : public QObject
{
    Q_OBJECT
public:
    Scene();
    virtual ~Scene();
public:
    // types
    typedef CGAL::Bbox_3 Bbox;
    
private:
    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>         Facet_Primitive;
    typedef CGAL::AABB_traits<Kernel, Facet_Primitive>                  Facet_Traits;
    typedef CGAL::AABB_tree<Facet_Traits>                               Facet_tree;
  
    typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron>     Edge_Primitive;
    typedef CGAL::AABB_traits<Kernel, Edge_Primitive>                   Edge_Traits;
    typedef CGAL::AABB_tree<Edge_Traits>                                Edge_tree;
    
    typedef qglviewer::ManipulatedFrame ManipulatedFrame;
    
    enum Cut_planes_types {
        NONE, UNSIGNED_FACETS, SIGNED_FACETS, UNSIGNED_EDGES, CUT_SEGMENTS
    };
  
public:
    QGLContext* context;
    void draw(QGLViewer*);
    void update_bbox();
    Bbox bbox() { return m_bbox; }
    ManipulatedFrame* manipulatedFrame() const { return m_frame; }
    void initGL(Viewer *viewer);

private:
    // member data
    QOpenGLFunctions_2_1 *gl;
    Bbox m_bbox;
    Polyhedron *m_pPolyhedron;
    std::list<Point> m_points;
    std::list<Segment> m_segments;
    std::vector<Segment> m_cut_segments;
    bool ready_to_cut;

    // distance functions (simple 2D arrays)
    Color_ramp m_red_ramp;
    Color_ramp m_blue_ramp;
    Color_ramp m_thermal_ramp;
    FT m_max_distance_function;
    typedef std::pair<Point,FT> Point_distance;
    Point_distance m_distance_function[100][100];
  
    // frame
    ManipulatedFrame* m_frame;
    bool m_view_plane;
    int m_grid_size;
    bool m_fast_distance;

    // An aabb_tree indexing polyhedron facets/segments
    Facet_tree m_facet_tree;
    Edge_tree m_edge_tree;
    
    Cut_planes_types m_cut_plane;
    bool are_buffers_initialized;
  
private:
    // utility functions
    Vector random_vector();
    Ray random_ray(const Bbox& bbox);
    Line random_line(const Bbox& bbox);
    Point random_point(const Bbox& bbox);
    Plane random_plane(const Bbox& bbox);
    Segment random_segment(const Bbox& bbox);
    FT random_in(const double a,const double b);
    Plane frame_plane() const;
    Aff_transformation frame_transformation() const;
    FT bbox_diag() const;
    void build_facet_tree();
    void build_edge_tree();
    void clear_internal_data();
    void update_grid_size();

    template <typename Tree>
    void compute_distance_function(const Tree& tree);
    
    template <typename Tree>
    void sign_distance_function(const Tree& tree);

    //Shaders elements

    int poly_vertexLocation;
    int tex_Location;
    int points_vertexLocation;
    int lines_vertexLocation;
    int mvpLocation;
    int tex_mvpLocation;
    int fLocation;
    int tex_fLocation;
    int colorLocation;



    std::vector<float> pos_points;
    std::vector<float> pos_grid;
    std::vector<float> pos_lines;
    std::vector<float> pos_poly;
    std::vector<float> pos_plane;
    std::vector<float> pos_cut_segments;
    std::vector<float> tex_map;
    GLuint textureId;

    Texture *texture;
    GLint sampler_location;
    QOpenGLBuffer buffers[10];
    QOpenGLVertexArrayObject vao[10];
    QOpenGLShaderProgram tex_rendering_program;
    QOpenGLShaderProgram rendering_program;
    void initialize_buffers();
    void compute_elements(int mode);
    void attrib_buffers(QGLViewer*);
    void compile_shaders();
    void compute_texture(int, int, Color_ramp, Color_ramp);

public:
    // file menu
    int open(QString filename);

    // edit menu
    void clear_points() { m_points.clear(); changed(); }
    void clear_segments() { m_segments.clear(); changed(); }
    void clear_cutting_plane();
    
    // fast distance setter
    void set_fast_distance(bool b) { m_fast_distance = b; update_grid_size(); }

    // algorithms
    void generate_edge_points(const unsigned int nb_points);
    void generate_inside_points(const unsigned int nb_points);
    void generate_boundary_points(const unsigned int nb_points);
    void generate_boundary_segments(const unsigned int nb_slices);
    void generate_points_in(const unsigned int nb_points,
        const double min, const double max);

    // algorithms/refine
    void refine_loop();
    void refine_bisection(const FT max_sqlen);

    // distance functions 
    void signed_distance_function();
    void unsigned_distance_function();
    void unsigned_distance_function_to_edges();
    void cut_segment_plane();

    // toggle view options
    void toggle_view_points();
    void toggle_view_segments();
    void toggle_view_poyhedron();
    void toggle_view_plane();

    // view options
    bool m_view_points;
    bool m_view_segments;
    bool m_view_polyhedron;

    // benchmarks 
    enum {DO_INTERSECT,
        ANY_INTERSECTION,
        NB_INTERSECTIONS,
        ALL_INTERSECTIONS,
        ANY_INTERSECTED_PRIMITIVE,
        ALL_INTERSECTED_PRIMITIVES};
    void bench_memory();
    void bench_construction();
    void bench_distances_vs_nbt();
    void bench_intersections_vs_nbt();
    void benchmark_intersections(const double duration);
    std::size_t nb_digits(const std::size_t value);

    template <class Query>
    void bench_intersection(Facet_tree& tree,const int function,const double duration,
        const char *query_name, const std::vector<Query>& queries, const int nb_queries);
    void bench_intersections(Facet_tree& tree, const double duration, const int function,
        const char *function_name, const std::vector<Ray>& rays,
        const std::vector<Line>& lines, const std::vector<Plane>& planes,
        const std::vector<Segment>& segments, const int nb_queries);

    // distance benchmarks
    enum {SQ_DISTANCE,
        CLOSEST_POINT,
        CLOSEST_POINT_AND_PRIMITIVE_ID};
    void benchmark_distances(const double duration);
    void bench_closest_point(Facet_tree& tree,const double duration);
    void bench_squared_distance(Facet_tree& tree,const double duration);
    void bench_closest_point_and_primitive(Facet_tree& tree,const double duration);
    void bench_distance(Facet_tree& tree,const int function,const double duration);

    // cutting plane activation/deactivation
    void activate_cutting_plane();
    void deactivate_cutting_plane();

    //timer sends a top when all the events are finished
    void timerEvent(QTimerEvent *);

  
public slots:
    // cutting plane
    void cutting_plane(bool override = false);
    void changed();
}; // end class Scene

#endif // SCENE_H
