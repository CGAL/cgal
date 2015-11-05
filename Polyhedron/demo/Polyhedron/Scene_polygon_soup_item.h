#ifndef SCENE_POLYGON_SOUP_ITEM_H
#define SCENE_POLYGON_SOUP_ITEM_H
#include "Scene_polygon_soup_item_config.h"
#include "Scene_item.h"
#include "Viewer.h"
#include "Polyhedron_type.h"

#include <boost/foreach.hpp>
#include <boost/array.hpp>

#include <iostream>


struct Polygon_soup
{
    typedef Kernel::Point_3 Point_3;
    typedef std::vector<Point_3> Points;
    //vector containing 3 indices of points in Points
    typedef std::vector<std::size_t> Polygon_3;
    //vector containing a pair of indices of points in Points and a set of indices of Polygons
    //containing the edge.
    typedef std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t> > Edges_map;
    typedef boost::array<std::size_t, 2> Edge;
    typedef std::vector<Polygon_3> Polygons;
    typedef std::set<Edge> Edges;
    typedef Polygons::size_type size_type;
    Points points;
    Polygons polygons;
    Edges_map edges;
    Edges non_manifold_edges;
    bool display_non_manifold_edges;

    Polygon_soup():
        display_non_manifold_edges(false){}

    Polygon_soup* clone() const {
        Polygon_soup* result = new Polygon_soup();
        result->points = points;
        result->polygons = polygons;
        result->edges = edges;
        result->non_manifold_edges = non_manifold_edges;
        result->display_non_manifold_edges = display_non_manifold_edges;
        return result;
    }

    void clear() {
        points.clear();
        polygons.clear();
        edges.clear();
        non_manifold_edges.clear();
    }

    void fill_edges() {
        // Fill edges
        edges.clear();
        for(size_type i = 0; i < polygons.size(); ++i)
        {
            const size_type size = polygons[i].size();
            for(size_type j = 0; j < size; ++j) {
                const std::size_t& i0 = polygons[i][j];
                const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
                edges[std::make_pair(i0, i1)].insert(i);
                //         qDebug() << tr("edges[std::make_pair(%1, %2)].insert(%3). Size=%4")
                //           .arg(i0).arg(i1).arg(i).arg(edges[std::make_pair(i0, i1)].size());
            }
        }

        // Fill non-manifold edges
        non_manifold_edges.clear();
        for(size_type i = 0; i < polygons.size(); ++i)
        {
            const size_type size = polygons[i].size();
            for(size_type j = 0; j < size; ++j) {
                const std::size_t& i0 = polygons[i][j];
                const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];

                if(edges[std::make_pair(i0, i1)].size() +
                         edges[std::make_pair(i1, i0)].size() > 2)
                {
                    Edge edge;
                    edge[0] = i0;
                    edge[1] = i1;
                    if(i0 > i1) std::swap(edge[0], edge[1]);
                    non_manifold_edges.insert(edge);

                }
            }
        }
    }

    void inverse_orientation(const std::size_t index) {
        std::reverse(polygons[index].begin(), polygons[index].end());
    }
};


class Scene_polyhedron_item;

class SCENE_POLYGON_SOUP_ITEM_EXPORT Scene_polygon_soup_item 
        : public Scene_item
{
    typedef Kernel::Point_3 Point_3;

    Q_OBJECT
public:  
    Scene_polygon_soup_item();
    ~Scene_polygon_soup_item();

    Scene_polygon_soup_item* clone() const;
    bool load(std::istream& in);
    void load(Scene_polyhedron_item*);

    template <class Point, class Polygon>
    inline void load(const std::vector<Point>& points, const std::vector<Polygon>& polygons)
    {
        if(!soup)
            soup = new Polygon_soup;
        soup->clear();

        /// add points
        soup->points.reserve(points.size());
        BOOST_FOREACH(const Point& p, points)
                soup->points.push_back( Point_3(p[0], p[1], p[2]) );

        /// add polygons
        std::size_t nb_polygons=polygons.size();
        soup->polygons.resize(nb_polygons);
        for(std::size_t i=0; i<nb_polygons; ++i)
            soup->polygons[i].assign(polygons[i].begin(), polygons[i].end());

        /// fill non-manifold edges container
        //soup->fill_edges();
        oriented = false;

        Q_EMIT invalidate_buffers();
    }

    bool save(std::ostream& out) const;

    QString toolTip() const;

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const { return (m!=Gouraud && m!=PointsPlusNormals && m!=Splatting); } // CHECK THIS!
    // OpenGL drawing in a display list
    virtual void draw() const {}
    virtual void draw(Viewer_interface*) const;
    virtual void draw_points(Viewer_interface*) const;
    virtual void draw_edges(Viewer_interface* viewer) const;
    void invalidate_buffers();
    bool isFinite() const { return true; }
    bool isEmpty() const;
    Bbox bbox() const;

    void new_vertex(const double&, const double&, const double&);
    void new_triangle(const std::size_t, const std::size_t, const std::size_t);

    void init_polygon_soup(std::size_t nb_pts, std::size_t nb_polygons);

public Q_SLOTS:
    void shuffle_orientations();
    bool orient();
    bool exportAsPolyhedron(Polyhedron*);
    void inside_out();

    void setDisplayNonManifoldEdges(const bool);
    bool displayNonManifoldEdges() const;

private:
    typedef Polygon_soup::Polygons::const_iterator Polygons_iterator;
    Polygon_soup* soup;
    bool oriented;
    mutable std::vector<float> positions_poly;
    mutable std::vector<float> positions_lines;
    mutable std::vector<float> normals;
    mutable std::vector<float> positions_nm_lines;
    mutable std::size_t nb_nm_edges;
    mutable std::size_t nb_polys;
    mutable std::size_t nb_lines;
    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer) const;
    void compute_normals_and_vertices(void) const;
    void triangulate_polygon(Polygons_iterator ) const;
    mutable QOpenGLShaderProgram *program;

}; // end class Scene_polygon_soup_item

#endif // SCENE_POLYGON_SOUP_ITEM_H
