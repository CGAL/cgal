#ifndef SCENE_POLYGON_SOUP_ITEM_H
#define SCENE_POLYGON_SOUP_ITEM_H
#include "Scene_polygon_soup_item_config.h"
#include  <CGAL/Three/Scene_item.h>
#include "Polyhedron_type.h"

//#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include "CGAL/Surface_mesh/Surface_mesh.h"

#include <boost/foreach.hpp>
#include <boost/array.hpp>

#include <iostream>

struct Scene_polygon_soup_item_priv;
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
    typedef std::vector<CGAL::Color> Colors;
    typedef std::set<Edge> Edges;
    typedef Polygons::size_type size_type;
    Points points;
    Polygons polygons;
    Edges_map edges;
    Colors fcolors;
    Colors vcolors;
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
        : public CGAL::Three::Scene_item
{
    typedef Kernel::Point_3 Point_3;
    typedef Polygon_soup::Points Points;

    Q_OBJECT
public:  
    Scene_polygon_soup_item();
    ~Scene_polygon_soup_item();

    Scene_polygon_soup_item* clone() const;

    template <class Point, class Polygon>
    void load(const std::vector<Point>& points, const std::vector<Polygon>& polygons);

    bool load(std::istream& in);
    void load(Scene_polyhedron_item*);
    bool isDataColored();

    bool save(std::ostream& out) const;
    std::vector<CGAL::Color> getVColors() const;
    std::vector<CGAL::Color> getFColors() const;
    QString toolTip() const;

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const { return ( m!=PointsPlusNormals && m!=Splatting && m!=ShadedPoints); }
    // OpenGL drawing in a display list
    virtual void draw() const {}
    virtual void draw(CGAL::Three::Viewer_interface*) const;
    virtual void drawPoints(CGAL::Three::Viewer_interface*) const;
    virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
    void invalidateOpenGLBuffers();
    bool isFinite() const { return true; }
    bool isEmpty() const;
    void compute_bbox() const;

    void new_vertex(const double&, const double&, const double&);
    void new_triangle(const std::size_t, const std::size_t, const std::size_t);

    void init_polygon_soup(std::size_t nb_pts, std::size_t nb_polygons);

    const Points& points() const;
public Q_SLOTS:
    void shuffle_orientations();
    bool orient();
    bool exportAsPolyhedron(Polyhedron*);
    bool exportAsSurfaceMesh(CGAL::Surface_mesh<Point_3>*);
    void inside_out();

    void setDisplayNonManifoldEdges(const bool);
    bool displayNonManifoldEdges() const;

protected:
    friend struct Scene_polygon_soup_item_priv;
    Scene_polygon_soup_item_priv* d;

}; // end class Scene_polygon_soup_item

#endif // SCENE_POLYGON_SOUP_ITEM_H
