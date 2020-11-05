#ifndef SCENE_POLYGON_SOUP_ITEM_H
#define SCENE_POLYGON_SOUP_ITEM_H
#include "Scene_polygon_soup_item_config.h"
#include  <CGAL/Three/Scene_item_rendering_helper.h>
#include "SMesh_type.h"

#include <boost/array.hpp>

#include <iostream>

struct Scene_polygon_soup_item_priv;
struct Polygon_soup
{
    typedef EPICK::Point_3 Point_3;
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


class Scene_surface_mesh_item;

class SCENE_POLYGON_SOUP_ITEM_EXPORT Scene_polygon_soup_item
        : public CGAL::Three::Scene_item_rendering_helper
{
    Q_OBJECT
public:
    typedef EPICK::Point_3 Point_3;
    typedef Polygon_soup::Points Points;
    typedef Polygon_soup::Polygons Polygons;
    typedef Polygon_soup::Edges Edges;
    typedef Polygon_soup::Edge Edge;

    Scene_polygon_soup_item();
    ~Scene_polygon_soup_item();

    Scene_polygon_soup_item* clone() const Q_DECL_OVERRIDE;

    template <class Point, class Polygon>
    void load(const std::vector<Point>& points, const std::vector<Polygon>& polygons);

    template <class Point, class Polygon>
    void load(const std::vector<Point>& points, const std::vector<Polygon>& polygons,
              const std::vector<CGAL::Color>& fcolors,
              const std::vector<CGAL::Color>& vcolors);

    bool load(std::istream& in);
    void load(Scene_surface_mesh_item*);
    bool isDataColored();

    bool save(std::ostream& out) const;
    std::vector<CGAL::Color> getVColors() const;
    std::vector<CGAL::Color> getFColors() const;
    QString toolTip() const Q_DECL_OVERRIDE;

    // Indicate if rendering mode is supported
    virtual bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE{ return ( m!=PointsPlusNormals && m!=ShadedPoints); }
    // OpenGL drawing in a display list
    virtual void draw() const Q_DECL_OVERRIDE{}
    virtual void draw(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    virtual void drawPoints(CGAL::Three::Viewer_interface*) const Q_DECL_OVERRIDE;
    virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE;
    void invalidateOpenGLBuffers() Q_DECL_OVERRIDE;
    bool isFinite() const Q_DECL_OVERRIDE{ return true; }
    bool isEmpty() const Q_DECL_OVERRIDE;
    void compute_bbox() const Q_DECL_OVERRIDE;

    void new_vertex(const double&, const double&, const double&);
    void new_triangle(const std::size_t, const std::size_t, const std::size_t);

    void init_polygon_soup(std::size_t nb_pts, std::size_t nb_polygons);

    const Points& points() const;
    const Polygons& polygons() const;
    const Edges& non_manifold_edges() const;
    void initializeBuffers(CGAL::Three::Viewer_interface *) const Q_DECL_OVERRIDE;
    void computeElements() const Q_DECL_OVERRIDE;
    //statistics
    enum STATS {
      NB_VERTICES = 0,
      NB_FACETS,
      IS_PURE_TRIANGLE,
      IS_PURE_QUAD,
      NB_DEGENERATED_FACES,
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
public Q_SLOTS:
    void shuffle_orientations();
    bool orient();
    bool exportAsSurfaceMesh(SMesh*);
    void inside_out();
    void repair(bool erase_dup, bool req_same_orientation);

    void setDisplayNonManifoldEdges(const bool);
    bool displayNonManifoldEdges() const;
    void itemAboutToBeDestroyed(Scene_item *item) Q_DECL_OVERRIDE;

protected:
    friend struct Scene_polygon_soup_item_priv;
    Scene_polygon_soup_item_priv* d;

}; // end class Scene_polygon_soup_item

#endif // SCENE_POLYGON_SOUP_ITEM_H
