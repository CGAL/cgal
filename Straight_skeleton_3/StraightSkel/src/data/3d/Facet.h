/**
 * @file   data/3d/Facet.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_FACET_H
#define DATA_3D_FACET_H

#include "data/2d/ptrs.h"
#include "data/3d/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class Facet : public std::enable_shared_from_this<Facet> {
public:
    virtual ~Facet();

    static FacetSPtr create();

    static FacetSPtr create(unsigned int num_vertices, VertexSPtr vertices[]);

    /**
     * Each edge may have a facet to the left and a facet to the right.
     * Left and right is defined by a view outside the polyhedron.
     * In case an edge is used the first time, it is assumed that the left
     * side of the edge points to the interior of the facet.
     * In case the edge is already part of a facet,
     * it is assumed that the right side of the edge points to
     * the interior of the facet that is created.
     */
    static FacetSPtr create(unsigned int num_edges, EdgeSPtr edges[]);

    FacetSPtr clone() const;

    void addVertex(VertexSPtr vertex);
    bool removeVertex(VertexSPtr vertex);

    void addEdge(EdgeSPtr edge);
    bool removeEdge(EdgeSPtr edge);

    /**
     * Searches for an edge to the given facet.
     * Does not work when there is more than one edge adjacent to both facets.
     * @deprecated
     */
    EdgeSPtr findEdge(FacetSPtr facet) const;

    std::list<EdgeSPtr> findEdges(FacetSPtr facet) const;

    bool containsVertex(VertexSPtr vertex) const;
    bool containsEdge(EdgeSPtr edge) const;

    void sortVertices();
    void sortEdges();

    void addTriangle(TriangleSPtr triangle);
    bool removeTriangle(TriangleSPtr triangle);

    PolyhedronSPtr getPolyhedron() const;
    void setPolyhedron(PolyhedronSPtr polyhedron);
    std::list<FacetSPtr>::iterator getPolyhedronListIt() const;
    void setPolyhedronListIt(std::list<FacetSPtr>::iterator list_it);

    FacetDataSPtr getData() const;
    void setData(FacetDataSPtr data);
    bool hasData() const;

    std::list<VertexSPtr>& vertices();
    std::list<EdgeSPtr>& edges();
    std::list<TriangleSPtr>& triangles();

    /**
     * counter clockwise from outside
     */
    FacetSPtr next(VertexSPtr vertex) const;
    FacetSPtr prev(VertexSPtr vertex) const;

    void merge(FacetSPtr facet);

    int getID() const;
    void setID(int id);

    /**
     * The direction of the normal points to the outside.
     */
    Plane3SPtr getPlane() const;
    void setPlane(Plane3SPtr plane);

    /**
     * First vertices have to form a triangle that is inside.
     */
    bool initPlane();
    Plane3SPtr plane();

    bool makeFirstConvex();

    data::_2d::PolygonSPtr toPolygon();

    std::string toString() const;

protected:
    Facet();
    std::list<VertexSPtr> vertices_;
    std::list<EdgeSPtr> edges_;
    std::list<TriangleSPtr> triangles_;
    PolyhedronWPtr polyhedron_;
    std::list<FacetSPtr>::iterator polyhedron_list_it_;
    FacetDataSPtr data_;
    int id_;
    Plane3SPtr plane_;
};

} }

#endif /* DATA_3D_FACET_H */
