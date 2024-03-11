/**
 * @file   data/3d/Edge.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_EDGE_H
#define DATA_3D_EDGE_H

#include "data/3d/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class Edge : public std::enable_shared_from_this<Edge>  {
public:
    virtual ~Edge();

    static EdgeSPtr create(VertexSPtr src, VertexSPtr dst);

    EdgeSPtr clone() const;

    VertexSPtr getVertexSrc() const;
    void setVertexSrc(VertexSPtr src);
    std::list<EdgeWPtr>::iterator getVertexSrcListIt() const;
    void setVertexSrcListIt(std::list<EdgeWPtr>::iterator list_it);

    VertexSPtr getVertexDst() const;
    void setVertexDst(VertexSPtr dst);
    std::list<EdgeWPtr>::iterator getVertexDstListIt() const;
    void setVertexDstListIt(std::list<EdgeWPtr>::iterator list_it);

    FacetSPtr getFacetL() const;
    void setFacetL(FacetSPtr facet);
    std::list<EdgeSPtr>::iterator getFacetLListIt() const;
    void setFacetLListIt(std::list<EdgeSPtr>::iterator list_it);

    FacetSPtr getFacetR() const;
    void setFacetR(FacetSPtr facet);
    std::list<EdgeSPtr>::iterator getFacetRListIt() const;
    void setFacetRListIt(std::list<EdgeSPtr>::iterator list_it);

    PolyhedronSPtr getPolyhedron() const;
    void setPolyhedron(PolyhedronSPtr polyhedron);
    std::list<EdgeSPtr>::iterator getPolyhedronListIt() const;
    void setPolyhedronListIt(std::list<EdgeSPtr>::iterator list_it);

    EdgeDataSPtr getData() const;
    void setData(EdgeDataSPtr data);
    bool hasData() const;

    Segment3SPtr segment() const;
    Line3SPtr line() const;

    FacetSPtr other(FacetSPtr facet) const;
    VertexSPtr src(FacetSPtr facet_l) const;
    VertexSPtr dst(FacetSPtr facet_l) const;
    FacetSPtr left(VertexSPtr vertex_src) const;
    FacetSPtr right(VertexSPtr vertex_src) const;
    EdgeSPtr next(FacetSPtr facet) const;
    EdgeSPtr prev(FacetSPtr facet) const;

    /**
     * counter clockwise from outside
     */
    EdgeSPtr next(VertexSPtr vertex) const;
    EdgeSPtr prev(VertexSPtr vertex) const;

    /**
     * Swaps src and dst vertex and left and right facet
     */
    void invert();

    /**
     * Splits this edge into 2 edges.
     * The destination vertex of this edge is set to the given middle vertex.
     * It returns newly created edge.
     */
    EdgeSPtr split(VertexSPtr middle);

    /**
     * More than just a simple set method.
     */
    void replaceVertexSrc(VertexSPtr vertex_src);
    void replaceVertexDst(VertexSPtr vertex_dst);
    void replaceFacetL(FacetSPtr facet_l);
    void replaceFacetR(FacetSPtr facet_r);

    bool hasSameFacets(EdgeSPtr edge);

    int getID() const;
    void setID(int id);

    double angle() const;
    bool isReflex() const;

    double angleTo(EdgeSPtr edge) const;

    std::string toString() const;

protected:
    Edge(VertexSPtr src, VertexSPtr dst);
    Edge(const Edge& edge);
    VertexSPtr vertex_src_;
    std::list<EdgeWPtr>::iterator vertex_src_list_it_;
    VertexSPtr vertex_dst_;
    std::list<EdgeWPtr>::iterator vertex_dst_list_it_;
    FacetWPtr facet_l_;
    std::list<EdgeSPtr>::iterator facet_l_list_it_;
    FacetWPtr facet_r_;
    std::list<EdgeSPtr>::iterator facet_r_list_it_;
    PolyhedronWPtr polyhedron_;
    std::list<EdgeSPtr>::iterator polyhedron_list_it_;
    EdgeDataSPtr data_;

    int id_;
};

} }

#endif /* DATA_3D_EDGE_H */
