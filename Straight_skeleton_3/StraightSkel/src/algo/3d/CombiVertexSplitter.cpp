/**
 * @file   algo/3d/CombiVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-07-26
 */

#include "algo/3d/CombiVertexSplitter.h"

#include "debug.h"
#include "typedefs_thread.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/SkelVertexData.h"
#include "data/3d/skel/SkelFacetData.h"
#include "util/Configuration.h"
#include <map>
#include <sstream>

namespace algo { namespace _3d {

CombiVertexSplitter::CombiVertexSplitter() {
    type_ = AbstractVertexSplitter::COMBI_VERTEX_SPLITTER;
    selected_combi_ = 0;
    util::ConfigurationSPtr config = util::Configuration::getInstance();
    if (config->isLoaded()) {
        selected_combi_ = config->getInt(
                "algo_3d_CombiVertexSplitter", "selected_combi");
    }
}

CombiVertexSplitter::~CombiVertexSplitter() {
    // intentionally does nothing
}

CombiVertexSplitterSPtr CombiVertexSplitter::create() {
    CombiVertexSplitterSPtr result =
            CombiVertexSplitterSPtr(new CombiVertexSplitter());
    return result;
}

vec2i CombiVertexSplitter::createSplit(int begin, int end) {
    vec2i result(new int[2]);
    result[0] = begin;
    result[1] = end;
    return result;
}

int CombiVertexSplitter::compareSplits(vec2i split1, vec2i split2) {
    int result = 0;
    if (split1[0] < split2[0] ||
            (split1[0] == split2[0] && split1[1] < split2[1])) {
        result = 1;
    } else if (split1[0] > split2[0] ||
            (split1[0] == split2[0] && split1[1] > split2[1])) {
        result = -1;
    }
    return result;
}

std::vector<int> CombiVertexSplitter::initLabels(unsigned int degree) {
    std::vector<int> result;
    for (unsigned int i = 0; i < degree; i++) {
        result.push_back(i);
    }
    return result;
}

std::vector<int> CombiVertexSplitter::splitLabels(std::vector<int>& labels, vec2i split) {
    std::vector<int> result;
    int begin = split[0];
    int end = split[1];
    bool inside = false;
    std::vector<int>::iterator it = labels.begin();
    while (it != labels.end()) {
        std::vector<int>::iterator it_current = it;
        int label = *it++;
        if (label == begin) {
            inside = true;
        }
        if (inside) {
            result.push_back(label);
            if (label != begin && label != end) {
                it = labels.erase(it_current);
            }
        }
        if (label == end) {
            inside = false;
        }
        if (inside) {
            if (it == labels.end()) {
                it = labels.begin();
            }
        }
    }
    return result;
}

std::list<vec2i> CombiVertexSplitter::createSingleSplitCombinations(std::vector<int> labels) {
    std::list<vec2i> result;
    unsigned int degree = labels.size();
    for (unsigned int i = 0; i < degree-1; i++) {
        for (unsigned int j = i+2; j < degree; j++) {
            if (i == 0 && j == degree-1) {
                continue;
            }
            vec2i split = createSplit(labels[i], labels[j]);
            result.push_back(split);
        }
    }
    return result;
}

std::list<combi> CombiVertexSplitter::appendSplitCombinations(
        combi history, std::list<vec2i> splits) {
    std::list<combi> result;
    if (history.size() == 0) {
        std::list<vec2i>::iterator it_splits = splits.begin();
        while (it_splits != splits.end()) {
            vec2i split = *it_splits++;
            combi combination;
            combination.push_back(split);
            result.push_back(combination);
        }
    } else {
        vec2i last_split = history.back();
        std::list<vec2i>::iterator it_splits = splits.begin();
        while (it_splits != splits.end()) {
            vec2i split = *it_splits++;
            if (compareSplits(last_split, split) > 0) {
                combi combination(history);
                combination.push_back(split);
                result.push_back(combination);
            }
        }
    }
    return result;
}

std::list<combi> CombiVertexSplitter::mergeCombinations(combi history,
        std::list<combi> combis1, std::list<combi> combis2) {
    std::list<combi> result;
    unsigned int history_size = history.size();
    std::list<combi>::iterator it_combis1 = combis1.begin();
    while (it_combis1 != combis1.end()) {
        combi combi1 = *it_combis1++;
        std::list<combi>::iterator it_combis2 = combis2.begin();
        while (it_combis2 != combis2.end()) {
            combi combi2 = *it_combis2++;
            combi combi_merged(history);
            std::vector<vec2i>::iterator it_combi1 = combi1.begin();
            std::vector<vec2i>::iterator it_combi2 = combi2.begin();
            for (unsigned int i = 0; i < history_size; i++) {
                it_combi1++;
                it_combi2++;
            }
            while (it_combi1 != combi1.end() || it_combi2 != combi2.end()) {
                if (it_combi1 == combi1.end()) {
                    combi_merged.push_back(*it_combi2++);
                    continue;
                }
                if (it_combi2 == combi2.end()) {
                    combi_merged.push_back(*it_combi1++);
                    continue;
                }
                vec2i split1 = *it_combi1;
                vec2i split2 = *it_combi2;
                if (compareSplits(split1, split2) > 0) {
                    combi_merged.push_back(split1);
                    it_combi1++;
                } else {
                    combi_merged.push_back(split2);
                    it_combi2++;
                }
            }
            result.push_back(combi_merged);
        }
    }
    return result;
}

std::list<combi> CombiVertexSplitter::generateCombinationsRec(
        combi history, std::vector<int> labels) {
    std::list<combi> result;
    std::list<vec2i> splits = createSingleSplitCombinations(labels);
    std::list<combi> combis = appendSplitCombinations(history, splits);
    if (labels.size() <= 4) {
        result = combis;
    } else {
        std::list<combi>::iterator it_combis = combis.begin();
        while (it_combis != combis.end()) {
            combi split_combi = *it_combis++;
            std::list<combi> combis_next;
            vec2i split = split_combi.back();
            std::vector<int> labels1(labels);
            std::vector<int> labels2 = splitLabels(labels1, split);
            if (labels1.size() > 3 && labels2.size() > 3) {
                std::list<combi> combis1 = generateCombinationsRec(split_combi, labels1);
                std::list<combi> combis2 = generateCombinationsRec(split_combi, labels2);
                combis_next = mergeCombinations(split_combi, combis1, combis2);
            } else if (labels1.size() > 3) {
                combis_next = generateCombinationsRec(split_combi, labels1);
            } else if (labels2.size() > 3) {
                combis_next = generateCombinationsRec(split_combi, labels2);
            }
            result.insert(result.end(), combis_next.begin(), combis_next.end());
        }
    }
    return result;
}

std::list<combi> CombiVertexSplitter::generateAllCombinations(unsigned int degree) {
    combi history;
    std::vector<int> labels = initLabels(degree);
    std::list<combi> result = generateCombinationsRec(history, labels);
    DEBUG_VAR(degree);
    std::list<combi>::iterator it_combis = result.begin();
    while (it_combis != result.end()) {
        combi combination = *it_combis++;
        DEBUG_VAL(combiToString(combination));
    }
    return result;
}

PolyhedronSPtr CombiVertexSplitter::copyVertex(VertexSPtr vertex) {
    PolyhedronSPtr result = Polyhedron::create();
    VertexSPtr vertex_c = vertex->clone();
    result->addVertex(vertex_c);
    FacetSPtr facet = vertex->firstFacet();
    EdgeSPtr edge_first = vertex->findEdge(facet);
    EdgeSPtr edge;
    EdgeSPtr edge_prev;
    EdgeSPtr edge_prev_c;
    EdgeSPtr edge_first_c;
    FacetSPtr facet_next = facet;
    FacetSPtr facet_c;
    SkelFacetDataSPtr facet_c_data;
    while (edge != edge_first) {
        if (!edge) {
            edge = edge_first;
        }
        EdgeSPtr edge_c;
        VertexSPtr vertex_dst_c;
        if (edge->getVertexSrc() == vertex) {
            vertex_dst_c = edge->getVertexDst()->clone();
        } else if (edge->getVertexDst() == vertex) {
            vertex_dst_c = edge->getVertexSrc()->clone();
        }
        result->addVertex(vertex_dst_c);
        edge_c = Edge::create(vertex_c, vertex_dst_c);
        if (edge == edge_first) {
            edge_first_c = edge_c;
        }
        result->addEdge(edge_c);
        if (edge_prev_c) {
            facet_c = Facet::create();
            facet_c->setPlane(facet->plane());
            edge_prev_c->setFacetL(facet_c);
            facet_c->addEdge(edge_prev_c);
            edge_c->setFacetR(facet_c);
            facet_c->addEdge(edge_c);
            facet_c_data = SkelFacetData::create(facet_c);
            facet_c_data->setFacetOrigin(facet);
            if (facet->hasData()) {
                facet_c_data->setSpeed(std::dynamic_pointer_cast<SkelFacetData>(
                        facet->getData())->getSpeed());
            }
            result->addFacet(facet_c);
        }
        edge_prev_c = edge_c;
        edge_prev = edge;
        edge = edge->next(vertex);
        facet = facet_next;
        facet_next = edge->other(facet);
    }
    facet_c = Facet::create();
    facet_c->setPlane(facet->plane());
    edge_prev_c->setFacetL(facet_c);
    facet_c->addEdge(edge_prev_c);
    edge_first_c->setFacetR(facet_c);
    facet_c->addEdge(edge_first_c);
    facet_c_data = SkelFacetData::create(facet_c);
    facet_c_data->setFacetOrigin(facet);
    if (facet->hasData()) {
        facet_c_data->setSpeed(std::dynamic_pointer_cast<SkelFacetData>(
                facet->getData())->getSpeed());
    }
    result->addFacet(facet_c);
    return result;
}

PolyhedronSPtr CombiVertexSplitter::splitVertex(VertexSPtr vertex, combi combination) {
    assert(vertex->degree() - 3 == combination.size());
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    std::list<VertexSPtr> vertices_tosplit;
    vertices_tosplit.push_back(vertex);

    std::vector<FacetSPtr> facets;
    FacetSPtr facet = vertex->firstFacet();
    EdgeSPtr edge_first = vertex->findEdge(facet);
    EdgeSPtr edge;
    while (edge != edge_first) {
        if (!edge) {
            edge = edge_first;
        }
        facets.push_back(facet);
        edge = edge->next(vertex);
        facet = edge->other(facet);
    }

    std::list<EdgeSPtr> edges_toremove;
    for (unsigned int i = 0; i < combination.size(); i++) {
        vec2i split = combination[i];
        FacetSPtr facet_right = facets[split[0]];
        FacetSPtr facet_left = facets[split[1]];
        std::list<VertexSPtr>::iterator it_v = vertices_tosplit.begin();
        while (it_v != vertices_tosplit.end()) {
            std::list<VertexSPtr>::iterator it_current = it_v;
            VertexSPtr vertex = *it_v++;
            if (vertex->containsFacet(facet_right) &&
                    vertex->containsFacet(facet_left)) {
                vertices_tosplit.erase(it_current);
                VertexSPtr vertex2 = vertex->split(facet_left, facet_right);
                if (facet_left->getPlane() == facet_right->getPlane()) {
                    EdgeSPtr edge = vertex->findEdge(vertex2);
                    edges_toremove.push_back(edge);
                }
                if (vertex->hasData()) {
                    SkelVertexDataSPtr data = std::dynamic_pointer_cast<SkelVertexData>(
                            vertex->getData());
                    SkelVertexDataSPtr data2 = SkelVertexData::create(vertex2);
                    data2->setNode(data->getNode());
                }
                if (vertex->degree() > 3) {
                    vertices_tosplit.push_back(vertex);
                }
                if (vertex2->degree() > 3) {
                    vertices_tosplit.push_back(vertex2);
                }
                break;
            }
        }
    }

    std::list<EdgeSPtr>::iterator it_e = edges_toremove.begin();
    while (it_e != edges_toremove.end()) {
        EdgeSPtr edge = *it_e++;
        VertexSPtr vertex_src = edge->getVertexSrc();
        VertexSPtr vertex_dst = edge->getVertexDst();
        FacetSPtr facet_l = edge->getFacetL();
        FacetSPtr facet_r = edge->getFacetR();
        facet_l->removeEdge(edge);
        facet_r->removeEdge(edge);
        polyhedron->removeEdge(edge);
        if (facet_l != facet_r) {
            facet_l->merge(facet_r);
            polyhedron->removeFacet(facet_r);
        }

        // remove vertices of degree 2
        for (unsigned int i = 0; i < 2; i++) {
            VertexSPtr vertex;
            if (i == 0) {
                vertex = vertex_src;
            } else if (i == 1) {
                vertex = vertex_dst;
            }
            VertexSPtr vertex_merged_src = vertex->prev(facet_l);
            VertexSPtr vertex_merged_dst = vertex->next(facet_l);
            FacetSPtr facet_r = facet_l->next(vertex);
            std::list<EdgeWPtr>::iterator it_ew = vertex->edges().begin();
            while (it_ew != vertex->edges().end()) {
                EdgeWPtr edge_wptr = *it_ew++;
                if (!edge_wptr.expired()) {
                    EdgeSPtr edge_toremove(edge_wptr);
                    facet_l->removeEdge(edge_toremove);
                    facet_r->removeEdge(edge_toremove);
                    polyhedron->removeEdge(edge_toremove);
                }
            }
            std::list<FacetWPtr>::iterator it_fw = vertex->facets().begin();
            while (it_fw != vertex->facets().end()) {
                FacetWPtr facet_wptr = *it_fw++;
                if (!facet_wptr.expired()) {
                    FacetSPtr facet(facet_wptr);
                    facet->removeVertex(vertex);
                }
            }
            polyhedron->removeVertex(vertex);
            EdgeSPtr edge_merged = Edge::create(vertex_merged_src, vertex_merged_dst);
            edge_merged->setFacetL(facet_l);
            edge_merged->setFacetR(facet_r);
            facet_l->addEdge(edge_merged);
            facet_r->addEdge(edge_merged);
            polyhedron->addEdge(edge_merged);
        }
    }
    return polyhedron;
}

PolyhedronSPtr CombiVertexSplitter::apply(PolyhedronSPtr poly_split, VertexSPtr vertex) {
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    NodeSPtr node;
    if (vertex->hasData()) {
        SkelVertexDataSPtr data =
                std::dynamic_pointer_cast<SkelVertexData>(vertex->getData());
        node = data->getNode();
    }
    std::map<VertexSPtr, VertexSPtr> vertices;
    std::list<VertexSPtr>::iterator it_v = poly_split->vertices().begin();
    while (it_v != poly_split->vertices().end()) {
        VertexSPtr vertex_ps = *it_v++;
        if (vertex_ps->degree() > 1) {
            VertexSPtr vertex_vs = Vertex::create(vertex_ps->getPoint());
            vertices[vertex_ps] = vertex_vs;
            polyhedron->addVertex(vertex_vs);
            SkelVertexDataSPtr data = SkelVertexData::create(vertex_vs);
            data->setNode(node);
        }
    }

    std::list<EdgeSPtr>::iterator it_e = poly_split->edges().begin();
    while (it_e != poly_split->edges().end()) {
        EdgeSPtr edge_ps = *it_e++;
        VertexSPtr vertex_ps_src = edge_ps->getVertexSrc();
        VertexSPtr vertex_ps_dst = edge_ps->getVertexDst();
        FacetSPtr facet_ps_l = edge_ps->getFacetL();
        FacetSPtr facet_ps_r = edge_ps->getFacetR();
        if (vertex_ps_dst->degree() == 1) {
            EdgeSPtr edge_vs;
            std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
            while (it_ve != vertex->edges().end()) {
                EdgeWPtr edge_wptr = *it_ve++;
                if (!edge_wptr.expired()) {
                    EdgeSPtr edge(edge_wptr);
                    if (edge->getVertexSrc()->getPoint() == vertex_ps_dst->getPoint() ||
                            edge->getVertexDst()->getPoint() == vertex_ps_dst->getPoint()) {
                        edge_vs = edge;
                        break;
                    }
                }
            }
            VertexSPtr vertex_vs = vertices[vertex_ps_src];
            if (edge_vs->getVertexSrc() == vertex) {
                edge_vs->replaceVertexSrc(vertex_vs);
            } else {
                edge_vs->replaceVertexDst(vertex_vs);
            }
            if (!edge_vs->getFacetL()->containsVertex(vertex_vs)) {
                edge_vs->getFacetL()->addVertex(vertex_vs);
            }
            if (!edge_vs->getFacetR()->containsVertex(vertex_vs)) {
                edge_vs->getFacetR()->addVertex(vertex_vs);
            }
        } else if (vertex_ps_src->degree() == 1) {
            EdgeSPtr edge_vs;
            std::list<EdgeWPtr>::iterator it_ve = vertex->edges().begin();
            while (it_ve != vertex->edges().end()) {
                EdgeWPtr edge_wptr = *it_ve++;
                if (!edge_wptr.expired()) {
                    EdgeSPtr edge(edge_wptr);
                    if (edge->getVertexSrc()->getPoint() == vertex_ps_src->getPoint() ||
                            edge->getVertexDst()->getPoint() == vertex_ps_src->getPoint()) {
                        edge_vs = edge;
                        break;
                    }
                }
            }
            VertexSPtr vertex_vs = vertices[vertex_ps_dst];
            if (edge_vs->getVertexSrc() == vertex) {
                edge_vs->replaceVertexSrc(vertex_vs);
            } else {
                edge_vs->replaceVertexDst(vertex_vs);
            }
            if (!edge_vs->getFacetL()->containsVertex(vertex_vs)) {
                edge_vs->getFacetL()->addVertex(vertex_vs);
            }
            if (!edge_vs->getFacetR()->containsVertex(vertex_vs)) {
                edge_vs->getFacetR()->addVertex(vertex_vs);
            }
        } else {
            EdgeSPtr edge_vs = Edge::create(vertices[vertex_ps_src], vertices[vertex_ps_dst]);
            SkelFacetDataSPtr data_l =
                    std::dynamic_pointer_cast<SkelFacetData>(
                    facet_ps_l->getData());
            SkelFacetDataSPtr data_r =
                    std::dynamic_pointer_cast<SkelFacetData>(
                    facet_ps_r->getData());
            FacetSPtr facet_vs_l = data_l->getFacetOrigin();
            FacetSPtr facet_vs_r = data_r->getFacetOrigin();
            edge_vs->setFacetL(facet_vs_l);
            edge_vs->setFacetR(facet_vs_r);
            facet_vs_l->addEdge(edge_vs);
            facet_vs_r->addEdge(edge_vs);
            polyhedron->addEdge(edge_vs);
        }
    }

    std::list<FacetWPtr>::iterator it_f = vertex->facets().begin();
    while (it_f != vertex->facets().end()) {
        FacetWPtr facet_wptr = *it_f++;
        if (!facet_wptr.expired()) {
            FacetSPtr facet(facet_wptr);
            facet->removeVertex(vertex);
        }
    }
    polyhedron->removeVertex(vertex);
    return polyhedron;
}

PolyhedronSPtr CombiVertexSplitter::splitVertex(VertexSPtr vertex) {
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    WriteLock l(polyhedron->mutex());
    vertex->sort();
    std::list<combi> combinations = generateAllCombinations(vertex->degree());
    std::list<combi> combinations_valid;
    std::list<PolyhedronSPtr> polys_split;
    std::list<combi>::iterator it_combi = combinations.begin();
    while (it_combi != combinations.end()) {
        combi combination = *it_combi++;
        PolyhedronSPtr poly_c = copyVertex(vertex);
        VertexSPtr vertex_c = poly_c->vertices().front();
        splitVertex(vertex_c, combination);
        if (checkSplitted(poly_c)) {
            DEBUG_VAL("Valid split-combination found: " << combiToString(combination));
            combinations_valid.push_back(combination);
            polys_split.push_back(poly_c);
        }
    }

    unsigned int selected_combi = selected_combi_ % combinations_valid.size();
    it_combi = combinations_valid.begin();
    std::list<PolyhedronSPtr>::iterator it_p = polys_split.begin();
    for (unsigned int i = 0; i < combinations_valid.size(); i++) {
        combi combination_valid = *it_combi++;
        PolyhedronSPtr poly_split = *it_p++;
        if (i == selected_combi) {
            DEBUG_VAL("Selected split-combination: " << combiToString(combination_valid));
            apply(poly_split, vertex);
            break;
        }
    }
    return polyhedron;
}

std::string CombiVertexSplitter::combiToString(combi combination) {
    std::stringstream sstr;
    sstr << "{ ";
    for (unsigned int i = 0; i < combination.size(); i++) {
        vec2i split = combination[i];
        if (i > 0) {
            sstr << ", ";
        }
        sstr << "{" << split[0] << ", " << split[1] << "}";
    }
    sstr << " }";
    return sstr.str();
}

std::string CombiVertexSplitter::toString() const {
    std::stringstream sstr;
    sstr << "CombiVertexSplitter(" << selected_combi_ << ")";
    return sstr.str();
}

} }
