// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   algo/3d/CombiVertexSplitter.h
 * @author Gernot Walzl
 * @date   2012-07-26
 */

#ifndef ALGO_3D_COMBIVERTEXSPLITTER_H
#define ALGO_3D_COMBIVERTEXSPLITTER_H

#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractVertexSplitter.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <boost/shared_array.hpp>
#include <list>
#include <vector>
#include <string>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

typedef boost::shared_array<int> vec2i;
typedef std::vector<vec2i> combi;

class CombiVertexSplitter : public AbstractVertexSplitter {
public:
    virtual ~CombiVertexSplitter();

    static CombiVertexSplitterSPtr create();

    static vec2i createSplit(int begin, int end);
    static int compareSplits(vec2i split1, vec2i split2);
    static std::vector<int> initLabels(unsigned int degree);
    static std::vector<int> splitLabels(std::vector<int>& labels, const vec2i& split);
    static std::list<vec2i> createSingleSplitCombinations(const std::vector<int>& labels);
    static std::list<combi> appendSplitCombinations(
            const combi& history, const std::list<vec2i>& splits);
    static std::list<combi> mergeCombinations(
            const combi& history,
            const std::list<combi>& combis1,
            const std::list<combi>& combis2);
    static std::list<combi> generateCombinationsRec(
            const combi& history, const std::vector<int>& labels);
    static std::list<combi> generateAllCombinations(unsigned int degree);

    static PolyhedronSPtr copyVertex(VertexSPtr vertex);
    static PolyhedronSPtr splitVertex(VertexSPtr vertex, const combi& combination);
    static PolyhedronSPtr apply(PolyhedronSPtr poly_split, VertexSPtr vertex);

    virtual PolyhedronSPtr splitVertex(VertexSPtr vertex);

    static std::string combiToString(const combi& combination);

    virtual std::string toString() const;

protected:
    CombiVertexSplitter();

    unsigned int selected_combi_;
};

} }

#endif /* ALGO_3D_COMBIVERTEXSPLITTER_H */
