// Copyright (c) 2019-2020 X, The Moonshot Factory (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later
//
//
// Author(s)     : Pierre Alliez pierre.alliez@inria.fr
//               : Michael Hemmer mhsaar@gmail.com
//               : Cedric Portaneri cportaneri@gmail.com
//
#ifndef CGAL_ALPHA_WRAP_2_DEMO_TYPES_H
#define CGAL_ALPHA_WRAP_2_DEMO_TYPES_H

#define CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
#define CGAL_AW2_COMPUTE_AND_STORE_STEINER_INFO_AT_GATE_CREATION

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/alpha_wrap_2.h>

#include "pslg_2.h"
#include "pslg_io.h"
#include "conversion_utils.h"

namespace AW2 = CGAL::Alpha_wraps_2;
namespace AW2i = CGAL::Alpha_wraps_2::internal;

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = EPICK::FT;
using Point_2 = EPICK::Point_2;
using Segment_2 = EPICK::Segment_2;
using Vector_2 = EPICK::Vector_2;
using Line_2 = EPICK::Line_2;
using Pslg_component = AW2i::Component<EPICK>;
using Pslg = AW2i::Pslg<EPICK>;

using Segment_Oracle = AW2i::Segment_soup_oracle<EPICK>;
using Oracle = AW2i::Point_set_oracle<EPICK, Segment_Oracle>;
using Alpha_wrapper = AW2i::Alpha_wrapper_2<Oracle>;

using Triangulation = Alpha_wrapper::Triangulation;
using Vertex_handle = Triangulation::Vertex_handle;
using Edge = Triangulation::Edge;
using Face_handle = Triangulation::Face_handle;

using Alpha_PQ = Alpha_wrapper::Alpha_PQ;
using Gate = Alpha_wrapper::Gate;


#endif  // CGAL_ALPHA_WRAP_2_DEMO_TYPES_H
