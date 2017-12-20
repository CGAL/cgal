// Copyright (c) 2013-2015  The University of Western Sydney, Australia.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s): Quincy Tse, Weisheng Si

/*! \file gnuplot_output_2.h
 *
 * This header implements the function that can generate data and script files for plotting
 * graphs by Gnuplot. This function requires that graphs be represented by boost::adjacency_list
 * with its template parameter VertexProperties set to CGAL::Point_2.
 */

#ifndef GNUPLOT_OUTPUT_2_H
#define GNUPLOT_OUTPUT_2_H

#include <CGAL/license/Cone_spanners_2.h>


#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

/*  ------ Declarations go first, then implementations follow.   ------ */

/*!
*  \ingroup PkgConeBasedSpanners
*  \brief Output a set of files used by Gnuplot to plot `g`.
*
*  The files that are generated for Gnuplot are:
*  (1) prefix.v (vertex list)
*  (2) prefix.plt (Gnuplot script), This script will read
*      prefix.v as input to plot the vertex list. The edge list is also
*      included in this script.
*
*  Notes:
*  (1) If these files already exists, this function will overwrite these
*      files.
*  (2) Parallel and self-edges cannot be plotted.
*
*  \tparam Graph  The type of the graph to be plotted. For this function to work,
*                 the graph type must be `boost::adjacency_list` with `CGAL::Point_2`
*                 as the `VertexProperties`.
*
*  \param g       A `boost::adjacency_list` graph with `CGAL::Point_2` as the VertexProperties to be plotted
*  \param prefix  The prefix of the output files names
*/
template <typename Graph>
void gnuplot_output_2 (const Graph& g, const std::string& prefix);

/*
*  \brief Compiles a multi-lined %string to draw the edges in \p g by Gnuplot.
*  Compiles an edge list in the following format:
*
*  set arrow from (start x, start y) to (end x, end y)
*   ...
*
*  NOTE: For undirected graphs, use "set style arrow nohead"; for directed graphs,
*        use "set style arrow head"
*
*  \param g  A boost::adjacency_list graph with CGAL::Point_2 as the VertexProperties to be plotted
*  \return   The edge list string.
*/
template <typename Graph>
std::string gnuplot_edge_list (const Graph& g);

/*
*  \brief Compiles a multi-lined %string representing the vertices in \p g.
*
* Compiles a vertex list in the following format:
*  x  y
*  x  y
*  ...
*
*  \param g  A boost::adjacency_list graph with CGAL::Point_2 as the VertexProperties to be plotted
*  \return   The vertex list string.
*/
template <typename Graph>
std::string gnuplot_vertex_list (const Graph& g);

/* This struct is defined to use partial specialization to generate arrow styles differently for
 * directed and undirected graphs.
 * Note: Need to use structs because C++ before 11 doesn't allow partial specialisation
 * for functions
 */
template <typename Graph, typename Directedness=typename Graph::directed_selector>
struct Gnuplot_edges_2;

/* -------   IMPLEMENTATIONS  ------- */

template <typename Graph>
std::string gnuplot_edge_list (const Graph& g)
{
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;  // Use fixed floating-point notation

    typename Graph::edge_iterator eit, ee;
    for (boost::tie(eit, ee) = boost::edges(g); eit != ee; ++eit) {
        typename Graph::vertex_descriptor src = boost::source(*eit, g);
        typename Graph::vertex_descriptor end = boost::target(*eit, g);
        ss << "set arrow from ";
        ss << to_double(g[src].x()) << "," << to_double(g[src].y());
        ss << " to ";
        ss << to_double(g[end].x()) << "," << to_double(g[end].y());
        ss << " as 1" << std::endl;
    }
    return ss.str();
}

// Common regardless of whether g is directed.
template <typename Graph>
std::string gnuplot_vertex_list(const Graph& g) {
    std::stringstream ss;
    ss.precision(3);
    ss << std::fixed;

    typename Graph::vertex_iterator vit, ve;
    for (boost::tie(vit, ve) = boost::vertices(g); vit != ve; ++vit) {
        ss << to_double(g[*vit].x()) << "  " << to_double(g[*vit].y()) << std::endl;
    }
    return ss.str();
}

template <typename Graph>
void gnuplot_output_2 (const Graph& g, const std::string& prefix)
{
    // Generate the vertex list to the .v file
    std::ofstream fs((prefix + ".v").c_str(),
                     std::ofstream::out | std::ofstream::trunc);
    fs << gnuplot_vertex_list(g);
    fs.close();
    std::cout << prefix << ".v" << " is generated. " << std::endl;

    // Generate the Gnuplot Script to the .plt file
    fs.open((prefix + ".plt").c_str(),
            std::ofstream::out | std::ofstream::trunc);
    fs << "set term ";
    //  Choose one:
    fs << "wxt ";
    // fs << "postscript eps ";

    fs << "font \", 9\" enhanced" << std::endl;

    //  Uncomment if eps:
    // fs << "set output \"" << prefix << ".eps\"" << std::endl;

    fs << "set title" << std::endl;
    fs << "set xlabel  # when no options, clear the xlabel" << std::endl;
    fs << "set ylabel" << std::endl;
    fs << "unset key" << std::endl;
    fs << "set size square" << std::endl;
    fs << "unset xtics" << std::endl;
    fs << "unset ytics" << std::endl;
    fs << "unset border" << std::endl;

    /* Uncomment if you need the following:
    ss << "set xtics" << std::endl;
    ss << "set ytics" << std::endl;
    ss << "set border" << std::endl;
    ss << "set grid xtics ytics" << std::endl;
    */

    fs << std::endl;
    // Specific part that depends on whether g is directed
    fs << Gnuplot_edges_2<Graph>::gnuplot_edge_script(g);
    fs << std::endl;

    // plot the vertices
    fs << "plot \"" << prefix << ".v\" with points pt 7 ps 0.8 lt rgb \"blue\"" << std::endl;

    //  Uncomment if wxt and also want eps output:
    //  fs << "set term postscript eps " << std::endl;
    //  fs << "set output \"" << prefix << ".eps\"" << std::endl;
    //  fs << "replot" << std::endl;

    fs.close();
    std::cout << prefix << ".plt" << " is generated. " << std::endl;
}

// directed graphs
/* Writing edge list to the gnuplot script for directed graphs */
template <typename Graph>
struct Gnuplot_edges_2<Graph, boost::directedS> {

    // Uses "set style arrow 1 head" to draw directed edges
    // Edges are written directly into the script file.
    static std::string gnuplot_edge_script(const Graph& g)
    {
        std::stringstream ss;

        ss << "set style arrow 1 head filled lc rgb \"black\"" << std::endl;
        ss << std::endl;
        ss << "# edges" << std::endl;
        ss << gnuplot_edge_list(g);
        ss << "# end of edges" << std::endl;

        return ss.str();
    }
};

// Bidirectional graph, the same as the directed graph.
/* Writing edge list to the gnuplot script for bidirectional graphs */
template <typename Graph>
struct Gnuplot_edges_2<Graph, boost::bidirectionalS> : public Gnuplot_edges_2<Graph, boost::directedS> {};

// Undirected graphs
/* Writing edge list to the gnuplot script for undirected graphs */
template <typename Graph>
struct Gnuplot_edges_2<Graph, boost::undirectedS> {

    static std::string gnuplot_edge_script(const Graph& g)
    {
        std::stringstream ss;
        ss << "set style arrow 1 nohead lc rgb \"black\"" << std::endl;
        ss << std::endl;
        ss << "# edges" << std::endl;
        ss << gnuplot_edge_list(g);
        ss << "# end of edges" << std::endl;

        return ss.str();
    }

};

}  // namespace CGAL

#endif
