// Copyright (c) 2015,2016 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_H
#define CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Image_3.h>
#include <cstdlib>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <boost/container/deque.hpp>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/cstdint.hpp>
#ifdef CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
#  include <iostream>
#  include <boost/format.hpp>
#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
template <typename PointsOutputIterator,
          typename DomainsOutputIterator,
          typename TransformOperator,
          typename Construct_point,
          typename Image_word_type>
void
search_for_connected_components_in_labeled_image(const CGAL::Image_3& image,
                                                 PointsOutputIterator it,
                                                 DomainsOutputIterator dom_it,
                                                 TransformOperator transform,
                                                 Construct_point point,
                                                 Image_word_type)
{
  const std::size_t nx = image.xdim();
  const std::size_t ny = image.ydim();
  const std::size_t nz = image.zdim();
  const std::size_t size = nx * ny * nz;

  typedef boost::uint16_t uint;

  if(nx > 65535 || ny > 65535 || nz > 65535)
  {
    CGAL_error_msg("The dimensions of the image must be lower than 2^16");
  }

  typedef typename TransformOperator::result_type Label;

  std::vector<bool> visited(size, false);
  std::vector<bool> second_pass(size, false);
  typedef boost::tuple<uint, uint, uint, uint> Indices;
  typedef std::queue<Indices, boost::container::deque<Indices> > Indices_queue;
  typedef std::vector<Indices>  Border_vector;

#ifdef CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
  int number_of_connected_components = 0;
#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
  std::size_t voxel_index = 0;
  for(uint k=0; k<nz; k++)
    for(uint j=0; j<ny; j++)
      for(uint i=0; i<nx; i++)
      {
        using CGAL::IMAGEIO::static_evaluate;

        if(visited[voxel_index] | second_pass[voxel_index]) {
          ++voxel_index;
          continue;
        }
        const Label current_label =
          transform(static_evaluate<Image_word_type>(image.image(),
                                                     voxel_index));
        *dom_it++ = current_label;
        if(current_label == Label()) {
          visited[voxel_index] = true;
          second_pass[voxel_index] = true;
          ++voxel_index;
          continue;
        }

#ifdef CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
        // if we reach here, (i, j, k) is a new connected component
        ++number_of_connected_components;
        std::cerr << boost::format("Found new connected component (#%5%) "
                                   "at voxel (%1%, %2%, %3%), value=%4%, volume id=%6%\n")
          % i % j % k
          % (long)static_evaluate<Image_word_type>(image.image(), i, j, k)
          % number_of_connected_components
          % (int)current_label;
#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE

        int nb_voxels = 0;

        Indices_queue queue;
        Indices indices(i, j ,k, 0);
        queue.push(indices);

        Border_vector border;

        /*
         * First pass is a BFS to retrieve all the connected component, and
         * its border.
         * Second pass is a BFS initialized with all voxel of the border.
         * The last voxel of that BFS is used as the seed.
         */
        int pass = 1; // pass will be equal to 2 in second pass

        Indices bbox_min = indices;
        Indices bbox_max = indices;

        while(!queue.empty()) // walk through the connected component
        {
          Indices indices = queue.front();
          queue.pop();

          // warning: those indices i, j and k are local to the while loop
          const uint i = boost::get<0>(indices);
          const uint j = boost::get<1>(indices);
          const uint k = boost::get<2>(indices);
          const uint depth = boost::get<3>(indices);

          const size_t offset = i + nx * (j + ny * k);
          const int m = (visited[offset] ? 1 : 0) + (second_pass[offset] ? 2 : 0);
          if(m < pass)
          {
            if(pass == 1 )
            {
              visited[offset] = true;
              second_pass[offset] = false;
              ++nb_voxels;
              boost::get<0>(bbox_min) = (std::min)(i, boost::get<0>(bbox_min));
              boost::get<0>(bbox_max) = (std::max)(i, boost::get<0>(bbox_max));
              boost::get<1>(bbox_min) = (std::min)(j, boost::get<1>(bbox_min));
              boost::get<1>(bbox_max) = (std::max)(j, boost::get<1>(bbox_max));
              boost::get<2>(bbox_min) = (std::min)(k, boost::get<2>(bbox_min));
              boost::get<2>(bbox_max) = (std::max)(k, boost::get<2>(bbox_max));
            } else
            {
              CGAL_assertion(pass == 2);
              visited[offset] = false;
              second_pass[offset] = true;
            }

            static const int neighbors_offset[6][3] = { { +1,  0,  0 },
                                                        { -1,  0,  0 },
                                                        {  0, +1,  0 },
                                                        {  0, -1,  0 },
                                                        {  0,  0, +1 },
                                                        {  0,  0, -1 } };
            bool voxel_is_on_border = false;

            // Visit neighbors.
            // (i_n, j_n, k_n) are indices of neighbors.
            for(int n = 0; n < 6; ++n)
            {
              const ptrdiff_t i_n = i + neighbors_offset[n][0];
              const ptrdiff_t j_n = j + neighbors_offset[n][1];
              const ptrdiff_t k_n = k + neighbors_offset[n][2];
              if(i_n < 0 || i_n >= static_cast<ptrdiff_t>(nx) ||
                 j_n < 0 || j_n >= static_cast<ptrdiff_t>(ny) ||
                 k_n < 0 || k_n >= static_cast<ptrdiff_t>(nz))
              {
                voxel_is_on_border = true;
                continue;
              }
              else
              {
                const std::size_t offset_n = i_n + nx * (j_n + k_n * ny);
                if(transform(static_evaluate<Image_word_type>(image.image(),
                                                              offset_n))
                   == current_label)
                {
                  const int m_n = (visited[offset_n] ? 1 : 0) +
                    (second_pass[offset_n] ? 2 : 0);
                  if(m_n < pass) {
                    Indices indices(uint(i_n), uint(j_n), uint(k_n), uint(depth+1));
                    queue.push(indices);
                  }
                }
                else
                  voxel_is_on_border = true;
              }
            } // end for neighbors

            if(pass == 1 && voxel_is_on_border)
              border.push_back(indices);
          } // end if voxel not already visited

          if(queue.empty()) {
            if(pass == 1)
            { // End of first pass. Begin second pass with the voxels of
              // the border.
              for(typename Border_vector::const_iterator
                    border_it = border.begin(), border_end = border.end();
                  border_it != border_end; ++border_it)
                queue.push(*border_it);
              pass = 2;
            }
            else // end of second pass, return the last visited voxel
            {
//               if(nb_voxels >= 100)
              {
                *it++ = std::make_pair(point(i, j, k),
                                       depth+1);
#if CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE > 1
                std::cerr << boost::format("Found seed %5%, which is voxel "
                                           "(%1%, %2%, %3%), value=%4%\n")
                  % i % j % k
                  % (long)static_evaluate<Image_word_type>(image.image(), i, j, k)
                  % point(i, j, k);
#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE>1
              }
            }
          } // end if queue.empty()
        } // end while !queue.empty() (with local indices i, j, k)
#ifdef CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE
        std::cerr
          << boost::format("There was %1% voxel(s) in that component.\n"
                           "The bounding box is (%2% %3% %4%, %5% %6% %7%).\n"
                           "%8% voxel(s) on border\n")
          % nb_voxels
          % boost::get<0>(bbox_min) % boost::get<1>(bbox_min)
          % boost::get<2>(bbox_min)
          % boost::get<0>(bbox_max) % boost::get<1>(bbox_max)
          % boost::get<2>(bbox_max)
          % border.size();
#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_VERBOSE

        ++voxel_index;
      } // end for i,j,k
} // end function search_for_connected_components_in_labeled_image()

#endif // CGAL_MESH_3_SEARCH_FOR_CONNECTED_COMPONENTS_IN_LABELED_IMAGE_H
