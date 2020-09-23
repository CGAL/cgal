#include "Kernel_type.h"
#include "SMesh_type.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polygon_soup_item.h"
#include "Point_set_3.h"

#include <CGAL/compute_average_spacing.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Alpha_shape_mesher.h>
#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_smoother.h>


// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

typedef CGAL::Scale_space_surface_reconstruction_3<Kernel> ScaleSpace;
typedef CGAL::Scale_space_reconstruction_3::Advancing_front_mesher<Kernel> ScaleSpaceAFM;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<Kernel> ScaleSpaceASM;
typedef CGAL::Scale_space_reconstruction_3::Jet_smoother<Kernel> ScaleSpaceJS;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<Kernel> ScaleSpaceWPS;

void scale_space (const Point_set& points,
                  std::vector<Scene_polygon_soup_item*>& items,
                  bool jet_smoother,
                  unsigned int iterations,
                  unsigned int neighbors, unsigned int fitting, unsigned int monge,
                  unsigned int neighborhood_size, unsigned int samples,
                  bool advancing_front_mesher,
                  bool generate_smooth,
                  double longest_edge, double radius_ratio_bound, double beta_angle,
                  bool separate_shells, bool force_manifold)
{
  ScaleSpace reconstruct (points.points().begin(), points.points().end());

  double squared_radius = 0.;
  if (jet_smoother)
  {
    ScaleSpaceJS smoother(neighbors, fitting, monge);
    reconstruct.increase_scale(iterations, smoother);
    if (!advancing_front_mesher)
      squared_radius = CGAL::compute_average_spacing<Concurrency_tag> (points, neighbors);
  }
  else
  {
    ScaleSpaceWPS smoother(neighborhood_size, samples);
    reconstruct.increase_scale(iterations, smoother);
    squared_radius = smoother.squared_radius();

  }

  if (advancing_front_mesher)
  {
    ScaleSpaceAFM mesher (longest_edge, radius_ratio_bound, beta_angle);
    reconstruct.reconstruct_surface (mesher);

    Scene_polygon_soup_item* new_item
      = new Scene_polygon_soup_item ();
    new_item->setColor(Qt::lightGray);
    new_item->setRenderingMode(FlatPlusEdges);
    new_item->init_polygon_soup(points.size(), reconstruct.number_of_facets ());

    Scene_polygon_soup_item* smooth_item = NULL;
    if (generate_smooth)
    {
      smooth_item = new Scene_polygon_soup_item ();
      smooth_item->setColor(Qt::lightGray);
      smooth_item->setRenderingMode(FlatPlusEdges);
      smooth_item->init_polygon_soup(points.size(), reconstruct.number_of_facets ());
    }

    std::map<std::size_t, std::size_t> map_i2i;
    std::size_t current_index = 0;

    for (ScaleSpace::Facet_iterator it = reconstruct.facets_begin();
         it != reconstruct.facets_end(); ++ it)
    {
      for (unsigned int ind = 0; ind < 3; ++ ind)
      {
        if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
        {
          map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
          Point_3 p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
          new_item->new_vertex (p.x (), p.y (), p.z ());

          if (generate_smooth)
          {
            p = *(reconstruct.points_begin() + (*it)[ind]);
            smooth_item->new_vertex (p.x (), p.y (), p.z ());
          }
        }
      }
      new_item->new_triangle( map_i2i[(*it)[0]],
                              map_i2i[(*it)[1]],
                              map_i2i[(*it)[2]] );
      if (generate_smooth)
        smooth_item->new_triangle( map_i2i[(*it)[0]],
                                   map_i2i[(*it)[1]],
                                   map_i2i[(*it)[2]] );

    }

    items.push_back(new_item);
    if (generate_smooth)
      items.push_back (smooth_item);
  }
  else
  {
    ScaleSpaceASM mesher (squared_radius, separate_shells, force_manifold);
    reconstruct.reconstruct_surface (mesher);

    for( unsigned int sh = 0; sh < mesher.number_of_shells(); ++sh )
    {
      Scene_polygon_soup_item* new_item
        = new Scene_polygon_soup_item ();
      new_item->setColor(Qt::lightGray);
      new_item->setRenderingMode(FlatPlusEdges);
      new_item->init_polygon_soup(points.size(), mesher.number_of_triangles ());

      Scene_polygon_soup_item* smooth_item = NULL;
      if (generate_smooth)
      {
        smooth_item = new Scene_polygon_soup_item ();
        smooth_item->setColor(Qt::lightGray);
        smooth_item->setRenderingMode(FlatPlusEdges);
        smooth_item->init_polygon_soup(points.size(), mesher.number_of_triangles ());
      }

      std::map<unsigned int, unsigned int> map_i2i;
      unsigned int current_index = 0;

      for (ScaleSpaceASM::Facet_iterator it = mesher.shell_begin (sh);
           it != mesher.shell_end (sh); ++ it)
      {
        for (unsigned int ind = 0; ind < 3; ++ ind)
        {
          if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
          {
            map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
            Point_3 p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
            new_item->new_vertex (p.x (), p.y (), p.z ());

            if (generate_smooth)
            {
              p = *(reconstruct.points_begin() + (*it)[ind]);
              smooth_item->new_vertex (p.x (), p.y (), p.z ());
            }
          }
        }
        new_item->new_triangle( map_i2i[(*it)[0]],
                                map_i2i[(*it)[1]],
                                map_i2i[(*it)[2]] );
        if (generate_smooth)
          smooth_item->new_triangle( map_i2i[(*it)[0]],
                                     map_i2i[(*it)[1]],
                                     map_i2i[(*it)[2]] );

      }

      items.push_back (new_item);
      if (generate_smooth)
        items.push_back (smooth_item);
    }

    if (force_manifold)
    {
      std::ptrdiff_t num = std::distance( mesher.garbage_begin(  ),
                                          mesher.garbage_end(  ) );

      Scene_polygon_soup_item* new_item
        = new Scene_polygon_soup_item ();
      new_item->setColor(Qt::blue);
      new_item->setRenderingMode(FlatPlusEdges);
      new_item->init_polygon_soup(points.size(), num);

      Scene_polygon_soup_item* smooth_item = NULL;
      if (generate_smooth)
      {
        smooth_item = new Scene_polygon_soup_item ();
        smooth_item->setColor(Qt::blue);
        smooth_item->setRenderingMode(FlatPlusEdges);
        smooth_item->init_polygon_soup(points.size(), num);
      }

      std::map<std::size_t, std::size_t> map_i2i;

      std::size_t current_index = 0;
      for (ScaleSpaceASM::Facet_iterator it=mesher.garbage_begin(),
             end=mesher.garbage_end();it!=end;++it)
      {
        for (unsigned int ind = 0; ind < 3; ++ ind)
        {
          if (map_i2i.find ((*it)[ind]) == map_i2i.end ())
          {
            map_i2i.insert (std::make_pair ((*it)[ind], current_index ++));
            Point_3 p = points.point(*(points.begin_or_selection_begin() + (*it)[ind]));
            new_item->new_vertex (p.x (), p.y (), p.z ());

            if (generate_smooth)
            {
              p = *(reconstruct.points_begin() + (*it)[ind]);
              smooth_item->new_vertex (p.x (), p.y (), p.z ());
            }
          }

        }
        new_item->new_triangle( map_i2i[(*it)[0]],
                                map_i2i[(*it)[1]],
                                map_i2i[(*it)[2]] );
        if (generate_smooth)
          smooth_item->new_triangle( map_i2i[(*it)[0]],
                                     map_i2i[(*it)[1]],
                                     map_i2i[(*it)[2]] );
      }

      items.push_back (new_item);
      if (generate_smooth)
        items.push_back (smooth_item);

    }
  }
}

