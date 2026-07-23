
// Verify that every piece in the partition is simple and convex.
static bool partition_is_valid(const std::list<Polygon_2>& pieces)
{
   for (const Polygon_2& p : pieces)
   {
      if (!p.is_simple())   return false;
      if (!CGAL::is_convex_2(p.vertices_begin(), p.vertices_end()))
         return false;
   }
   return true;
}

void test_optimal_convex()
{
   Polygon_2              polygon;
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   assert(partition_polys.size() == 1 &&
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_convex_2(partition_polys.front().vertices_begin(),
                            partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_convex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   // Regression test for https://github.com/CGAL/cgal/issues/9322:
   // optimal_convex_partition_2 produced non-simple / non-convex pieces for
   // a cup-shaped octagon when certain vertices were used as the start vertex.
   // Rotations 0 (start at (4,4)), 1 (start at (6,6)), and 2 (start at (0,6))
   // were the originally failing cases.  Test all 8 rotations to be thorough.
   for (int rot = 0; rot < 8; ++rot)
   {
      partition_polys.clear();
      polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
      make_cup_octagon_issue_9322(polygon, rot);
      CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                       polygon.vertices_end(),
                                       std::back_inserter(partition_polys));
      assert(partition_is_valid(partition_polys));
   }
/*
   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_hilbert_polygon(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
*/
}
