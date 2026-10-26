#include <vector>

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
// The pieces for every choice of first vertex: always valid, and always the
// optimal number (issue #4398 depended on the first vertex).
void check_optimal_convex_count(const std::vector<Point_2>& points,
                                std::size_t expected)
{
   for (std::size_t first = 0; first < points.size(); ++first)
   {
      Polygon_2 polygon;
      for (std::size_t i = 0; i < points.size(); ++i)
         polygon.push_back(points[(first + i) % points.size()]);
      std::list<Polygon_2> partition_polys;
      CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                       polygon.vertices_end(),
                                       std::back_inserter(partition_polys));
      assert(partition_polys.size() == expected);
      assert(CGAL::convex_partition_is_valid_2(polygon.vertices_begin(),
                                               polygon.vertices_end(),
                                               partition_polys.begin(),
                                               partition_polys.end()));
   }
}

void test_optimal_convex_issue_4398()
{
   // The polygon of the issue: 3 pieces, but 4 when starting at (3,4).
   const double issue[][2] = {{3,4}, {1,5}, {0,4}, {2,2}, {5,0},
                              {4,3}, {6,0}, {6,6}, {5,7}, {4,5}};
   std::vector<Point_2> points;
   for (const auto& p : issue) points.push_back(Point_2(p[0], p[1]));
   check_optimal_convex_count(points, 3);

   // No three vertices collinear: 16 pieces, but 17 when starting at the
   // first vertex below.
   const double random41[][2] = {
      {-1076,-1296}, {-829,-1344}, {-818,-1378}, {-921,-1421}, {-895,-1459},
      {-743,-1575}, {-1070,-1463}, {-760,-1609}, {-588,-1769}, {-595,-1743},
      {-651,-1605}, {-586,-1529}, {-538,-1569}, {-447,-1626}, {-291,-1538},
      {-317,-1486}, {-268,-1377}, {-179,-1315}, {-312,-920}, {-426,-1091},
      {-502,-1162}, {-447,-930}, {-365,-997}, {-498,-876}, {-493,-885},
      {-545,-1280}, {-440,-1348}, {-365,-1462}, {-552,-1530}, {-648,-1361},
      {-528,-851}, {-583,-877}, {-685,-1341}, {-605,-913}, {-634,-809},
      {-645,-858}, {-810,-1078}, {-678,-844}, {-874,-994}, {-1060,-1231},
      {-1010,-1245}};
   points.clear();
   for (const auto& p : random41) points.push_back(Point_2(p[0], p[1]));
   check_optimal_convex_count(points, 16);
}

void test_optimal_convex()
{
   test_optimal_convex_issue_4398();

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
