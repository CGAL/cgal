template <class Polygon_2>
void make_monotone_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

/*
   polygon.push_back(Point_2(186, 101));
   polygon.push_back(Point_2(304, 101));
   polygon.push_back(Point_2(351, 101));
   polygon.push_back(Point_2(367, 122));
   polygon.push_back(Point_2(388, 151));
   polygon.push_back(Point_2(407, 179));
   polygon.push_back(Point_2(339, 208));
   polygon.push_back(Point_2(280, 220));
   polygon.push_back(Point_2(214, 204));
   polygon.push_back(Point_2(214, 178));
   polygon.push_back(Point_2(214, 158));
   polygon.push_back(Point_2(164, 156));
*/   
/*
   polygon.push_back(Point_2(149,  91));
   polygon.push_back(Point_2(309,  91));
   polygon.push_back(Point_2(351,  91));
   polygon.push_back(Point_2(384,  91));
   polygon.push_back(Point_2(384, 127));
   polygon.push_back(Point_2(384, 158));
   polygon.push_back(Point_2(384, 186));
   polygon.push_back(Point_2(220, 186));
   polygon.push_back(Point_2(152, 186));
*/
   polygon.push_back(Point_2(140, 226));
   polygon.push_back(Point_2(293, 226));
   polygon.push_back(Point_2(335, 226));
   polygon.push_back(Point_2(358, 226));
   polygon.push_back(Point_2(358, 264));
   polygon.push_back(Point_2(358, 286));
   polygon.push_back(Point_2(358, 315));
   polygon.push_back(Point_2(358, 334));
   polygon.push_back(Point_2(297, 334));
   polygon.push_back(Point_2(246, 334));
   polygon.push_back(Point_2(211, 334));
   polygon.push_back(Point_2(163, 334));
   polygon.push_back(Point_2(146, 307));
   polygon.push_back(Point_2(132, 278));
   polygon.push_back(Point_2(132, 259));
   polygon.push_back(Point_2(132, 249));
   polygon.push_back(Point_2(132, 236));
}

template <class Polygon_2>
void make_monotone_convex(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(239, 108));
   polygon.push_back(Point_2(358, 170));
   polygon.push_back(Point_2(387, 305));
   polygon.push_back(Point_2(340, 416));
   polygon.push_back(Point_2(215, 432));
   polygon.push_back(Point_2(106, 358));
   polygon.push_back(Point_2(106, 236));
   polygon.push_back(Point_2(138, 162));
}

template <class Polygon_2>
void make_nonmonotone(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(336, 218));
   polygon.push_back(Point_2(444, 290));
   polygon.push_back(Point_2(416, 402));
   polygon.push_back(Point_2(282, 331));
   polygon.push_back(Point_2(366, 290));
   polygon.push_back(Point_2(314, 245));
   polygon.push_back(Point_2(209, 269));
   polygon.push_back(Point_2(202, 371));
   polygon.push_back(Point_2( 75, 261));
   polygon.push_back(Point_2(148, 267));
   polygon.push_back(Point_2(166, 231));
   polygon.push_back(Point_2(176, 256));
   polygon.push_back(Point_2(215, 152));
   polygon.push_back(Point_2(145, 134));
   polygon.push_back(Point_2(286, 123));
   polygon.push_back(Point_2(226, 197));
   polygon.push_back(Point_2(274, 208));
   polygon.push_back(Point_2(375, 125));
}
