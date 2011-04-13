template <class Polygon_2>
void make_convex_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

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
void make_nonconvex_w_collinear_points(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(129, 147));
   polygon.push_back(Point_2(232, 147));
   polygon.push_back(Point_2(344, 147));
   polygon.push_back(Point_2(365, 234));
   polygon.push_back(Point_2(239, 251));
   polygon.push_back(Point_2(249, 210));
   polygon.push_back(Point_2(190, 212));
   polygon.push_back(Point_2(205, 323));
   polygon.push_back(Point_2(104, 323));
   polygon.push_back(Point_2(104, 279));
   polygon.push_back(Point_2(104, 243));
}

template <class Polygon_2>
void make_nonconvex(Polygon_2& polygon)
{
   typedef typename Polygon_2::Point_2  Point_2;

   polygon.push_back(Point_2(391, 374));
   polygon.push_back(Point_2(240, 431));
   polygon.push_back(Point_2(252, 340));
   polygon.push_back(Point_2(374, 320));
   polygon.push_back(Point_2(289, 214));
   polygon.push_back(Point_2(134, 390));
   polygon.push_back(Point_2( 68, 186));
   polygon.push_back(Point_2(154, 259));
   polygon.push_back(Point_2(161, 107));
   polygon.push_back(Point_2(435, 108));
   polygon.push_back(Point_2(208, 148));
   polygon.push_back(Point_2(295, 160));
   polygon.push_back(Point_2(421, 212));
   polygon.push_back(Point_2(441, 303));
}
