#ifndef TEST_ARC_TRAITS_H
#define TEST_ARC_TRAITS_H

#if TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read an arc point */
template < typename T_Traits , typename stream >
bool read_arc_point(stream & is, typename T_Traits::Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = typename T_Traits::Point_2(lp);
  return true;
}

bool is_deg_1(char c)
{
  return (c=='z' || c=='Z') || (c=='y' || c=='Y') || (c=='x' || c=='X') ||
         (c=='w' || c=='W') || (c=='v' || c=='V') || (c=='l' || c=='L');
}

bool is_deg_2(char c)
{
  return (c=='b' || c=='B') || (c=='c' || c=='C') ||
         (c=='d' || c=='D') || (c=='e' || c=='E');
}

#if TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
template <class stream>
Circular_kernel::Line_arc_2 read_line(char type,stream & is)
{
  if (type == 'z' || type == 'Z')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circle_2 c_temp1,c_temp2;
    bool b1,b2;
    is >> l_temp >> c_temp1 >> b1 >> c_temp2 >> b2;
    return Circular_kernel::Line_arc_2(l_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'y' || type == 'Y')
  {
    Circular_kernel::Line_2 l_temp,l_temp1,l_temp2;
    is >> l_temp >> l_temp1 >> l_temp2;
    return Circular_kernel::Line_arc_2(l_temp,l_temp1,l_temp2);
  }
  else if (type == 'x' || type == 'X')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> l_temp >> p0 >> p1;
    //std::cout << "got here l_temp p0 p1 " << l_temp << " " << p0 << " " << p1 << std::endl;
    return Circular_kernel::Line_arc_2(l_temp,p0,p1);
  }
  else if (type == 'w' || type == 'W' || type == 'l' || type == 'L')
  {
    Circular_kernel::Point_2 p0,p1;
    is >> p0 >> p1;
    return Circular_kernel::Line_arc_2(p0,p1);
  }
  else if (type == 'v' || type == 'V')
  {
    Circular_kernel::Segment_2 seg;
    is >> seg;
    return Circular_kernel::Line_arc_2(seg);
  }
  std::cout << "should never happen Line_arc_2 " << type <<std::endl;
  return Circular_kernel::Line_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
template <class stream>
Circular_kernel::Circular_arc_2 read_arc(char type,stream & is)
{
  if (type == 'b' || type == 'B')
  {
    Circular_kernel::Circle_2 c_temp,c_temp1,c_temp2;
    bool b1,b2;
    is >> c_temp >> c_temp1 >> b1 >> c_temp2 >> b2 ;
    return Circular_kernel::Circular_arc_2(c_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'c' || type == 'C')
  {
    Circular_kernel::Circle_2 c_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> c_temp >> p0 >> p1;
    return Circular_kernel::Circular_arc_2(c_temp,p0,p1);
  }
  else if (type == 'd' || type == 'D')
  {
    Circular_kernel::Circle_2 c_temp;
    Circular_kernel::Line_2 l_temp1,l_temp2;
    bool b1,b2;
    is >> c_temp >> l_temp1 >> b1 >> l_temp2 >> b2;
    return Circular_kernel::Circular_arc_2(c_temp,l_temp1,b1,l_temp2,b2);
  }
  else if (type == 'e' || type == 'E')
  {
    Circular_kernel::Circular_arc_2 a_temp;
    Circular_kernel::Circle_2 c_temp;
    bool b1,b2;
    is >> a_temp >> b1 >> c_temp >> b2;
    return Circular_kernel::Circular_arc_2(a_temp,b1,c_temp,b2);
  }
  std::cout << "should never happen Circular_arc_2" << std::endl;
  return Circular_kernel::Circular_arc_2(); //should never happen
}
#endif

#endif

#if TEST_TRAITS == LINE_ARC_TRAITS

/*! Read a line arc point */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_point(stream & is, Traits::Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_xcurve(stream & is, Traits::X_monotone_curve_2 & xc)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    xc=read_line(type,is);
    return true;
  }
  return false;

}

/*! Read a general line arc curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is, Traits::Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    cv=read_line(type,is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS

/*! Read a circular arc point */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_point(stream & is,
           Traits::Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_xcurve(stream & is,Traits::X_monotone_curve_2 & xc)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_2(type))
  {
    xc=read_arc(type,is);
    return true;
  }
  return false;
}

/*! Read a general circular curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is, Traits::Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    Circular_kernel::Circle_2 c_temp;
    is >> c_temp;
    cv=Circular_kernel::Circular_arc_2(c_temp);
    return true;
  }
  else if (is_deg_2(type))
  {
    cv=read_arc(type,is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read a circular-line arc point */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_point(stream & is, Traits::Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_xcurve(stream & is,Traits::X_monotone_curve_2 & xc)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    xc = read_line(type,is);
    return true;
  }
  else if (is_deg_2(type))
  {
    xc = Traits::X_monotone_curve_2(read_arc(type,is));
    return true;
  }
  return false;
}

/*! Read a general circular-line curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is, Traits::Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    Circular_kernel::Circle_2 c_temp;
    is >> c_temp;
    cv=Traits::Curve_2(c_temp);
    return true;
  }
  else if (is_deg_1(type))
  {
    cv = Traits::Curve_2(read_line(type,is));
    return true;
  }
  else if (is_deg_2(type))
  {
    cv = read_arc(type,is);
    return true;
  }
  return false;

}

#endif

#endif
