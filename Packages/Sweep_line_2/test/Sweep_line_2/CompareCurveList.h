#ifndef _COMPARE_CURVE_LIST_H
#define _COMPARE_CURVE_LIST_H

#include <list>

template <class CurveList, class Traits>
class CompareCurveList
{

public:

  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_curve_2 X_curve_2;

  typedef typename CurveList::iterator CurveListIter;

  CompareCurveList() {}

  bool IsIdentical(CurveListIter begin1, CurveListIter end1,
		   CurveListIter begin2, CurveListIter end2)
  {
    Traits traits;

    CurveList list1;
    CurveListIter iter = begin1;
    while ( iter != end1 )
    {
      list1.push_back(*iter);
      //std::cout << *iter << "\n";
      ++iter;
    }
    //getc(stdin);
    int count1 = list1.size();

    int count2 = 0;
    iter = begin2;
    while ( iter != end2 ) {
      count2++;
      ++iter;
    }
    std::cout << "counts: " << count1 << " " << count2 << "\n";
    if ( count1 != count2 )
    {
      std::cout << "The list are not of the same lengths ("
		<< count1 << "," << count2 << ")\n";
      return false;
    }

    iter = begin2;
    while ( iter != end2 )
    {
      //std::cout << "looking for " << *iter << " \n" ;
      CurveListIter i;
      for ( i = list1.begin();
	    i != list1.end() ; ++i )
      {
	//std::cout << "\tcomparing to 1" << *i << " to " << *iter << "\n";
	if ( traits.curve_is_same(*i, *iter))
	{
	  //std::cout << *i << " is the same as  " << *iter << "\n";
	  break;
	}
      }
      if ( i == list1.end() ) return false;
      ++iter;
    }
    return true;
  }
};


template <class PointList, class Traits>
class ComparePointList
{

public:

  typedef typename Traits::Point_2 Point_2;
  typedef typename PointList::iterator PointListIter;

  ComparePointList() {}

  bool IsIdentical(PointListIter begin1, PointListIter end1,
		   PointListIter begin2, PointListIter end2)
  {
    Traits traits;

    PointList list1;
    PointListIter iter = begin1;
    while ( iter != end1 )
    {
      list1.push_back(*iter);
      //std::cout << *iter << "\n";
      ++iter;
    }
    //getc(stdin);
    int count1 = list1.size();

    int count2 = 0;
    iter = begin2;
    while ( iter != end2 ) {
      count2++;
      ++iter;
    }
    std::cout << "counts: " << count1 << " " << count2 << "\n";
    if ( count1 != count2 )
    {
      std::cout << "The list are not of the same lengths ("
		<< count1 << "," << count2 << ")\n";
      return false;
    }

    iter = begin2;
    while ( iter != end2 )
    {
      //std::cout << "looking for " << *iter << " \n" ;
      PointListIter i;
      for ( i = list1.begin();
	    i != list1.end() ; ++i )
      {
	//std::cout << "\tcomparing to 1" << *i << " to " << *iter << "\n";
	if ( traits.point_is_same(*i, *iter))
	{
	  //std::cout << *i << " is the same as  " << *iter << "\n";
	  break;
	}
      }
      if ( i == list1.end() ) return false;
      ++iter;
    }
    return true;
  }
};

#endif
