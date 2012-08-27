// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sebastien Loriot, Nicolas Carrez

template<class Kernel>
double give_weight(const CGAL::Weighted_point<typename Kernel::Point_2,typename Kernel::FT>& P){return CGAL::to_double(P.weight());}
template<class Kernel>
double give_x(const CGAL::Weighted_point<typename Kernel::Point_2,typename Kernel::FT>& P){return CGAL::to_double(P.point().x());}
template<class Kernel>
double give_y(const CGAL::Weighted_point<typename Kernel::Point_2,typename Kernel::FT>& P){return CGAL::to_double(P.point().y());}

template<class Kernel>
double give_weight(const typename Kernel::Point_2& ){return 0;}
template<class Kernel>
double give_x(const typename Kernel::Point_2& P){return CGAL::to_double(P.x());}
template<class Kernel>
double give_y(const typename Kernel::Point_2& P){return CGAL::to_double(P.y());}


template <class Kernel, class Regular,class input_DS>
void k_delaunay(Regular& rt,input_DS& input_wpt,int order){
  typedef typename Kernel::Point_2 Point_2;
  typedef typename CGAL::Weighted_point<typename Kernel::Point_2,typename Kernel::FT>  Weighted_point_2;
  
  std::vector<typename input_DS::iterator> Current_sel;//DS that will contain all possible combinaisons of k points (iterator), where k is the order
  typename input_DS::iterator it_wpt = input_wpt.begin();
  typename input_DS::iterator stop_combi = input_wpt.end();
  for(int i=0;i<order-1;++i){       //First fill the DS with the k fist wpoints
    Current_sel.push_back(it_wpt);//Useful to know when all combinaisons have been treated
    ++it_wpt;
  }
  --it_wpt;
  Current_sel.push_back(it_wpt);

  for(int i=0;i<order;++i){    //Fix end point for combinaison searching
    --stop_combi;
  }
  do{
    //To create a new combination
    typename std::vector<typename input_DS::iterator>::iterator it_it_wpt = Current_sel.end()-1;//take last element of selection
    ++(*it_it_wpt);//take next element pointed by last element of the selection
    if(*it_it_wpt==input_wpt.end()){//if we reach the end
      --(*it_it_wpt);
      while((*it_it_wpt)-1==*(it_it_wpt-1)){//looking for  a hole in the sequence
        --it_it_wpt;
      }
      --it_it_wpt;
      ++(*it_it_wpt);//take next point at previously selected locus
      ++it_it_wpt;
      for(;it_it_wpt!=Current_sel.end();++it_it_wpt){//move next pointers to points just next their predecessor
        *it_it_wpt=*(it_it_wpt-1)+1;
      }
    }
    //Create the weighted point associated to the current selection of k points
    it_it_wpt = Current_sel.begin();
    double weight = 0;
    double pt_x = 0;
    double pt_y = 0;
    for(;it_it_wpt!=Current_sel.end();++it_it_wpt){
      pt_x = pt_x + give_x<Kernel>((**it_it_wpt));
      pt_y = pt_y + give_y<Kernel>((**it_it_wpt));      
      weight = weight + order * give_weight<Kernel>((**it_it_wpt));
      //substract form the weight the sum of the squared distances between each pair of wpoints selected
      for(typename std::vector<typename input_DS::iterator>::iterator le_WptI_cgal0 = it_it_wpt+1 ;le_WptI_cgal0!=Current_sel.end();++le_WptI_cgal0){
        weight = weight - CGAL::to_double(CGAL::squared_distance(**le_WptI_cgal0,**it_it_wpt));
      }
    }
    weight = weight /(double) (order*order);
    pt_x = pt_x /(double) order;
    pt_y = pt_y /(double) order;
    rt.insert(Weighted_point_2(Point_2(pt_x,pt_y),weight));
  }while(*(Current_sel.begin())!=stop_combi);
}
