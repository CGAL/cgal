// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_Convex_Polygon_2.C
// package       : bops (2.2)
// source        : include/CGAL/bops_Convex_Polygon_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_CONVEX_POLYGON_2_C
#define CGAL_BOPS_CONVEX_POLYGON_2_C

#include <CGAL/bops_Convex_Polygon_2.h>

CGAL_BEGIN_NAMESPACE

#define ORIGINAL // now use this
#ifdef ORIGINAL
template <class R_type, class Container>
    Polygon_2<Polygon_traits_2<R_type>, Container>
       Convex_Intersection(
           Polygon_2<Polygon_traits_2<R_type>, Container> P,
           Polygon_2<Polygon_traits_2<R_type>, Container> Q)


{
   typedef Polygon_2<Polygon_traits_2<R_type>,Container> My_Polygon;
   typedef Point_2<R_type> My_Point;
   typedef Segment_2<R_type> My_Segment;
   typedef Vector_2<R_type> My_Vector;

   enum tInFlag {unknown, Pin, Qin};  

   int pSize, qSize;                  // size of polygons
   int aAdv, bAdv;                    // number of iterations over each polygon
   My_Point ipoint,targ;              // variables for intersection
   My_Polygon::Vertex_circulator pCir, qCir, pCir1, qCir1;
                                      // circulators over polygons 
   Orientation bHA, aHB; // orientation of polygons
   typedef typename R_type::FT FT;
   FT cross;         // cross product
   tInFlag inflag;                    // which polygon is inside
   My_Vector aVec, bVec;              // vectors for computation
   My_Segment s1,s2,iseg;             // segments to check intersection
   Object result;                // result of intersection of s1 and s2
   std::list<My_Point> pt_list;            // points forming solution-polynom
   int lastaAdv, lastbAdv;            // variables to avoid endless loop

   //check correct orientation for algorithm
   if (!(P.orientation()==COUNTERCLOCKWISE)) P.reverse_orientation();
   if (!(Q.orientation()==COUNTERCLOCKWISE)) Q.reverse_orientation(); 

   // initialize variables
   lastaAdv = lastbAdv = 0;
   aAdv = bAdv = 0;
   pSize = P.size(); qSize = Q.size();
   inflag = unknown;
   pCir = pCir1 = P.vertices_circulator ();
   qCir = qCir1 = Q.vertices_circulator ();
   pCir1 = pCir; qCir1 = qCir;
   pCir1--; qCir1--;   

   do
     {
       // create vectors
       aVec = *pCir-*pCir1; bVec=*qCir-*qCir1;
       // calculate 'cross-product'
       cross = aVec[0]*bVec[1]-aVec[1]*bVec[0];
       // get orientation
       bHA = orientation(*pCir1, *pCir, *qCir);
       aHB = orientation(*qCir1, *qCir, *pCir);
       // create segments for intersection
       Segment_2<R_type> s1(*pCir1,*pCir);
       Segment_2<R_type> s2(*qCir1,*qCir);
       //check for intersection
       if (do_intersect(s1,s2))
	 {
	   // get type of intersection
	   result = intersection(s1,s2);
	   // if result is a point, insert in result-polynom
	   if (assign(ipoint, result))
	     {
	       // if first intersection start algorithm
	       if (inflag == unknown) aAdv=bAdv=0;
	       if (aHB==LEFTTURN) {inflag=Pin;}
	       else if (aHB==RIGHTTURN) {inflag=Qin;}
	       else
		 {
       		   // aHB is collinear
		   if (*pCir == *qCir)
		     {
		         //the endpoints of the segments are the same
		         //check which polynom is inside
			 if (orientation(*pCir,*(qCir++),*(pCir++))==RIGHTTURN) {inflag = Qin;}
			 else
			   {inflag=Pin;}
                         pCir--;qCir--;
		       } // end of pCir==qCir
		   else if (ipoint == *pCir)
		       {
			 //intersection-point is end of P, i.e. it lies on Q 
			 if (orientation(ipoint,*(pCir++),*qCir)==RIGHTTURN)
			   {inflag=Pin;}
			 else  
                           {inflag = Qin;}
			 pCir--;
		       } // end of ipoint==pCir
		   else    
		       {
                          
			 //point of intersection is end of Q, i.e. it lies on P
                         if (orientation(ipoint,*pCir,*(qCir++))==LEFTTURN)
			     {inflag=Qin;}
			 else
			     {inflag=Pin;}  
			 qCir--;
		       } // end of default      
		 } //end of else
	       if (!pt_list.empty())
		   {if ((pt_list.back()!=ipoint)&&(pt_list.front()!=ipoint)) pt_list.push_back(ipoint);}
	       else
		 pt_list.push_back(ipoint);

	     } // end of point intersection
	   else
	     {
	       // segment intersection              
               if (assign(iseg,result))
		 {
		   // make sure that pCir or qCir is target
                   if ((iseg.source()==*pCir) || (iseg.source()==*qCir))
		     {
		      targ=iseg.source();
                     }
		   else
		     {targ=iseg.target();}
		   // if P and Q have the same endpoint
		   if (*pCir == *qCir)
		     {
		       if (orientation(*pCir,*(qCir++),*(pCir++))==RIGHTTURN)
		         { inflag = Qin;}
		       else
		         { inflag = Pin;}
		       pCir--;qCir--;
		     } // end of if (*pCir == *qCir)
		   else if (targ==*pCir)
		     {
		       // segment of intersection is on end of P
                       if ((pt_list.back()!=targ)&&(pt_list.front()!=targ)) pt_list.push_back(targ);
		       inflag = Pin;
		     } //end of else if
		   else
		     {
     		       if ((pt_list.back()!=targ)&&(pt_list.front()!=targ)) pt_list.push_back(targ);
		       inflag = Qin;
		     } // end of else
               } //end of segment intersection
	     } //end of else for intersection
	 } //end of intersection
       
       if (cross >= FT(0) )
	 {
	   if (bHA == LEFTTURN)
 	     {
	       if ((inflag == Pin)&&(pt_list.back()!=*pCir)&&(pt_list.front()!=*pCir)) pt_list.push_back(*pCir);
	       aAdv++; pCir++; pCir1++;
	       // abortion check
	       lastaAdv++; lastbAdv=0;
	     }
	   else
	     { 
	       if ((inflag == Qin)&&(pt_list.back()!=*qCir)&&(pt_list.front()!=*qCir)) pt_list.push_back(*qCir);
	       bAdv++; qCir++; qCir1++;
	       // abortion check
	       lastbAdv++; lastaAdv=0;
	     }
	 }
       else
	 {
	   if (aHB == LEFTTURN)
	     {
	       if ((inflag == Qin)&&(pt_list.back()!=*qCir)&&(pt_list.front()!=*qCir)) pt_list.push_back(*qCir);
	       bAdv++; qCir++; qCir1++;
	       // abortion check
	       lastbAdv++; lastaAdv=0;

	     }
	   else
	     {  
	       if ((inflag == Pin)&&(pt_list.back()!=*pCir)&&(pt_list.front()!=*pCir)) pt_list.push_back(*pCir);
	       aAdv++; pCir++; pCir1++;
	       // abortion check
	       lastaAdv++; lastbAdv=0;
	     }
	 }
     } //end of do-loop
   while (((aAdv < pSize) || (bAdv < qSize))&&(lastaAdv <= pSize)&&(lastbAdv<=qSize));
   //create result polygon
   My_Polygon final(pt_list.begin(),pt_list.end());

   if (inflag == unknown)
     {
       //polygons do not intersect
       if (P.bounded_side(*qCir)==ON_BOUNDED_SIDE) return Q;
       else if (Q.bounded_side(*pCir)==ON_BOUNDED_SIDE) return P;
       else return final;
     }
   else
     {
       return final;
     }
}

#else

template <class R_type, class Container>
Polygon_2<Polygon_traits_2<R_type>, Container>
Convex_Intersection(
           Polygon_2<Polygon_traits_2<R_type>, Container> P,
           Polygon_2<Polygon_traits_2<R_type>, Container> Q
)
{
   typedef Polygon_2<Polygon_traits_2<R_type>,Container> My_Polygon;
   typedef Point_2<R_type> My_Point;
   typedef Segment_2<R_type> My_Segment;
   typedef Vector_2<R_type> My_Vector;
   typedef typename R_type::FT FT;

   enum tInFlag {unknown, Pin, Qin};  

   int pSize, qSize;                  // size of polygons
   int aAdv, bAdv;                    // number of iterations over each polygon
   My_Point ipoint,targ;              // variables for intersection
   My_Polygon::Vertex_circulator pCir, qCir, pCir1, qCir1;
                                      // circulators over polygons 
   Orientation bHA, aHB; // orientation of polygons
   FT cross;         // cross product
   tInFlag inflag;                    // which polygon is inside
   My_Vector aVec, bVec;              // vectors for computation
   My_Segment s1,s2,iseg;             // segments to check intersection
   Object result;                // result of intersection of s1 and s2
   std::list<My_Point> pt_list;            // points forming solution-polynom
   int lastaAdv, lastbAdv;            // variables to avoid endless loop

   //check correct orientation for algorithm
   if (!(P.orientation()==COUNTERCLOCKWISE)) P.reverse_orientation();
   if (!(Q.orientation()==COUNTERCLOCKWISE)) Q.reverse_orientation(); 

   // initialize variables
   lastaAdv = lastbAdv = 0;
   aAdv = bAdv = 0;
   pSize = P.size(); qSize = Q.size();
   inflag = unknown;
   pCir = pCir1 = P.vertices_circulator ();
   qCir = qCir1 = Q.vertices_circulator ();
   pCir1 = pCir; qCir1 = qCir;
   pCir1--; qCir1--;   

   do {
       // create vectors
       aVec = *pCir-*pCir1; bVec=*qCir-*qCir1;
       // calculate 'cross-product'
       cross = aVec[0]*bVec[1]-aVec[1]*bVec[0];
       // get orientation
       bHA = orientation(*pCir1, *pCir, *qCir);
       aHB = orientation(*qCir1, *qCir, *pCir);
       // create segments for intersection
       Segment_2<R_type> s1(*pCir1,*pCir);
       Segment_2<R_type> s2(*qCir1,*qCir);
       //check for intersection
       if (do_intersect(s1,s2)) {
	   // get type of intersection
	   result = intersection(s1,s2);
	   // if result is a point, insert in result-polynom
	   if (assign(ipoint, result)) {
	       // if first intersection start algorithm
	       //if (inflag == unknown) aAdv=bAdv=0;
	       if (aHB==LEFTTURN) {inflag=Pin;}
	       else if (aHB==RIGHTTURN) {inflag=Qin;}
	       else {
       		   // aHB is collinear
		      if (*pCir == *qCir) {
		         //the endpoints of the segments are the same
		         //check which polynom is inside
				if (orientation(*pCir,*(qCir++),*(pCir++))==RIGHTTURN) {
				  inflag = Qin;
				}
				else {
				  inflag=Pin;
				}
				pCir--;qCir--;    
			   } // end of pCir==qCir
			   else if (ipoint == *pCir) {
				 //intersection-point is end of P, i.e. it lies on Q 
				 if (orientation(ipoint,*(pCir++),*qCir)==RIGHTTURN) {
					inflag=Pin;
				 }
				 else {
					inflag = Qin;
				 }
				 pCir--;
			   } // end of ipoint==pCir
			   else {          
				 //point of intersection is end of Q, i.e. it lies on P
        		 if (orientation(ipoint,*pCir,*(qCir++))==LEFTTURN) {
					inflag=Qin;
				 }
			     else {
					inflag=Pin;
				}
					qCir--;
			   } // end of default      
		  } //end of else
		  
	      if (!pt_list.empty() ) {
		  	if ((pt_list.back()!=ipoint)&&(pt_list.front()!=ipoint))
				pt_list.push_back(ipoint);
		  }
	      else
		 	pt_list.push_back(ipoint);

	   } // end of point intersection
	   else {
	       // segment intersection              
         if (assign(iseg,result)) {
		   // make sure that pCir or qCir is target
           targ= ((iseg.source()==*pCir) || (iseg.source()==*qCir)) ?
		      		iseg.source() : iseg.target();
			  
		   // if P and Q have the same endpoint
		   if (*pCir == *qCir) {
		   		if (orientation(*pCir,*(qCir++),*(pCir++))==RIGHTTURN) {
		   			inflag = Qin;
		   		}
		    	else {
				   inflag = Pin;
				}
		    	pCir--;qCir--;
		   } // end of if (*pCir == *qCir)
		   else if (targ==*pCir) {
		        // segment of intersection is on end of P
           		if ((pt_list.back()!=targ)&&(pt_list.front()!=targ))
					pt_list.push_back(targ);
		   		inflag = Pin;
		   } //end of else if
		   else {
     	   		if ((pt_list.back()!=targ)&&(pt_list.front()!=targ))
					pt_list.push_back(targ);
		        inflag = Qin;
		   } // end of else
         } //end of segment intersection
	   } //end of else for intersection
	 } //end of intersection
       
     if (cross >= FT(0)) {
	   if (bHA == LEFTTURN) {
	       if ((inflag == Pin)&&(pt_list.back()!=*pCir)&&(pt_list.front()!=*pCir))
		   		pt_list.push_back(*pCir);
	       aAdv++; pCir++; pCir1++;
	       // abortion check
	       lastaAdv++; lastbAdv=0;
	   }
	   else { 
	       if ((inflag == Qin)&&(pt_list.back()!=*qCir)&&(pt_list.front()!=*qCir))
		   		pt_list.push_back(*qCir);
	       bAdv++; qCir++; qCir1++;
	       // abortion check
	       lastbAdv++; lastaAdv=0;
	   }
	 }
     else {
	   if (aHB == LEFTTURN) {
	       if ((inflag == Qin)&&(pt_list.back()!=*qCir)&&(pt_list.front()!=*qCir))
		   		pt_list.push_back(*qCir);
	       bAdv++; qCir++; qCir1++;
	       // abortion check
	       lastbAdv++; lastaAdv=0;
	   }
	   else {  
	       if ((inflag == Pin)&&(pt_list.back()!=*pCir)&&(pt_list.front()!=*pCir))
		   		pt_list.push_back(*pCir);
	       aAdv++; pCir++; pCir1++;
	       // abortion check
	       lastaAdv++; lastbAdv=0;
	   }
	 }
   } //end of do-loop
   while (((aAdv < pSize) || (bAdv < qSize))&&(lastaAdv <= pSize)&&(lastbAdv<=qSize));
   //create result polygon
   My_Polygon final(pt_list.begin(),pt_list.end());

   if (inflag == unknown) {
       //polygons do not intersect
       if (P.bounded_side(*qCir)==ON_BOUNDED_SIDE)
	   		return Q;
       else if (Q.bounded_side(*pCir)==ON_BOUNDED_SIDE)
	   		return P;
       else
	   		return final;
   }
   else {
       return final;
   }
}
#endif

template< class I>
Bops_Convex_Polygons_2<I>::InFlag
Bops_Convex_Polygons_2<I>::handleIntersectionPoint(
	const Point& pt, const Orientation& aHB
)
{
	InFlag inflag;
	if ( aHB == LEFTTURN ) 
		inflag= Pin;
	else if ( aHB == RIGHTTURN )
		inflag= Qin;
	else { // aHB is collinear
		if (*_pCir == *_qCir) {
			//the endpoints of the segments are the same
			//check which polynom is inside
			inflag= is_leftturn(*_pCir, *(_qCir++), *(_pCir++)) ? Pin : Qin;
			_pCir--;
			_qCir--;
		} // end of pCir==qCir
		else if (pt == *_pCir) {
			//intersection-point is end of P, i.e. it lies on Q 
			inflag= is_leftturn(pt, *(_pCir++), *_qCir) ? Qin : Pin;
			_pCir--;
		} // end of ipoint==pCir
		else {          
			//point of intersection is end of Q, i.e. it lies on P
			inflag= is_leftturn(pt, *_pCir, *(_qCir++)) ? Qin : Pin;
			_qCir--;
		} // end of default      
	} //end of else
	insert(pt);
	return inflag;
}



template< class I>
Bops_Convex_Polygons_2<I>::InFlag
Bops_Convex_Polygons_2<I>::handleIntersectionSegment(const Segment& seg)
{
	InFlag inflag;
	// make sure that pCir or qCir is target
	_Point_2 targ= (seg.source()==*_pCir || seg.source()==*_qCir) ?
			seg.source() : seg.target();
	// if P and Q have the same endpoint
	if (*_pCir == *_qCir) {
		inflag= is_leftturn(*_pCir,*(_qCir++),*(_pCir++)) ? Pin : Qin; // is_rightturn
		   _pCir--; _qCir--;
	} // end of if (*pCir == *qCir)
	else { // segment of intersection is on end of P or Q
		inflag= (targ==*_pCir) ? Pin : Qin;
		insert(targ);
	} // end of else
	return inflag;
}


template< class I>
void Bops_Convex_Polygons_2<I>::mainProcedure()
{
	Orientation bHA, aHB; 		// orientation of polygons
	NT cross;         			// cross product

	Segment s1,s2,iseg;       // segments to check intersection
	Point ipoint, targ;       // variables for intersection
	Vector aVec, bVec;        // vectors for computation
	CGAL::Object result;              // result of intersection of s1 and s2

	do {
		// create vectors
		aVec= *_pCir-*_pCir1;
		bVec= *_qCir-*_qCir1;
		// calculate 'cross-product'
		cross = aVec[0]*bVec[1]-aVec[1]*bVec[0];
		// get orientation
		bHA = orientation(*_pCir1, *_pCir, *_qCir);
		aHB = orientation(*_qCir1, *_qCir, *_pCir);
		// create segments for intersection
		Segment s1(*_pCir1,*_pCir);
		Segment s2(*_qCir1,*_qCir);
		//check for intersection
		if (do_intersect(s1,s2)) {
			// get type of intersection
			result = intersection(s1,s2);
			// if result is a point, insert in result-polynom
			if ( assign(ipoint, result) )
				_inflag= handleIntersectionPoint(ipoint, aHB);
			else if ( assign(iseg,result) ) // segment intersection   
				_inflag= handleIntersectionSegment( iseg );
			else // empty intersection
			;
		} //end of intersection

		bool is_in_A;
		if (cross >= NT(0)) 
			is_in_A=  (bHA == LEFTTURN) ? true : false;
		else 
			is_in_A= !(aHB == LEFTTURN) ? true : false;
		
		if( is_in_A == true ) {
			conditional_insert(Pin, _inflag, *_pCir);
			advancePolygonA();
		}
		else { // otherwise is_in_B == true
			conditional_insert(Qin, _inflag, *_qCir);
			advancePolygonB();
		}
		
		
   	} //end of do-loop
   	while ( !wentAround() );
	return;
}

CGAL_END_NAMESPACE

#endif // CGAL_BOPS_CONVEX_POLYGON_2_C
