
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <list>
#include <fstream>
#include <iterator>
#include <stdio.h>


#include <CGAL/polyap_traits.h>
#include <CGAL/polyap_fct.h>
 

typedef double   DT; 			// base data type
typedef CGAL::Cartesian<DT> K; 	// kernel type

typedef CGAL::DistTypeTraits<K>::DataT DataType;
typedef CGAL::Point_2<K> Point2;

typedef	std::list<Point2> DataContainer;

#include"test1.h"

int main()
{
	FILE *fis;

//	fis=fopen("file.dat","r");
	fis=fopen("kk10.dat","r");

	DataContainer il;
	DataContainer ol;

	double x,y;
	DataType error;
	unsigned long m_final;

	unsigned ni=0;

	while(!feof(fis))
	{
		fscanf(fis,"%lf %lf",&x,&y);
	
		il.push_back(Point2(x,y));
		ni++;
	}
	fclose(fis);


	std::cout <<"No. points input curve= " << ni<< std::endl<< std::endl;

	CGAL::MaxSquaredEuclideanDistanceError_ls<K> MSED_ls;
	std::cout<<"Max Squared Euclidean Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,MSED_ls);

	CGAL::MaxSquaredEuclideanDistanceError_seg<K> MSED_seg;
	std::cout<<"Max Squared Euclidean Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,MSED_seg);

	CGAL::MaxLinfinityError_ls<K> MIED_ls;
	std::cout<<"Max L Infinity Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,MIED_ls);

	CGAL::MaxLinfinityError_seg<K> MIED_seg;
	std::cout<<"Max L Infinity Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,MIED_seg);

	CGAL::MaxManhattanError_ls<K> MMED_ls;
	std::cout<<"Max Manhattan Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,MMED_ls);

	CGAL::MaxManhattanError_seg<K> MMED_seg;
	std::cout<<"Max Manhattan Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,MMED_seg);

	CGAL::MaxVerticalError<K> MVED;
	std::cout<<"Max Vertical Distance Error"<<std::endl<< std::endl;
	compute(il,ol,MVED);

	CGAL::MaxBoundVol<K> MBVD;
	std::cout<<"Max Bounding Volumes Distance Error"<<std::endl<< std::endl;
//	compute(il,ol,MBVD);

	CGAL::SumSquaredEuclideanDistanceError_ls<K> SSEDE_ls;
	std::cout<<"Sum Squared Euclidean Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,SSEDE_ls);

	CGAL::SumSquaredEuclideanDistanceError_seg<K> SSEDE_seg;
	std::cout<<"Sum Squared Euclidean Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,SSEDE_seg);

	CGAL::SumLinfinityError_ls<K> SIDE_ls;
	std::cout<<"Sum L Infinity Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,SIDE_ls);

	CGAL::SumLinfinityError_seg<K> SIDE_seg;
	std::cout<<"Sum L Infinity Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,SIDE_seg);

	CGAL::SumManhattanError_ls<K> SMDE_ls;
	std::cout<<"Sum Manhattan Distance Error ls"<<std::endl<< std::endl;
	compute(il,ol,SMDE_ls);

	CGAL::SumManhattanError_seg<K> SMDE_seg;
	std::cout<<"Sum Manhattan Distance Error seg"<<std::endl<< std::endl;
	compute(il,ol,SMDE_seg);

	CGAL::SumVerticalError<K> SVDE;
	std::cout<<"Sum Vertical Distance Error"<<std::endl<< std::endl;
	compute(il,ol,SVDE);

	CGAL::SumSquaredEuclideanDistanceError_inc<K> SSEDE_inc;
	std::cout<<"Sum Squared Euclidean Distance Error Incremental Perez & Vidal Algorithm"<<std::endl;
	compute_DP(il,ol,SSEDE_inc);

	CGAL::SumSquaredVerticalDistanceError_inc<K> SSVDE_inc;
	std::cout<<"Sum Squared Vertical Distance Error Incremental Perez & Vidal Algorithm"<<std::endl;
	compute_DP(il,ol,SSVDE_inc);



/////////////////////////////////////////////////////////////////////////////////////////////
	CGAL::DouglasPeuckerHullTraits_Split<K> DPHTS;

	std::cout<<"Douglas-Peucker Path Hull Distance Split"<<std::endl;
	std::cout<<"    Rec Split Approx bounded-E (min-#) recursive"<<std::endl;

	ol.erase(ol.begin(),ol.end());
	
	error = K::make_FT(K::RT(25), K::RT(10)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::recSplitApprox(il.begin(),il.end(),&m_final,error,std::back_inserter(ol),DPHTS);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
	CGAL::DouglasPeuckerHullTraits_Build<K> DPHTB;

	std::cout<<"Douglas-Peucker Path Hull Distance Rebuild"<<std::endl;
	std::cout<<"    Rec Split Approx bounded-E (min-#) recursive"<<std::endl;

	ol.erase(ol.begin(),ol.end());
	
	error = K::make_FT(K::RT(25), K::RT(100)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::recSplitApprox(il.begin(),il.end(),&m_final,error,std::back_inserter(ol),DPHTB);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Graph Toussaint Approx bounded-E (min-#) for Euclidean Dist"<<std::endl;

	ol.erase(ol.begin(),ol.end());
	
	error = K::make_FT(K::RT(25), K::RT(100)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::GraphToussaintApprox<CGAL::DistTypeTraits<K> >(il.begin(),il.end(),&m_final,error,std::back_inserter(ol));

	std::cout <<"      Output No. points= " << m_final<< std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Graph Toussaint and Binary Search Approx bounded-# (min-E) for Euclidean Dist"<<std::endl;

	ol.erase(ol.begin(),ol.end());
	
	m_final=80;
	std::cout <<"      Input No.Points= " <<m_final<< std::endl;

	CGAL::BinarySearchApprox<CGAL::DistTypeTraits<K> >(il.begin(),il.end(),m_final,&error,std::back_inserter(ol));

	std::cout <<"      Output error= " << error<< std::endl<< std::endl;




	
	return 1;
}


