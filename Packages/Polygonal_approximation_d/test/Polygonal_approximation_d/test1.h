

template<class DataCt,class Dist>
void compute(DataCt in,DataCt out,Dist d_traits)
{
	DataType error;
	unsigned long m_final;


//////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Dyn Prog Approx min-#"<<std::endl;

	out.erase(out.begin(),out.end());

	error = K::make_FT(K::RT(10), K::RT(4)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::dynProgApprox(in.begin(),in.end(),&m_final,error,std::back_inserter(out),d_traits);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

//////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Dyn Prog Approx min-E"<<std::endl;

	out.erase(out.begin(),out.end());
	
	m_final=20;
	std::cout <<"      Input No.Points= " <<m_final<< std::endl;

	CGAL::dynProgApprox(in.begin(),in.end(),m_final,&error,std::back_inserter(out),d_traits);

	std::cout <<"      Output error= " << error<< std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Rec Split Approx bounded-E (min-#)"<<std::endl;

	out.erase(out.begin(),out.end());
	
	error = K::make_FT(K::RT(25), K::RT(100)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::recSplitApprox(in.begin(),in.end(),&m_final,error,std::back_inserter(out),d_traits);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Rec Split Approx bounded-E (min-#) iterative"<<std::endl;

	out.erase(out.begin(),out.end());
	
	error = K::make_FT(K::RT(25), K::RT(100)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::recSplitApprox_itr(in.begin(),in.end(),&m_final,error,std::back_inserter(out),d_traits);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

//////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Rec Split Approx bounded-# (min-E)"<<std::endl;

	out.erase(out.begin(),out.end());
	
	m_final=20;
	std::cout <<"      Input No.Points= " <<m_final<< std::endl;

	CGAL::recSplitApprox(in.begin(),in.end(),m_final,&error,std::back_inserter(out),d_traits);

	std::cout <<"      Output error= " << error<< std::endl<< std::endl;

//////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Graph Search Approx min-#"<<std::endl;

	out.erase(out.begin(),out.end());

	error = K::make_FT(K::RT(10), K::RT(4)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	CGAL::GraphSearchApprox(in.begin(),in.end(),&m_final,error,std::back_inserter(out),d_traits);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

//////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Graph Search Approx min-E"<<std::endl;

	out.erase(out.begin(),out.end());
	
	m_final=20;
	std::cout <<"      Input No.Points= " <<m_final<< std::endl;

	CGAL::GraphSearchApprox(in.begin(),in.end(),m_final,&error,std::back_inserter(out),d_traits);

	std::cout <<"      Output error= " << error<< std::endl;
}


template<class DataCt,class Dist>
void compute_DP(DataCt in,DataCt out,Dist d_traits)
{
	DataType error;
	unsigned long m_final;


	//////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Dyn Prog Approx min-#"<<std::endl;

	out.erase(out.begin(),out.end());

	error = K::make_FT(K::RT(10), K::RT(4)); // define a KT type using the two RT type components
	std::cout <<"      Input error= " << error<< std::endl;

	dynProgApprox(in.begin(),in.end(),&m_final,error,std::back_inserter(out),d_traits);

	std::cout <<"      Output No. points= " << m_final<< std::endl;

//////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"    Dyn Prog Approx min-E"<<std::endl;

	out.erase(out.begin(),out.end());
	
	m_final=20;
	std::cout <<"      Input No.Points= " <<m_final<< std::endl;

	dynProgApprox(in.begin(),in.end(),m_final,&error,std::back_inserter(out),d_traits);

	std::cout <<"      Output error= " << error<< std::endl;

}
