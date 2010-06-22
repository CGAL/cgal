#include <CGAL/Random.h>
#include <CGAL/IO/Dxf_variant_reader.h>

template <class CK,class ArcContainer>
void _bench_input_dfx(ArcContainer& arcs,char* Dfxfile){
	
	typedef typename CK::Circular_arc_2                                  Circular_arc_2;
 	typedef typename CK::Line_arc_2                                      Line_arc_2;
	std::ifstream fin;
	
	fin.open (Dfxfile);
	CGAL::variant_load<CK, Circular_arc_2, Line_arc_2>(fin, std::back_inserter(arcs));
	fin.close();
}

template <class CK,class ArcContainer>
void _bench_input_grid(ArcContainer& ac){
 	
	typedef typename CK::Point_2                    Point_2;
 	typedef typename CK::Circle_2                    Circle_2;
  	
	int x = 5;
  	int y = 5;
  	int r = 5;
  	
	for(int i = 0; i < 20; i++){
    		for(int j = 0; j <20; j++){
      			ac.push_back(typename CK::Circular_arc_2 (  Circle_2( Point_2(x + j*r, y + i*r), r*r)));
    		}
  	}
}

template <class CK,class ArcContainer>
void _bench_input_grid2(ArcContainer& ac){
 	
	typedef typename CK::Point_2                    Point_2;
 	typedef typename CK::Circle_2                    Circle_2;

 	int x = 10;
 	int y = 10;
 	int r = 10;
 
  	for(int i = 0; i < 20; i++){
    		for(int j = 0; j < 20; j++){
     			ac.push_back(typename CK::Circular_arc_2 ( Circle_2( Point_2(x + j*r, y + i*r), r*r) ));
    		}
  	}

  	x += r;
  	y += r;
  	
	for(int i = 0; i < 20; i++){
    		for(int j = 0; j < 20; j++){
      			ac.push_back(typename CK::Circular_arc_2 (  Circle_2( Point_2(x + j*r, y + i*r), r*r)));
    		}
  	}

}

template <class CK,class ArcContainter>
void _bench_input_random(ArcContainter& arcs)
{
	std::ofstream fout("random.inp");
 	typedef typename CK::Point_2                    Point_2;
 	typedef typename CK::Circle_2                    Circle_2;
 
 	CGAL::Random generatorOfgenerator;
  	int random_seed = generatorOfgenerator.get_int(0, 123456);
  	std::cout << "random_seed = " << random_seed << std::endl; 
 	CGAL::Random theRandom(random_seed);

	double random_max = 10;
	double random_min = -10;

    	for(int i = 0; i < 100 ; i++){
     		double x = theRandom.get_double(random_min,random_max);
     		double y = theRandom.get_double(random_min,random_max);
     		double r = theRandom.get_double(0.1,random_max);
     		fout << x << " " << y << " " << r << " " ;
    		arcs.push_back(typename CK::Circular_arc_2 (  Circle_2( Point_2(x,y), r*r) ));
    	}
	
	fout.close();
	std::cout << arcs.size() << std::endl;

}

template <class CK,class ArcContainter>
void _bench_input_file(ArcContainter& arcs){
	double x,y, r;
	typedef typename CK::Point_2                    Point_2;
	typedef typename CK::Circle_2                    Circle_2;
      
      	std::ifstream fin("random.inp");
      	for(int i = 0; i < 100 ; i++){
		fin >> x >> y >> r ;
         	arcs.push_back(typename CK::Circular_arc_2 (  Circle_2( Point_2(x,y), r*r)));
      	}
      	
	fin.close();
}

template <class CK,class ArcContainter>
void _bench_input_rotation(ArcContainter& ac){ 

	typedef typename CK::Point_2                    Point_2;
  	typedef typename CK::Circle_2                    Circle_2;
 	
	int x = 10;
        int y = 10;
 	int r = 10;
	int R =0;
	
	for(int k = 0; k < 5; k++){
 		for(int i = 0; i < 5; i++){
    			for(int j = 0; j < 5; j++){
				R = (r+k)*(r+k);
     				ac.push_back(typename CK::Circular_arc_2 (  Circle_2( Point_2(x + i*r, y + j*r),R) ));
     			}
 		}
	}

}

