

#include <map>
#include <iostream>
#include <cmath>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
//#include <CGAL/Periodic_4_hyperbolic_offset.h>
//#include <CGAL/Periodic_4_hyperbolic_offset_holder.h>
#include <CGAL/Aff_transformation_2.h>

#include <CGAL/Hyperbolic_octagon_word_4.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel                       R;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<R>          K;
typedef K::FT                                                                   FT;
typedef K::Point_2                                                              Point;
typedef K::Segment_2                                                            Segment_2;
typedef K::Arc_2                                                                Arc_2;
typedef K::Line_segment_2                                                       Line_segment_2;
typedef K::Circle_2                                                             Circle;
typedef Hyperbolic_octagon_translation_matrix<K>                                Octagon_translation_matrix;

using namespace std;

const string names[] = { "a", 
                         "\\bar{b}", 
                         "c", 
                         "\\bar{d}", 
                         "\\bar{a}", 
                         "b", 
                         "\\bar{c}", 
                         "d" };



void recurr(vector<Octagon_translation_matrix>& v, vector<string>& w, 
            vector<Octagon_translation_matrix> g, int depth = 1) {
    if (depth > 1) {
        
        recurr(v, w, g, depth-1);

        vector< Octagon_translation_matrix > tmp;
        vector<string> tmpw;
        for (int i = 0; i < v.size(); i++) {
            tmp.push_back(v[i]);
            tmpw.push_back(w[i]);
        }

        for (int i = 0; i < tmp.size(); i++) {
            for (int j = 0; j < g.size(); j++) {
                v.push_back(tmp[i]*g[j]);
                w.push_back(tmpw[i]+names[j]);
            }
        }
    } else {
        for (int i = 0; i < g.size(); i++) {
            v.push_back(g[i]);
            w.push_back(names[i]);
        }
    }
}


/*

void recurr(vector<Index_sequence>& idx, int depth = 1) {
    if (depth > 1) {
        recurr(idx, depth - 1);

        vector<Index_sequence> tmp;
        for (int i = 0; i < idx.size(); i++) {
            tmp.push_back(idx[i]);
        }
        for (int i = 0; i < tmp.size(); i++) {
            for (int j = 0; j < 8; j++) {
                Index_sequence tt = tmp[i];
                tt.push_back(j);
                idx.push_back(tt);
            }
        }
    } else {
        for (int i = 0; i < 8; i++) {
            Index_sequence tmp;
            tmp.push_back(i);
            idx.push_back(tmp);
        }
    }
}

*/

/*
void export_words_tex(vector<string> w, vector<Hyperbolic_octagon_translation_matrix> v, string filename) {

    ofstream fcout;
    fcout.open(filename);
    
    fcout << "\\documentclass{report}"  << endl;
    fcout << endl;
    fcout << "\\usepackage{fullpage}"   << endl;
    fcout << "\\usepackage{a4wide}"     << endl;
    fcout << "\\usepackage{longtable}"  << endl;
    fcout << endl;
    fcout << "\\setlength{\\parindent}{0cm}" << endl;
    fcout << endl;
    fcout << "\\begin{document}" << endl;
    fcout << endl;
    fcout << "\\begin{center}"          << endl;
    fcout << "  \\begin{longtable}{l|l|l|l}"    << endl;
    for (int i = 0; i < w.size()-1; i++) {
        fcout << "      " << i << " & $" << w[i] << "$ & " << v[i].length() << " & " << v[i] << " \\\\" << endl;
    }
    int i = w.size()-1;
    fcout << "      " << i << " & " << w[i] << " & " << v[i].length() << " & " << v[i] << endl;
    fcout << "  \\end{longtable}"       << endl;
    fcout << "\\end{center}"            << endl;
    fcout << endl;
    fcout << "\\end{document}" << endl;

    fcout.close();
}


void export_words_tex(map<Hyperbolic_octagon_translation_matrix,string> m, string filename) {

    ofstream fcout;
    fcout.open(filename);
    
    fcout << "\\documentclass{report}"  << endl;
    fcout << endl;
    fcout << "\\usepackage{fullpage}"   << endl;
    fcout << "\\usepackage{a4wide}"     << endl;
    fcout << "\\usepackage{longtable}"  << endl;
    fcout << endl;
    fcout << "\\setlength{\\parindent}{0cm}" << endl;
    fcout << endl;
    fcout << "\\begin{document}" << endl;
    fcout << endl;
    fcout << "\\begin{center}"          << endl;
    fcout << "  \\begin{longtable}{l|l|l|l}"    << endl;

    int i = 0;
    for (std::map<Hyperbolic_octagon_translation_matrix, string>::iterator it = m.begin(); it != m.end(); it++) {
        fcout << "      " << i++ << " & $" << it->second << "$ & " << it->first.length() << " & " << it->first << " \\\\" << endl;
    }

    fcout << "  \\end{longtable}"       << endl;
    fcout << "\\end{center}"            << endl;
    fcout << endl;
    fcout << "\\end{document}" << endl;

    fcout.close();
}

void export_words_tex(map<string,Hyperbolic_octagon_translation_matrix> m, string filename) {

    ofstream fcout;
    fcout.open(filename);
    
    fcout << "\\documentclass{report}"  << endl;
    fcout << endl;
    fcout << "\\usepackage{fullpage}"   << endl;
    fcout << "\\usepackage{a4wide}"     << endl;
    fcout << "\\usepackage{longtable}"  << endl;
    fcout << endl;
    fcout << "\\setlength{\\parindent}{0cm}" << endl;
    fcout << endl;
    fcout << "\\begin{document}" << endl;
    fcout << endl;
    fcout << "\\begin{center}"          << endl;
    fcout << "  \\begin{longtable}{l|l|l|l}"    << endl;

    int i = 0;
    for (std::map<string, Hyperbolic_octagon_translation_matrix>::iterator it = m.begin(); it != m.end(); it++) {
        fcout << "      " << i++ << " & $" << it->first << "$ & " << it->second.length() << " & " << it->second << " \\\\" << endl;
    }

    fcout << "  \\end{longtable}"       << endl;
    fcout << "\\end{center}"            << endl;
    fcout << endl;
    fcout << "\\end{document}" << endl;

    fcout.close();
}
*/


double edist2(double x, double y) {
    return (x*x + y*y);
}

double norm2(pair<double, double> p) {
    return sqrt(edist2(p.first, p.second));
}

double dot(pair<double, double> p1, pair<double, double> p2) {
    return (p1.first*p2.first + p1.second*p2.second);
}

double angle(pair<double, double> p1, pair<double, double> p2) {
    return acos( dot(p1, p2) / norm2(p1) / norm2(p2) );
}


double hdistorigin(pair<double, double> p) {
    double r = norm2(p);
    return 2*atanh(r);
}

double hdisthcos(pair<double, double> p1, pair<double, double> p2) {
    double hd1 = hdistorigin(p1);
    double hd2 = hdistorigin(p2);
    double theta = angle(p1, p2);

    return acosh( cosh(hd1)*cosh(hd2) - sinh(hd1)*sinh(hd2)*cos(theta) );
}


/*
Hyperbolic_octagon_translation_matrix offset_word_4(int index, vector<Hyperbolic_octagon_translation_matrix> gens) {
    Hyperbolic_octagon_translation_matrix r;
    for (int i = 0; i < 4; i++) {
        r = r * gens[index];
        index = (index + 5) % 8;
    }
    return r;
}
*/

int main(int argc, char** argv) {

    vector<Octagon_translation_matrix> g;
    get_generators(g);

    Point v0 (sqrt(sqrt(2.) + 1.) / 2. , -sqrt(sqrt(2.) - 1.) / 2.);       //  v_0
    Point s1 (sqrt(sqrt(2.) - 1.) / 2. ,  sqrt(sqrt(2.) - 1.) / 2.);       //  μ(s_1)

    Point av0 = g[DIRECTION_A_BAR].apply(v0);
    Point as1 = g[DIRECTION_A_BAR].apply(s1);
    Point mp (0.5*av0.x() + 0.5*as1.x(), 0.5*av0.y() + 0.5*as1.y());

    cout << mp << endl;


/*
    for (int i = 0; i < 8; i++) {
        Hyperbolic_octagon_translation_matrix m = offset_word_4(i, g);
        Point p = m.apply(Point(0,0));
        //cout << "With offset " << i << ": " << p << endl;
        cout << p.x() << ", " << p.y() << ";" << endl;
    }

*/

    /*
    std::vector<Point> pts;

    pts.push_back(Point(FT(  0.776887), FT( -0.321797))); // 23
    pts.push_back(Point(FT(  0.643594), FT(  0.000000))); // 22
    pts.push_back(Point(FT(  0.455090), FT(  0.455090))); // 21
    pts.push_back(Point(FT(  0.000000), FT(  0.643594))); // 13
    pts.push_back(Point(FT( -0.455090), FT(  0.455090))); //  4
    pts.push_back(Point(FT( -0.374741), FT(  0.155223))); //  6
    pts.push_back(Point(FT( -0.374741), FT( -0.155223))); //  5
    pts.push_back(Point(FT( -0.155223), FT( -0.374741))); //  9
    pts.push_back(Point(FT(  0.155223), FT( -0.374741))); // 14
    pts.push_back(Point(FT(  0.374741), FT( -0.155223))); // 18
    pts.push_back(Point(FT(  0.374741), FT(  0.155223))); // 19
    pts.push_back(Point(FT(  0.155223), FT(  0.374741))); // 15
    pts.push_back(Point(FT( -0.155223), FT(  0.374741))); // 10
    pts.push_back(Point(FT(  0.000000), FT(  0.000000))); // 12

    for (int i = 0; i < pts.size(); i++) {
        cout << pts[i] << endl;
    }
*/
/*
    pts.push_back(Point(FT( -0.776887), FT( -0.321797))); //  0
    pts.push_back(Point(FT( -0.776887), FT(  0.321797))); //  1
    pts.push_back(Point(FT( -0.643594), FT(  0.000000))); //  2
    pts.push_back(Point(FT( -0.455090), FT( -0.455090))); //  3
    pts.push_back(Point(FT( -0.321797), FT( -0.776887))); //  7
    pts.push_back(Point(FT( -0.321797), FT(  0.776887))); //  8
    pts.push_back(Point(FT(  0.000000), FT( -0.643594))); // 11
    pts.push_back(Point(FT(  0.321797), FT( -0.776887))); // 16
    pts.push_back(Point(FT(  0.321797), FT(  0.776887))); // 17
    pts.push_back(Point(FT(  0.455090), FT( -0.455090))); // 20
    pts.push_back(Point(FT(  0.776887), FT(  0.321797))); // 24
*/
    
    /*
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_A_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_B_BAR]*g[DIRECTION_A_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_C_BAR]*g[DIRECTION_B_BAR]*g[DIRECTION_A_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_D_BAR]*g[DIRECTION_C_BAR]*g[DIRECTION_B_BAR]*g[DIRECTION_A_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_B_BAR]*g[DIRECTION_C_BAR]*g[DIRECTION_D_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_C_BAR]*g[DIRECTION_D_BAR]).apply(pts[23]) << endl;
    cout << "Point[" << 23 << "] = " << pts[23] << ", image = " << (g[DIRECTION_D_BAR]).apply(pts[23]) << endl;
    */

    //pair<double, double> O(0,0);

    /*
    pair<double, double> O(-0.3218, 0.7769);

    Hyperbolic_octagon_translation_matrix t = g[DIRECTION_A];
    
    pair<double, double> image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_B];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_C];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_D];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_A_BAR];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_B_BAR];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_C_BAR];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    t = t*g[DIRECTION_D_BAR];
    image = t.apply(O.first, O.second);
    cout << "$" << t.label << "$ & " << t.length() << " & " << hdistorigin(image) << " & " << hdisthcos(O, image) << "\\\\" << endl;
    */

    /*
    typedef unsigned short int Int;
    std::vector<Int> v;
    v.push_back(3);
    v.push_back(6);
    v.push_back(1);

    Word seq(v);

    seq.append(4);

    cout << seq << endl;

    std::string str = seq.get_string();
    cout << str << endl;

    cout << "size() = " << seq.size() << endl;
    cout << "is_identity() = " << (seq.is_identity() ? "true" : "false") << endl;

    Word w;
    Dehn_reduce_word(w, seq);
    cout << "Reduced word: " << w << endl;*/

    /*
    vector<Hyperbolic_octagon_translation_matrix> gens;
    get_generators(gens);

    Point O(0,0);
    Point aO = gens[0].apply(O);
    Point ovbO = gens[1].apply(O);

    cout << "Origin:           " << O << endl;
    cout << "Image by a:       " << aO << endl;
    cout << "Image by \\bar{b}: " << ovbO << endl;

    Circle Poincare( Point(0, 0), 1*1 );


    Point  CenterA ( FT( sqrt((sqrt(2.) + 1.) / 2.) ), FT(0.) );
    FT     Radius2 ( (sqrt(2.) - 1.) / 2. );

    CGAL::Aff_transformation_2<K> rotate(CGAL::ROTATION, std::sqrt(0.5), std::sqrt(0.5));

    Circle DomainA ( CenterA,           Radius2 );
    Circle DomainBb( rotate(CenterA),   Radius2 );


    Point p(FT(atof(argv[1])), FT(atof(argv[2])));
    FT x(abs(p.x()));
    FT y(abs(p.y()));

    if (y > x) {
        FT tmp = x;
        x = y;
        y = tmp;
    }

    Point t(x, y);

    CGAL::Bounded_side bs1 = Poincare.bounded_side(t);
    CGAL::Bounded_side bs2 = DomainA.bounded_side(t);
    CGAL::Bounded_side bs3 = DomainBb.bounded_side(t);

    bool inside = true;
    inside = inside && (bs1 == CGAL::ON_BOUNDED_SIDE);
    inside = inside && (bs2 == CGAL::ON_UNBOUNDED_SIDE);
    inside = inside && (bs3 == CGAL::ON_UNBOUNDED_SIDE);

    cout << "Point " << p << " is " << (inside ? "" : "NOT ") << "inside the fundamental domain!" << endl;

    */
    //cout << "Segment is here! Getting Arc..." << endl;

    //if (Arc_2* arc = boost::get<Arc_2>(&seg)) {
    //    cout << "Yey!" << endl;
    //}

    //Arc_2 arc = boost::get<Arc_2>(seg);

    //cout << "We got the Arc!" << endl;

  
    //K::Side_of_oriented_circle_2 soc = K::Side_of_oriented_circle_2();
    //soc(seg, O);

/*
    vector<Hyperbolic_octagon_translation_matrix> g;
    get_generators(g);

    pair<double, double> p(0.23, -0.31);
    vector< pair<double, double> > img;
    vector<double>  dist;
    for (int i = 0; i < 8; i++) {
        img.push_back( g[i].apply(p.first, p.second) );
        dist.push_back( hdist(p, img[i]) );
    
        cout << "dist[" << i << "] = " << dist[i] << endl;

    }

    cout << endl;

    vector< pair<double, double> > img2;
    vector<double> dist2;
    for (int i = 0; i < 8; i++) {
        img2.push_back( g[0].apply(img[i].first, img[i].second) );
        dist2.push_back( hdist(img[i], img2[i]) );

        cout << "dist2[" << i << "] = " << dist2[i] << endl;
    }
*/

/*  
    vector<Index_sequence> cmb;
    recurr(cmb, atoi(argv[1]));

    vector<Hyperbolic_octagon_translation_matrix> trans;

    CGAL::Periodic_4_hyperbolic_offset_holder oh;

    cout << "Translating..." << endl;
    oh.Translate(trans, cmb);
    cout << "Done!" << endl;

    cout << "Working with Offset holder..." << endl;
    oh.Insert(cmb);
    cout << "Done!" << endl;

    vector<Hyperbolic_octagon_translation_matrix> tw;
    vector<Index_sequence> tiw;
    oh.Get_unique_words(tiw, tw);

    set<Hyperbolic_octagon_translation_matrix> sidx;
    for (int i = 0; i < trans.size(); i++) {
        sidx.insert(trans[i]);
    }

    cout << endl;
    cout << "Count of all words:        " << cmb.size() << endl;
    cout << "Count of Dehn subs:        " << oh.size() << endl;
    cout << "Size reduction:            " << (1.0 - (double)oh.size()/(double)cmb.size())*100.0 << "%" << endl;
    cout << endl;
    cout << "Count of all translations: " << trans.size() << endl;
    cout << "Count of set translations: " << sidx.size() << endl;
    cout << "Size reduction:            " << (1.0 - (double)sidx.size()/(double)trans.size())*100.0 << "%" << endl;


    cout << endl << "Duplicates" << endl << "-------------------------" << endl;
    for (int i = 0; i < tw.size() - 1; i++) {
        for (int j = i+1; j < tw.size(); j++) {
            if (tw[i] == tw[j]) {
                cout << "word[" << i << "] = " << tiw[i] << " and word[" << j << "] " << tiw[j] << ", matrix is " << tw[i] << endl;
            }
        }
    }
        

    */

    
    //vector<Hyperbolic_octagon_translation_matrix> g;
    //get_generators(g);
    
    /*
    
    Hyperbolic_octagon_translation_matrix id = (g[DIRECTION_A]*g[DIRECTION_B]*g[DIRECTION_C]*g[DIRECTION_D])*(g[DIRECTION_A_BAR]*g[DIRECTION_B_BAR]*g[DIRECTION_C_BAR]*g[DIRECTION_D_BAR]);
    cout << "Identity: " << id << endl;
    cout << "Label: " << id.label << endl;
    

    vector< pair<double, double> > img;
    for (int i = 0; i < 8; i++) {
        img.push_back( g[i].apply(0,0) );
        cout << "img[" << i << "] = " << img[i].first << ", " << img[i].second << endl;
    }

    */

    /*
    Hyperbolic_octagon_translation_matrix t = ( g[6]*g[1]*g[4] );
    cout << t << endl;
    Hyperbolic_octagon_translation_matrix s = ( g[7]*g[4]*g[1]*g[6]*g[3] );
    cout << s << endl;
    */

    /*
    Hyperbolic_octagon_translation_matrix t = (g[6]*g[1]*g[4]*g[7]*g[0]*g[3]);
    cout << t << endl;
    Hyperbolic_octagon_translation_matrix s = (g[0]*g[5]*g[2]*g[7]*g[4]*g[1]*g[6]*g[3]*g[6]*g[1]*g[4]*g[7]*g[1]*g[5]*g[0]*g[3]);
    cout << s << endl; 
    */

    /*
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_C);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_A_BAR); 
    idx.push_back(DIRECTION_B_BAR); 
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_D_BAR);
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_A_BAR);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_C);
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_D_BAR);
    idx.push_back(DIRECTION_B_BAR);
    */

    /* vector<Hyperbolic_octagon_translation_matrix> vi; */

    /*
    std::vector<int> idx;
    
    
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_C);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_A_BAR);
    idx.push_back(DIRECTION_B_BAR);
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_D_BAR);
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_B_BAR);
    idx.push_back(DIRECTION_A_BAR);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_B_BAR);
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_D_BAR);
    
    

    
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_C);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_A_BAR); 
    idx.push_back(DIRECTION_B_BAR); 
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_D_BAR);
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_C_BAR);
    idx.push_back(DIRECTION_A_BAR);
    idx.push_back(DIRECTION_D);
    idx.push_back(DIRECTION_C);
    idx.push_back(DIRECTION_B);
    idx.push_back(DIRECTION_A);
    idx.push_back(DIRECTION_D_BAR);
    idx.push_back(DIRECTION_B_BAR); 
    

    

    std::vector<int> dw;
    Dehn_reduce_word(dw, idx);

    CGAL::Periodic_4_hyperbolic_offset io(idx);
    CGAL::Periodic_4_hyperbolic_offset fo(dw);

    Word iw, fw;

    io.Get_word(iw);
    fo.Get_word(fw);

    cout << "Initial sequence: ";
    for (int i = 0; i < idx.size(); i++) {
        cout << idx[i] << " ";
    }
    cout << endl;

    cout << "Final sequence:   ";
    for (int i = 0; i < dw.size(); i++) {
        cout << dw[i] << " ";
    }
    cout << endl;

    cout << "Initial word: " << iw << endl;
    cout << "Final word:   " << fw << endl;

    */  

    /*
    cout << endl << "Initial word: ";
    Hyperbolic_octagon_translation_matrix t = vi[0].first;
    cout << vi[0].second << ", ";   
    for (int i = 1; i < vi.size(); i++) {
        cout << vi[i].second << ", ";
        t = t * vi[i].first;    
    }
    cout << endl << "Equivalent matrix: " << t << endl << endl;

    vector<Hyperbolic_octagon_translation_matrix> di;
    Apply_Dehn(di, vi);
    cout << "Dehn's word:  ";
    Hyperbolic_octagon_translation_matrix s = di[0].first;
    cout << di[0].second << ", ";   
    for (int i = 1; i < di.size(); i++) {
        cout << di[i].second << ", ";
        s = s * di[i].first;    
    }
    cout << endl << "Equivalent matrix: " << s << endl << endl;

    cout << "Equivalence achieved: " << (t == s ? "YES" : "NO") << endl;
    */

    return 0;
}


