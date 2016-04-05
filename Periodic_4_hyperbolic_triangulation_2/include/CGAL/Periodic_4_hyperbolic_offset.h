#ifndef CGAL_PERIODIC_4_HYPERBOLIC_OFFSET_H
#define CGAL_PERIODIC_4_HYPERBOLIC_OFFSET_H

#include <vector>
#include <string>
#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

namespace CGAL {

typedef std::vector<int> 														Index_sequence;
typedef std::pair<Hyperbolic_octagon_translation_matrix, Index_sequence>		Word;






class Periodic_4_hyperbolic_offset {

private:
	Index_sequence 												_index_sequence;
	Hyperbolic_octagon_translation_matrix 						_Hyperbolic_octagon_translation_matrix;
	Word 														_word;
	std::vector<Hyperbolic_octagon_translation_matrix>		 	_gens;

public:

	Periodic_4_hyperbolic_offset(Index_sequence iw) {
		get_generators(_gens);
		_index_sequence = iw;

		_Hyperbolic_octagon_translation_matrix = _gens[_index_sequence[0]];
		for (int i = 1; i < _index_sequence.size(); i++) {
			_Hyperbolic_octagon_translation_matrix = _Hyperbolic_octagon_translation_matrix * _gens[_index_sequence[i]];
		}

		_word = Word(_Hyperbolic_octagon_translation_matrix, _index_sequence);
	}

	void Get_word(Word& w) {
		w = Word(_word.first, _word.second);
	}


	pair<double, double> Apply_to(double x, double y) {
		return _Hyperbolic_octagon_translation_matrix.apply(x, y);
	}


};	// class Periodic_4_hyperbolic_offset



} // namespace CGAL


std::ostream& operator<<(std::ostream& os, const CGAL::Index_sequence& w)
{
   	for (int i = 0; i < w.size() - 1; i++) {
   		os << w[i] << " ";
   	}
   	os << w[w.size() - 1];
   	return os;
}

std::ostream& operator<<(std::ostream& os, const CGAL::Word& w)
{
	os << w.first << ", ";
   	for (int i = 0; i < w.second.size() - 1; i++) {
   		os << w.second[i] << " ";
   	}
   	os << w.second[w.second.size() - 1];
   	return os;
}


#endif  // CGAL_PERIODIC_4_HYPERBOLIC_OFFSET_H