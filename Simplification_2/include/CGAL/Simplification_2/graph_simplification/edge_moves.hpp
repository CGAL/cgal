// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------

#include "utils.h"
#include "precision_safe_tests.h"

//#define EDGEMOVE_OUTPUT(line) std::cout << line << std::endl
#define EDGEMOVE_OUTPUT(line) //

namespace CGAL::internal {

	namespace detail {

		template<typename Edge_handle>
		inline Edge_handle safe_next(Edge_handle e) {
			if (e->target()->degree() == 2) {
				return e->next();
			}
			else {
				return nullptr;
			}
		}

		template<typename Edge_handle>
		inline Edge_handle safe_prev(Edge_handle e) {
			if (e->source()->degree() == 2) {
				return e->prev();
			}
			else {
				return nullptr;
			}
		}

		enum VertexType {
			// unmovable vertex (degree != 2 and != 3)
			UNMOVABLE,

			// degree 2, move shortens adjacent edge
			DEG_TWO_SUPPORT,

			// degree 2, move lengthens adjacent edge
			DEG_TWO_NO_SUPPORT,

			// degree 3, the two other adjacent edges are aligned, so vertex can just shift: they define identical tracks
			DEG_THREE_SMOOTH,

			// degree 3, needs extra vertex of deg-2 that stays at the old location, track defined by immediate adjacent edge
			DEG_THREE_SUPPORT,

			// degree 3, needs extra vertex of deg-2 that shifts, track defined by immediate adjacent edge
			DEG_THREE_NO_SUPPORT,

			// degree 3, as previous case, but the immediate adjacent edge is aligned, track defined by other adjacent edge
			DEG_THREE_ALIGNED
		};

		template<class G> struct BaseMove {

			virtual ~BaseMove() {} // necessary for inheritance to work...

			G::Edge_handle edge;

			Number<typename G::Kernel> cost;
			int queue_index = -1;

			bool blocked_by_degzero;
			std::vector<typename G::Edge_handle> blocked_by;

			virtual bool is_single() const = 0;

			bool is_blocked() const {
				return blocked_by_degzero || !blocked_by.empty();
			}
		};

		template<class G> struct SingleMove : public BaseMove<G> {

			using K = G::Kernel;
			using Vertex_handle = G::Vertex_handle;
			using Edge_handle = G::Edge_handle;

			bool left;
			enum VertexType src_type, tar_type;
			bool remove_self, remove_previous, remove_next;
			bool merge_previous, merge_next;
			Polygon<K> swept;
			bool contract_merges_highdegrees;
			int waiting_index = -1;

			bool is_single() const override {
				return true;
			}

			bool is_blocked_by(Edge_handle candidate) const {
				return is_blocked_by(candidate, swept, true);
			}

			bool is_blocked_by(Edge_handle candidate, const Polygon<K>& swept_polygon, const bool contract) const {
				if (this->edge->common_vertex(candidate) != nullptr) {
					// shares a vertex with the moving edge, part of configuration
					return false;
				}

				Edge_handle prev = this->edge->prev(left);
				Vertex_handle exclude = prev->common_vertex(candidate);

				Edge_handle next = this->edge->next(!left);
				Vertex_handle next_exclude = next->common_vertex(candidate);

				if (next_exclude != nullptr) {
					if (exclude != nullptr) {
						assert(contract);
						// convex quad, just make sure candidate is not aligned with edge itself
						Vector<K> dir_edge(this->edge->curve());
						Vector<K> dir_cand(candidate->curve());
						return safe_test::aligned_or_opposite_vectors(dir_edge, dir_cand);
					}
					else {
						exclude = next_exclude;
						// exclude comes from next, make sure it doesnt align with edge itself on contraction
						if (contract && this->remove_next) {
							// moving left: edge to next makes a left turn. test if next to candidate also makes a left turn
							// moving right: edge to next makes a right turn, test if next to candidate also makes a right turn
							if (left == CGAL::left_turn(this->edge->target()->point(), exclude->point(), candidate->other(exclude)->point())) {
								// same turns
								Vector<K> dir_edge(this->edge->curve());
								Vector<K> dir_cand(candidate->curve());
								return safe_test::aligned_or_opposite_vectors(dir_edge, dir_cand);
							}
						}
					}
				}
				else if (contract && exclude != nullptr && this->remove_previous) {
					// exclude comes from prev, make sure it doesnt align with edge itself on contraction
						// moving left: edge to prev makes a right turn. test if prev to candidate also makes a right turn
						// moving right: edge to prev makes a left turn, test if prev to candidate also makes a left turn
					if (left == CGAL::right_turn(this->edge->target()->point(), exclude->point(), candidate->other(exclude)->point())) {
						// same turns
						Vector<K> dir_edge(this->edge->curve());
						Vector<K> dir_cand(candidate->curve());
						return safe_test::aligned_or_opposite_vectors(dir_edge, dir_cand);
					}
				}

				if (exclude != candidate->source() && safe_test::convex_contains(candidate->source()->point(), swept_polygon)) {
					return true;
				}

				if (exclude != candidate->target() && safe_test::convex_contains(candidate->target()->point(), swept_polygon)) {
					return true;
				}

				const Segment<K> cc = candidate->curve();
				for (const Segment<K> e : swept_polygon.edges()) {
					auto is = CGAL::intersection(e, cc);
					if (is) {
						if (exclude != nullptr && std::holds_alternative<Point<K>>(is.value())) {
							// make sure its not the endpoint (we could have just touched the corner of swept)
							if (!safe_test::same_point(exclude->point(), std::get<Point<K>>(is.value()))) {
								return true;
							}
						}
						else {
							// anything beyond a point must be a proper intersect
							return true;
						}
					}
				}

				// no intersections were found, not blocking
				return false;
			}

			bool can_compensate_for(const SingleMove<G>& contract, const int reduc, const Number<K> area) const {
				if (this->is_blocked()) {
					return false;
				}

				if (!movable()) {
					// cannot move
					return false;
				}

				Number<K> swept = swept_area();
				if (area > swept) {
					// too small
					return false;
				}
				else if (swept == area) {
					// contract
					if (!this->contract_merges_highdegrees) {
						// not allowed to contract
						return false;
					}

					if (reduc + decrease_on_contract() <= 0) {
						// no progress
						return false;
					}
				}
				else {
					// move only
					if (increase_on_move() >= reduc) {
						// no progress
						return false;
					}
				}

				if (contract.edge->common_vertex(this->edge) != nullptr) {
					// same or neighboring edges
					return false;
				}

				Edge_handle next = safe_next(this->edge);
				if (next != nullptr && next->common_vertex(contract.edge) != nullptr) {
					// sharing next edge (== prev edge of contract)
					if (contract.remove_previous) {
						// but contract is relying on its removal
						return false;
					}
					else {
						// make sure its not convex-reflex sequence
						bool left_at_source = CGAL::left_turn(this->edge->source()->point(), next->source()->point(), next->target()->point());
						bool right_at_target = CGAL::right_turn(next->source()->point(), next->target()->point(), contract.edge->target()->point());

						if (left_at_source == right_at_target) {
							// different turns, disallow, and leave to ComboMove
							return false;
						}
					}
				}

				Edge_handle prev = safe_prev(this->edge);
				if (prev != nullptr && prev->common_vertex(contract.edge) != nullptr) {
					// sharing previous edge (== next edge of contract)
					if (contract.remove_next) {
						// but contract is relying on its removal
						return false;
					}
					else {
						// make sure its not convex-reflex sequence
						bool left_at_source = CGAL::left_turn(contract.edge->source()->point(), prev->source()->point(), prev->target()->point());
						bool right_at_target = CGAL::right_turn(prev->source()->point(), prev->target()->point(), this->edge->target()->point());

						if (left_at_source == right_at_target) {
							// different turns, disallow, and leave to CombinedMove
							return false;
						}
					}
				}

				return true;
			}

			bool movable() const {
				return src_type != UNMOVABLE && tar_type != UNMOVABLE;
			}

			bool contractable() const {
				return movable() && !contract_merges_highdegrees && decrease_on_contract() > 0;
			}

			int increase_on_move() const {
				int inc = 0;
				switch (src_type) {
				case DEG_THREE_NO_SUPPORT:
				case DEG_THREE_SUPPORT:
				case DEG_THREE_ALIGNED:
					inc++;
					break;
				}
				switch (tar_type) {
				case DEG_THREE_NO_SUPPORT:
				case DEG_THREE_SUPPORT:
				case DEG_THREE_ALIGNED:
					inc++;
					break;
				}
				return inc;
			}

			int decrease_on_contract() const {
				int dec = 0;
				if (remove_next) {
					dec++;
				}
				if (remove_self) {
					dec++;
				}
				if (remove_previous) {
					dec++;
				}
				if (merge_next) {
					dec++;
				}
				if (merge_previous) {
					dec++;
				}
				return dec - increase_on_move();
			}

			void update() {
				swept.clear();

				Vertex_handle a = nullptr, d = nullptr;
				Vertex_handle a_alt = nullptr, d_alt = nullptr;

				Vertex_handle b = this->edge->source();
				Vertex_handle c = this->edge->target();

				Line<K> baseline = Line<K>(b->point(), c->point());

				// determine start type
				switch (b->degree()) {
				case 2: {
					a = b->prev();
					if (left == CGAL::left_turn(a->point(), b->point(), c->point())) {
						src_type = detail::DEG_TWO_SUPPORT;
					}
					else {
						src_type = detail::DEG_TWO_NO_SUPPORT;
					}
					break;
				}
				case 3: {
					size_t edge_index = b->find_incident_index(this->edge);
					a = b->neighboring_incident_edge(edge_index, left)->other(b);
					a_alt = b->neighboring_incident_edge(edge_index, !left)->other(b);

					if (safe_test::point_on_line(a->point(), baseline)) {
						src_type = detail::DEG_THREE_ALIGNED;
					}
					else if (left == CGAL::left_turn(a->point(), b->point(), c->point())) {
						if (safe_test::collinear(a->point(), b->point(), a_alt->point())) {
							src_type = detail::DEG_THREE_SMOOTH;
						}
						else {
							src_type = detail::DEG_THREE_SUPPORT;
						}
					}
					else {
						src_type = detail::DEG_THREE_NO_SUPPORT;
					}
					break;
				}
				default:
					src_type = detail::UNMOVABLE;
					break;
				}

				// determine end type
				switch (c->degree()) {
				case 2: {
					d = c->next();
					if (left == CGAL::left_turn(b->point(), c->point(), d->point())) {
						tar_type = detail::DEG_TWO_SUPPORT;
					}
					else {
						tar_type = detail::DEG_TWO_NO_SUPPORT;
					}
					break;
				}
				case 3: {
					size_t edge_index = c->find_incident_index(this->edge);
					d = c->neighboring_incident_edge(edge_index, !left)->other(c);
					d_alt = c->neighboring_incident_edge(edge_index, left)->other(c);
					if (safe_test::point_on_line(d->point(), baseline)) {
						tar_type = detail::DEG_THREE_ALIGNED;
					}
					else if (left == CGAL::left_turn(b->point(), c->point(), d->point())) {
						if (safe_test::collinear(c->point(), d->point(), d_alt->point())) {
							tar_type = detail::DEG_THREE_SMOOTH;
						}
						else {
							tar_type = detail::DEG_THREE_SUPPORT;
						}
					}
					else {
						tar_type = detail::DEG_THREE_NO_SUPPORT;
					}
					break;
				}
				default:
					tar_type = detail::UNMOVABLE;
					break;
				}

				if (!movable()) {
					return;
				}

				// how far can we move along the previous track before the previous edge vanishes
				Number<K> prev_dist_squared;
				switch (src_type) {
				case DEG_TWO_SUPPORT:
				case DEG_THREE_SMOOTH:
				case DEG_THREE_SUPPORT:
					prev_dist_squared = CGAL::squared_distance(baseline, a->point());
					break;
				default:
					prev_dist_squared = -1;
					break;
				}

				// how far can we move along the next track before the next edge vanishes
				Number<K> next_dist_squared;
				switch (tar_type) {
				case DEG_TWO_SUPPORT:
				case DEG_THREE_SMOOTH:
				case DEG_THREE_SUPPORT:
					next_dist_squared = CGAL::squared_distance(baseline, d->point());
					break;
				default:
					next_dist_squared = -1;
					break;
				}

				// how far can we move until the edge itself vanishes
				Line<K> prev_track = Line<K>(b->point(), (src_type == detail::DEG_THREE_ALIGNED ? a_alt : a)->point());
				Line<K> next_track = Line<K>(c->point(), (tar_type == detail::DEG_THREE_ALIGNED ? d_alt : d)->point());

				Number<K> edge_dist_squared = -1;
				Point<K> edge_lim;

				{
					const auto is = CGAL::intersection(prev_track, next_track);
					if (is) {
						if (const Point<K>* pt = std::get_if<Point<K>>(&*is)) {
							edge_lim = *pt;
							if (left && baseline.oriented_side(edge_lim) != CGAL::ON_NEGATIVE_SIDE) {
								edge_dist_squared = CGAL::squared_distance(baseline, edge_lim);
							}
							else if (!left && baseline.oriented_side(edge_lim) != CGAL::ON_POSITIVE_SIDE) {
								edge_dist_squared = CGAL::squared_distance(baseline, edge_lim);
							}
						}
					}
				}

				// which is the limiting factor? and which edges then get removed?
				if (prev_dist_squared < 0 && next_dist_squared < 0) {
					// no support at all
					// TODO: refine, to allow converging unsupported moves?
					src_type = detail::UNMOVABLE;
					tar_type = detail::UNMOVABLE;
					return;
				}

				assert(prev_dist_squared >= 0 || next_dist_squared >= 0); // NB: if TODO above is resolved, this may fail

				Number<K> min_dist_squared;
				if (prev_dist_squared < 0) {
					min_dist_squared = next_dist_squared;
				}
				else if (next_dist_squared < 0) {
					min_dist_squared = prev_dist_squared;
				}
				else {
					min_dist_squared = CGAL::min(prev_dist_squared, next_dist_squared);
				}

				if (edge_dist_squared >= 0) {
					min_dist_squared = CGAL::min(min_dist_squared, edge_dist_squared);
				}

				remove_self = safe_test::close(edge_dist_squared, min_dist_squared);
				remove_previous = safe_test::close(prev_dist_squared, min_dist_squared);
				remove_next = safe_test::close(next_dist_squared, min_dist_squared);

				// determine if the previous edge allows merging
				if (remove_previous && b->degree() == 2 && a->degree() == 2) {
					if (remove_self) {
						// check if double-prev can merge with next
						merge_previous = safe_test::aligned_vectors(a->point() - a->prev()->point(),
							d->point() - c->point());
					}
					else {
						// check alignment with self
						merge_previous = safe_test::aligned_vectors(a->point() - a->prev()->point(),
							c->point() - b->point());
					}
				}
				else {
					merge_previous = false;
				}

				// determine if the next edge allows merging
				if (remove_next && c->degree() == 2 && d->degree() == 2) {

					if (remove_self) {
						// check if double-next can merge with prev
						merge_next = safe_test::aligned_vectors(d->next()->point() - d->point(),
							b->point() - a->point());
					}
					else {
						// check alignment with self
						merge_next = safe_test::aligned_vectors(d->next()->point() - d->point(),
							c->point() - b->point());
					}
				}
				else {
					merge_next = false;
				}

				// construct the swept area
				if (remove_self) {
					swept.push_back(edge_lim);
				}
				else if (remove_previous) {
					swept.push_back(a->point());
				}
				else {
					assert(remove_next);
					// force using the actual intersection, we (should) have eliminated exactly equal directions
					const auto is = CGAL::intersection(Line<K>(d->point(), c->point() - b->point()), prev_track);
					swept.push_back(*std::get_if<Point<K>>(&*is));
				}
				swept.push_back(b->point());
				swept.push_back(c->point());
				if (remove_self) {
					// do nothing, triangular area
				}
				else if (remove_next) {
					swept.push_back(d->point());
				}
				else {
					assert(remove_previous);
					// force using the actual intersection, we (should) have eliminated exactly equal directions
					const auto is = CGAL::intersection(Line<K>(a->point(), c->point() - b->point()), next_track);
					swept.push_back(*std::get_if<Point<K>>(&*is));
				}

				// precision, be sure we dont want to remove self...  
				if (!remove_self && safe_test::same_point(source_destination(), target_destination())) {
					remove_self = true;
					swept.erase(--swept.vertices_end());
				}

				if (remove_self && b->degree() != 2 && c->degree() != 2) {
					contract_merges_highdegrees = true;
				}
				else if (remove_previous && a->degree() != 2 && b->degree() != 2) {
					contract_merges_highdegrees = true;
				}
				else if (remove_next && c->degree() != 2 && d->degree() != 2) {
					contract_merges_highdegrees = true;
				}
				else {
					contract_merges_highdegrees = false;
				}
			}

			Number<Inexact> base_length() const {
				return std::sqrt(approximate(CGAL::squared_distance(source(), target())));
			}
			Number<Inexact> end_length() const {
				return swept.size() == 3 ? 0 : std::sqrt(approximate(CGAL::squared_distance(source_destination(), target_destination())));
			}
			Vector<Inexact> move_direction() const {
				Vector<K> vec = edge_vector();
				Vector<Inexact> dir = approximate(left ? vec.perpendicular(CGAL::COUNTERCLOCKWISE) : vec.perpendicular(CGAL::CLOCKWISE));
				return dir / std::sqrt(dir.squared_length());
			}
			Vector<K> edge_vector() const {
				return this->edge->target()->point() - this->edge->source()->point();
			}
			Point<K> source() const {
				return swept[1];
			}
			Point<K> source_destination() const {
				return swept[0];
			}
			Point<K> target() const {
				return swept[2];
			}
			Point<K> target_destination() const {
				return swept[swept.size() == 3 ? 0 : 3];
			}
			Vector<K> source_move() const {
				return source_destination() - source();
			}
			Vector<K> target_move() const {
				return target_destination() - target();
			}
			Number<K> swept_area() const {
				return CGAL::abs(swept.area());
			}

			std::optional<std::pair<Point<K>, Point<K>>> end_positions_for(const Number<K> area) const {

				Vector<Inexact> movedir = move_direction();
				Number<Inexact> len = base_length();
				Number<Inexact> end_len = end_length();
				Vector<Inexact> startTrackVec = approximate(source_move());
				Vector<Inexact> endTrackVec = approximate(target_move());
				Number<Inexact> lim_dist = CGAL::scalar_product(movedir, endTrackVec);
				Number<Inexact> r = (end_len - len) / lim_dist;

				// movearea = d * (len + d * r / 2)
				//          = d * len + d^2 * r/2
				Number<Inexact> d = utils::solveQuadraticEquationForSmallestPositive(r / 2.0, len, -approximate(area));
				if (d <= 0) {
					return std::nullopt;
				}
				Number<Inexact> ratio = CGAL::min(1.0, d / lim_dist);
				return std::pair<Point<K>, Point<K>>(source() + convert_kernel<K>(ratio * startTrackVec), target() + convert_kernel<K>(ratio * endTrackVec));
			}
		};

		template<class G> struct ComboMove : public BaseMove<G> {

			using K = G::Kernel;
			using Vertex_handle = G::Vertex_handle;
			using Edge_handle = G::Edge_handle;
			using Single = SingleMove<G>;

			bool executable;
			Single* prev_move;
			Single* next_move;
			Polygon<K> prev_swept;
			Polygon<K> next_swept;
			bool remove_self, merge_across;
			bool prev_partial, next_partial;


			bool is_single() const override {
				return false;
			}

			bool is_blocked_by(Edge_handle candidate) const {

				return prev_move->is_blocked_by(candidate, prev_swept, !prev_partial)
					|| next_move->is_blocked_by(candidate, next_swept, !next_partial);
			}

			Point<K> prev_source_destination() {
				return prev_swept[0];
			}

			Point<K> prev_target_destination() {
				return prev_swept[prev_swept.size() == 3 ? 0 : 3];
			}

			Point<K> next_source_destination() {
				return next_swept[0];
			}

			Point<K> next_target_destination() {
				return next_swept[next_swept.size() == 3 ? 0 : 3];
			}



			void update() {

				size_t edgeindex = this->edge->graph_index();

				prev_swept.clear();
				next_swept.clear();

				if (this->edge->source()->degree() != 2 || this->edge->target()->degree() != 2) {
					executable = false;
					return;
				}

				Edge_handle prev = this->edge->prev();
				Edge_handle next = this->edge->next();

				bool left_at_source = CGAL::left_turn(prev->source()->point(), this->edge->source()->point(), this->edge->target()->point());
				bool right_at_target = CGAL::right_turn(this->edge->source()->point(), this->edge->target()->point(), next->target()->point());

				if (left_at_source != right_at_target) {
					// same turn
					executable = false;
					return;
				}

				// TODO: use traits
				if (left_at_source) {
					prev_move = &(prev->data().left);
					next_move = &(next->data().right);
				}
				else {
					prev_move = &(prev->data().right);
					next_move = &(next->data().left);
				}

				if (!prev_move->movable() || !next_move->movable()) {
					executable = false;
					return;
				}

				// TODO: reduce approximate use?

				// now, determine how far both need to move until the common edge vanishes
				Vector<Inexact> edgedir = approximate(this->edge->target()->point() - this->edge->source()->point());
				Number<Inexact> edge_length = std::sqrt(edgedir.squared_length());
				edgedir /= edge_length;

				Vector<Inexact> prev_movedir = prev_move->move_direction();
				Number<Inexact> prev_len = prev_move->base_length();
				Number<Inexact> prev_len_end = prev_move->end_length();
				Vector<Inexact> prev_common_vector = approximate(prev_move->target_move());
				Number<Inexact> prev_s = 1.0 / CGAL::scalar_product(edgedir, prev_movedir); // NB: signswap between sides
				Number<Inexact> prev_lim_dist = CGAL::scalar_product(prev_movedir, prev_common_vector);
				Number<Inexact> prev_r = (prev_len_end - prev_len) / prev_lim_dist;

				Vector<Inexact> next_movedir = next_move->move_direction();
				Number<Inexact> next_len = next_move->base_length();
				Number<Inexact> next_len_end = next_move->end_length();
				Vector<Inexact> next_common_vector = approximate(next_move->source_move());
				Number<Inexact> next_s = -1.0 / CGAL::scalar_product(edgedir, next_movedir); // NB: signswap between sides
				Number<Inexact> next_lim_dist = CGAL::scalar_product(next_movedir, next_common_vector);
				Number<Inexact> next_r = (next_len_end - next_len) / next_lim_dist;

				// edge left/right of length _len expands at rate _r, covering the common edge at rate _s
				// so, moving a distance _d sweeps area
				//    _d * (_len + _d * _r / 2)
				// and covers the common edge for
				//    _d * _s;
				//
				// we need the swept areas to match, and the common edge to be fully covered
				//    p_d * p_s + n_d * n_s = e_len
				//    p_d * (p_len + p_d * p_r / 2) = n_d * (n_len + n_d * n_r / 2)
				// the former gives
				//    p_d = (e_len - n_d * n_s) / p_s
				// substituting in the latter and rewriting to quadratic form
				//     n_d^2 * (n_s ^ 2 * p_r / (2 * p_s ^ 2) - n_r / 2)
				//     + n_d (- n_s * p_len / p_s - e_len *  n_s * p_r / p_s^2 - n_len)
				//     + (e_len * p_len / p_s + e_len^2 * p_r / (2 * p_s^2)  )
				//     = 0
				// find smallest positive solution for n_d
				Number<Inexact> a = prev_r * next_s * next_s / (2 * prev_s * prev_s) - next_r / 2.0;
				Number<Inexact> b = -prev_len * next_s / prev_s - next_len - prev_r * next_s * edge_length / (prev_s * prev_s);
				Number<Inexact> c = prev_len * edge_length / prev_s + prev_r * edge_length * edge_length / (2 * prev_s * prev_s);

				Number<Inexact> next_d;
				if (safe_test::close(a, 0) && safe_test::close(b, 0) && safe_test::close(c, 0)) {
					// we'll get division by zero / NaNs... but if c is also zero, then things just line up perfectly...
					next_d = next_lim_dist;
				}
				else {
					next_d = utils::solveQuadraticEquationForSmallestPositive(a, b, c); // NB: will be <= 0 if there is no positive solution
				}
				Number<Inexact> prev_d = (edge_length - next_d * next_s) / prev_s;

				// we need to make sure that after the move, no (approximately) zero-length edges remain
				// there are a couple of cases to consider:
				// if _d > _lim_dist, then the move contracts before the edges meet
				// if _d = _lim_dist, then the move contracts as the edges meet
				// if _d < _lim_dist, then the move only has to do a partial move, for the edges to meet
				// now, if ONE of the edges contracts before meeting, then this effectively is just a pair of edge moves
				// otherwise, we need to account for the meeting point, and the possibility of either move still contracting
				prev_swept = prev_move->swept;
				next_swept = next_move->swept;

				if (next_d <= 0 || prev_d > prev_lim_dist + M_EPSILON || next_d > next_lim_dist + M_EPSILON) {
					// just a pair

					Number<K> prev_area = CGAL::abs(prev_move->swept.area());
					Number<K> next_area = CGAL::abs(next_move->swept.area());
					if (prev_area < next_area) {
						// prev contracts
						if (!prev_move->contractable()) {
							executable = false;
							return;
						}
						prev_partial = false;

						// and possibly next
						auto startend = next_move->end_positions_for(prev_area);
						if (!startend || safe_test::same_point(next_swept[0], startend->first)) {
							if (!next_move->contractable()) {
								executable = false;
								return;
							}
							next_partial = false;
						}
						else {
							next_partial = true;
							next_swept[0] = startend->first;
							if (next_swept.size() < 4) {
								next_swept.push_back(startend->second);
							}
							else {
								next_swept[3] = startend->second;
							}
						}

					}
					else {
						// next contract
						if (!next_move->contractable()) {
							executable = false;
							return;
						}
						next_partial = false;

						// and possibly prev
						auto startend = prev_move->end_positions_for(next_area);
						if (!startend || safe_test::same_point(prev_swept[0], startend->first)) {
							if (!prev_move->contractable()) {
								executable = false;
								return;
							}
							prev_partial = false;
						}
						else {
							prev_partial = true;
							prev_swept[0] = startend->first;
							if (prev_swept.size() < 4) {
								prev_swept.push_back(startend->second);
							}
							else {
								prev_swept[3] = startend->second;
							}
						}
					}

					// handling precision issues: if one is trying to remove the common edge, we should just let it be handled by the combined move
					remove_self = (!prev_partial && prev_move->remove_next) || (!next_partial && next_move->remove_previous);
				}
				else {
					// construct the meeting point
					remove_self = true;

					Number<K> prev_ratio = convert_kernel<K>(CGAL::min(1.0, prev_d / prev_lim_dist));
					Point<K> merge_point = this->edge->source()->point() + prev_ratio * convert_kernel<K>(prev_common_vector);
					Point<K> prev_other_end = prev_move->swept[1] + prev_ratio * prev_move->source_move();

					//std::cout << " meeting point " << merge_point << std::endl;
					prev_partial = !safe_test::same_point(merge_point, prev_swept[prev_swept.size() == 3 ? 0 : 3])
						&& !safe_test::same_point(merge_point, prev_other_end)
						&& !safe_test::same_point(prev_other_end, prev_move->source_destination());
					if (prev_partial) {
						prev_swept[0] = prev_other_end;
						if (prev_swept.size() < 4) {
							prev_swept.push_back(merge_point);
						}
						else {
							prev_swept[3] = merge_point;
						}
					}
					else if (!prev_move->contractable()) {
						executable = false;
						return;
					}

					Number<K> next_ratio = convert_kernel<K>(CGAL::min(1.0, next_d / next_lim_dist));
					Point<K> next_other_end = next_move->swept[2] + next_ratio * next_move->target_move();

					next_partial = !safe_test::same_point(merge_point, next_swept[0])
						&& !safe_test::same_point(merge_point, next_other_end)
						&& !safe_test::same_point(next_other_end, next_move->target_destination());
					if (next_partial) {

						next_swept[0] = merge_point;
						if (next_swept.size() < 4) {
							next_swept.push_back(next_other_end);
						}
						else {
							next_swept[3] = next_other_end;
						}
					}
					else if (!next_move->contractable()) {
						executable = false;
						return;
					}
				}

				if (remove_self) {
					Vector<K> prevdir_at_join;
					if (prev_swept.size() == 3) {
						prevdir_at_join = prev_swept[0] - prev_swept[1];
					}
					else {
						prevdir_at_join = prev_swept[2] - prev_swept[1];
					}

					Vector<K> nextdir_at_join;
					if (next_swept.size() == 3) {
						nextdir_at_join = next_swept[0] - next_swept[2];
					}
					else {
						nextdir_at_join = next_swept[2] - next_swept[1];
					}

					//std::cout << prevdir_at_join << "  ||  " << nextdir_at_join << std::endl;
					//std::cout << " -> det = " << CGAL::determinant(prevdir_at_join, nextdir_at_join) << std::endl;
					merge_across = safe_test::aligned_vectors(prevdir_at_join, nextdir_at_join);
					//std::cout << " merge_across = " << merge_across << std::endl;
				}
				else {
					merge_across = false;
				}

				int reduc = prev_partial ? -prev_move->increase_on_move() : prev_move->decrease_on_contract();
				reduc += next_partial ? -next_move->increase_on_move() : next_move->decrease_on_contract();
				if (remove_self) reduc++;
				if (merge_across) reduc++;

				executable = reduc > 0;
			}

			bool is_executable() const {
				return executable;
			}

			Number<K> swept_area() const {
				return CGAL::abs(prev_swept.area());
			}
		};


		template <class K, bool H> struct EMData {

			SingleMove<EdgeMovesGraph<K, H>> left, right;
			ComboMove<EdgeMovesGraph<K, H>> combo;
			std::vector<BaseMove<EdgeMovesGraph<K, H>>*> blocking;
		};


		template <class K, bool H> struct EMPathData {
			using Single = SingleMove<EdgeMovesGraph<K, H>>;
			std::vector<Single*> left_waiting;
			std::vector<Single*> right_waiting;

			void add_waiting(Single* single) {
				if (single->left) {
					//std::cout << ">> left move waiting for edge " << *single->edge << ", index = " << left_waiting.size() << std::endl;
					single->waiting_index = left_waiting.size();
					left_waiting.push_back(single);
				}
				else {
					//std::cout << ">> right move waiting for edge " << *single->edge << ", index = " << right_waiting.size() << std::endl;
					single->waiting_index = right_waiting.size();
					right_waiting.push_back(single);
				}
			}

			void remove_waiting(Single* single) {
				if (single->left) {
					//std::cout << ">> left move no longer waiting for edge " << *single->edge << std::endl;
					Single* other = utils::swapRemove(single->waiting_index, left_waiting);
					if (other != nullptr) {
						other->waiting_index = single->waiting_index;
					}
					single->waiting_index = -1;
				}
				else {
					//std::cout << ">> right move no longer waiting for edge " << *single->edge << std::endl;
					Single* other = utils::swapRemove(single->waiting_index, right_waiting);
					if (other != nullptr) {
						other->waiting_index = single->waiting_index;
					}
					single->waiting_index = -1;
				}
			}
		};

		template<class G>
		struct MoveQueueTraits {

			using Element_handle = BaseMove<G>*;

			static void setIndex(Element_handle m, int id) {
				m->queue_index = id;
			}

			static int getIndex(Element_handle m) {
				return m->queue_index;
			}

			static int compare(Element_handle a, Element_handle b) {
				Number<G::Kernel> ac = a->cost;
				Number<G::Kernel> bc = b->cost;
				if (ac < bc) {
					return -1;
				}
				else if (ac > bc) {
					return 1;
				}
				else {
					return 0;
				}
			}
		};
	} // namespace detail

	template <detail::EMTraits EMT>
	template <bool LEFT>
	void EdgeMoves<EMT>::update_single(Edge_handle e, Single& single) {

		single.edge = e;
		single.blocked_by_degzero = false;
		for (Edge_handle b : single.blocked_by) {
			utils::listRemove<Move>(&single, EMT::data(b).blocking);
		}
		single.blocked_by.clear();
		single.left = LEFT;

		if (single.waiting_index >= 0) {
			EMT::data(single.edge->path()).remove_waiting(&single);
			single.waiting_index = -1;
		}

		single.update();

		if (single.movable()) {

			std::vector<Single*>& waiting = LEFT ? EMT::data(single.edge->path()).right_waiting : EMT::data(single.edge->path()).left_waiting;
			for (Single* other : waiting) {
				other->waiting_index = -1;
				queue.push(other);
			}
			waiting.clear();
		}

		if (single.contractable()) {
			single.cost = EMT::determineSingleCost(single);

			if (queue.contains(&single)) {
				queue.update(&single);
			}
			else {
				queue.push(&single);
			}
		}
		else {
			queue.remove(&single);
		}

	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::update_singles(Edge_handle e) {

		assert(graph.edge(e->graph_index()) == e);

		Data& data = EMT::data(e);
		update_single<true>(e, data.left);
		update_single<false>(e, data.right);
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::update_combo(Edge_handle e) {

		Data& data = EMT::data(e);

		data.combo.edge = e;
		data.combo.blocked_by_degzero = false;
		for (Edge_handle b : data.combo.blocked_by) {
			utils::listRemove<Move>(&data.combo, EMT::data(b).blocking);
		}
		data.combo.blocked_by.clear();

		data.combo.update();

		if (data.combo.is_executable()) {
			data.combo.cost = EMT::determineComboCost(data.combo);

			if (queue.contains(&data.combo)) {
				queue.update(&data.combo);
			}
			else {
				queue.push(&data.combo);
			}
		}
		else {
			queue.remove(&data.combo);
		}
	}

	template <detail::EMTraits EMT>
	EdgeMoves<EMT>::EdgeMoves(Graph& g, EdgeTree& sqt, VertexTree& pqt)
		: graph(g), sqt(sqt), pqt(pqt), edgeset(g) {
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::initialize(bool initSQT, bool initPQT) {

		// erase deg-2 aligned vertices
		assert(graph.is_initialized());
		assert(graph.can_perform_operation());

		graph.start_operation_group();
		size_t i = 0;
		while (i < graph.number_of_vertices()) {
			Vertex_handle v = graph.vertex(i);

			if (v->degree() == 2 && safe_test::collinear(v->prev()->point(), v->point(), v->next()->point())) {
				graph.merge_vertex(v);
			}
			else {
				++i;
			}
		}
		graph.end_operation_group();

		// initialize datastructures
		if (initSQT) {
			sqt.clear();
			for (Edge_handle e : graph.edges()) {
				sqt.insert(e);
			}
		}

		if (initPQT) {
			pqt.clear();
			for (Vertex_handle v : graph.vertices()) {
				if (v->degree() == 0) {
					pqt.insert(v);
				}
			}
		}

		// initialize move data
		for (Edge_handle e : graph.edges()) {
			update_singles(e);
		}
		for (Edge_handle e : graph.edges()) {
			update_combo(e);
		}

		assert(validate_state());
	}

	template <detail::EMTraits EMT>
	bool EdgeMoves<EMT>::test_topology(Single& single) {
		// test if its blocked
		Rectangle<Kernel> rect = utils::boxOf(single.swept);

		single.blocked_by_degzero = false;
		pqt.findContained(rect, [&single](Vertex_handle b) {
			if (!single.swept.has_on_unbounded_side(b->point())) {
				single.blocked_by_degzero = true;
			}
			});

		if (single.blocked_by_degzero)
			return false;

		sqt.findOverlapped(rect, [this, &single](Edge_handle b) {
			if (single.is_blocked_by(b)) {
				EMT::data(b).blocking.push_back(&single);
				single.blocked_by.push_back(b);
			}
			});

		return single.blocked_by.empty();
	}

	template <detail::EMTraits EMT>
	bool EdgeMoves<EMT>::test_topology(Combo& combo) {
		// test if its blocked
		Rectangle<Kernel> rect = utils::boxOf(utils::boxOf(combo.prev_swept), utils::boxOf(combo.next_swept));

		combo.blocked_by_degzero = false;
		pqt.findContained(rect, [&combo](Vertex_handle b) {
			if (!combo.prev_swept.has_on_unbounded_side(b->point())
				|| !combo.next_swept.has_on_unbounded_side(b->point())) {
				combo.blocked_by_degzero = true;
			}
			});

		if (combo.blocked_by_degzero)
			return false;

		sqt.findOverlapped(rect, [this, &combo](Edge_handle b) {
			if (combo.is_blocked_by(b)) {
				EMT::data(b).blocking.push_back(&combo);
				combo.blocked_by.push_back(b);
			}
			});

		return combo.blocked_by.empty();
	}

	template <detail::EMTraits EMT>
	detail::SingleMove<typename EMT::Graph>* EdgeMoves<EMT>::find_compensate_move(Single& contract, Number<Kernel> area) {
		// find compensating move		

		using namespace detail;

		Edge_handle walkBck = safe_prev(contract.edge);
		if (walkBck != nullptr) {
			walkBck = safe_prev(walkBck);
		}
		Edge_handle walkFwd = safe_next(contract.edge);
		if (walkFwd != nullptr) {
			walkFwd = safe_next(walkFwd);
		}

		const int reduc = contract.decrease_on_contract();

		// TODO: partial blocked...?
		while (walkBck != nullptr || walkFwd != nullptr) {
			if (walkBck != nullptr) {
				// take a step back and test
				Single& candidate = contract.left ? EMT::data(walkBck).right : EMT::data(walkBck).left;
				if (candidate.can_compensate_for(contract, reduc, area)) {

					if (test_topology(candidate)) {
						return &candidate;
					}
					else {
						queue.remove(&candidate);
					}
				}

				walkBck = safe_prev(walkBck);

				if (walkBck == walkFwd) {
					// cycled around, stop
					break;
				}
			}

			if (walkFwd != nullptr) {
				// take a step forward and test
				Single& candidate = contract.left ? EMT::data(walkFwd).right : EMT::data(walkFwd).left;
				if (candidate.can_compensate_for(contract, reduc, area)) {

					if (test_topology(candidate)) {
						return &candidate;
					}
					else {
						queue.remove(&candidate);
					}
				}
				walkFwd = safe_next(walkFwd);

				if (walkBck == walkFwd) {
					// cycled around, stop
					break;
				}
			}

		}
		return nullptr;
	}

	template <detail::EMTraits EMT>
	std::optional<std::variant<detail::ComboMove<typename EMT::Graph>*, std::pair<detail::SingleMove<typename EMT::Graph>*, detail::SingleMove<typename EMT::Graph>*>, detail::SingleMove<typename EMT::Graph>*>>
		EdgeMoves<EMT>::findNextStep() {

		assert(graph.can_perform_operation());
		while (!queue.empty()) {
			Move* m = queue.peek();

			if (m->is_single()) {

				Single* sm = static_cast<Single*>(m);
				if (test_topology(*sm)) {

					Number<Kernel> area = sm->swept_area();

					if (area < M_EPSILON) { // TODO: can we avoid this check for exact mode?
						// tiny move, no compensation used
						return Operation(sm);
					}
					else {
						Single* compensate = find_compensate_move(*sm, area);
						if (compensate != nullptr) {
							return Operation(PairedSingles(sm, compensate));
						}
						else {
							EMT::data(sm->edge->path()).add_waiting(sm);
						}
					}
				}
				// if it arrived here: not a valid operation, pop it
				queue.pop();
			}
			else {

				Combo* combo = static_cast<Combo*>(m);

				if (test_topology(*combo)) {
					return Operation(combo);
				}
				else {
					// not a valid operation, pop it
					queue.pop();
				}
			}
		}

		// no move exists
		return std::nullopt;
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::start_move(Single& move) {
		using namespace detail;

		switch (move.src_type) {
		case DEG_TWO_SUPPORT:
		case DEG_TWO_NO_SUPPORT: {
			// skip
			break;
		}
		case DEG_THREE_SMOOTH: {
			// be sure to notice that the other edge will change
			//edgeset.add(move.edge->prev(!move.left));
			assert(edgeset.contains(move.edge->prev(!move.left)));
			break;
		}
		case DEG_THREE_SUPPORT: {
			Edge_handle other_prev = move.edge->prev(!move.left);
			bool add_inc = other_prev->source() == move.edge->source();
			Vertex_handle v = graph.subdivide_edge(other_prev, move.edge->source()->point(), !add_inc);
			edgeset.add(add_inc ? v->incoming() : v->outgoing());
			// NB: prev-edge on same face will be handled later
			break;
		}
		case DEG_THREE_ALIGNED:
		case DEG_THREE_NO_SUPPORT: {
			graph.subdivide_edge(move.edge, move.edge->source()->point(), false);
			// newedge is the prev edge, will be added later to edgeset
			edgeset.add(move.edge->prev());
			break;
		}
		default:
			assert(false); // unsupported vertex-moving type
		}

		switch (move.tar_type) {
		case DEG_TWO_SUPPORT:
		case DEG_TWO_NO_SUPPORT: {
			// skip
			break;
		}
		case DEG_THREE_SMOOTH: {
			// be sure to notice that the other edge will change
			//edgeset.add(move.edge->next(move.left));
			assert(edgeset.contains(move.edge->next(move.left)));
			break;
		}
		case DEG_THREE_SUPPORT: {
			Edge_handle other_next = move.edge->next(move.left);
			bool add_inc = other_next->source() == move.edge->target();
			Vertex_handle v = graph.subdivide_edge(other_next, move.edge->target()->point(), !add_inc);
			edgeset.add(add_inc ? v->incoming() : v->outgoing());
			// NB: next-edge on same face will be handled later
			break;
		}
		case DEG_THREE_ALIGNED:
		case DEG_THREE_NO_SUPPORT: {
			graph.subdivide_edge(move.edge, move.edge->target()->point(), true);
			// newedge is the next edge, will be added later to edgeset
			edgeset.add(move.edge->next());
			break;
		}
		default:
			assert(false); // unsupported vertex-moving type
		}
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::contract(Single& move) {

		EDGEMOVE_OUTPUT("Contracting " << *move.edge << (move.left ? " left" : " right") << " for area " << move.swept_area());

		start_move(move);

		Edge_handle prev = move.edge->prev(move.left);
		Edge_handle next = move.edge->next(!move.left);

		if (move.remove_self) {
			EDGEMOVE_OUTPUT("  remove self");
			Vertex_handle v;
			if (move.edge->source()->degree() != 2) {
				// target must be degree-2
				v = move.edge->source();
				graph.merge_with_next(move.edge);
			}
			else {
				// source is degree-2
				v = move.edge->target();
				graph.merge_with_prev(move.edge);
			}
			graph.move_vertex(v, move.source_destination());

			assert(!move.remove_previous || !move.remove_next);

			if (move.remove_previous) {
				if (v->degree() == 2) {
					EDGEMOVE_OUTPUT("  remove prev");
					graph.merge_with_next(prev);

					if (move.merge_previous) {
						EDGEMOVE_OUTPUT("  merge prev");
						graph.merge_with_next(next->prev());
					}
				}
				else {
					EDGEMOVE_OUTPUT("  remove prev, deg != 2");
					Vertex_handle pv = prev->other(v);
					assert(!move.merge_previous);
					assert(pv->degree() == 2);
					Edge_handle pe = graph.merge_vertex(pv, pv->outgoing() == prev); // make sure to erase prev. NB: combined with the move, the other edge at pv does not change
					graph.move_vertex(v, move.source_destination());
					//edgeset.add(pe);
				}
				//edgeset.add(next);
			}
			else if (move.remove_next) {
				if (v->degree() == 2) {
					EDGEMOVE_OUTPUT("  remove next");
					graph.merge_with_prev(next);

					if (move.merge_next) {
						EDGEMOVE_OUTPUT("  merge next");
						graph.merge_with_prev(prev->next());
					}
				}
				else {
					EDGEMOVE_OUTPUT("  remove next, deg != 2");
					Vertex_handle nv = next->other(v);
					assert(!move.merge_next);
					assert(nv->degree() == 2);
					Edge_handle ne = graph.merge_vertex(nv, nv->outgoing() == next); // make sure to erase next. NB: combined with the move, the other edge at nv does not change
					graph.move_vertex(v, move.target_destination());
					//edgeset.add(ne);
				}
				//edgeset.add(prev);
			}
			else {
				//edgeset.addAll(v->incident_edges());
			}
		}
		else {
			//edgeset.add(move.edge);

			if (move.remove_previous) {
				Vertex_handle v = move.edge->source();
				if (v->degree() == 2) {
					EDGEMOVE_OUTPUT("  remove prev");
					graph.merge_with_next(prev);
					if (move.merge_previous) {
						EDGEMOVE_OUTPUT("  merge prev");
						graph.merge_with_next(move.edge->prev());
					}
				}
				else {
					EDGEMOVE_OUTPUT("  remove prev, deg != 2");
					Vertex_handle pv = prev->other(v);
					assert(!move.merge_previous);
					assert(pv->degree() == 2);
					Edge_handle pe = graph.merge_vertex(pv, pv->outgoing() == prev); // make sure to erase prev. NB: combined with the move, the other edge at pv does not change
					graph.move_vertex(v, move.source_destination());
					//edgeset.add(pe);
				}
			}
			else {
				graph.move_vertex(move.edge->source(), move.source_destination());
				//edgeset.add(prev);
			}

			if (move.remove_next) {
				Vertex_handle v = move.edge->target();
				if (v->degree() == 2) {
					EDGEMOVE_OUTPUT("  remove next");
					graph.merge_with_prev(next);
					if (move.merge_next) {
						EDGEMOVE_OUTPUT("  merge next");
						graph.merge_with_prev(move.edge->next());
					}
				}
				else {
					EDGEMOVE_OUTPUT("  remove next, deg != 2");
					Vertex_handle nv = next->other(v);
					assert(!move.merge_next);
					assert(nv->degree() == 2);
					Edge_handle ne = graph.merge_vertex(nv, nv->outgoing() == next); // make sure to erase next. NB: combined with the move, the other edge at nv does not change
					graph.move_vertex(v, move.target_destination());
					//edgeset.add(ne);
				}
			}
			else {
				graph.move_vertex(move.edge->target(), move.target_destination());
				//edgeset.add(next);
			}
		}
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::move(Single& move, const Number<Kernel> area) {
		assert(area > 0);
		auto [s, t] = move.end_positions_for(area).value();
		this->move(move, s, t);
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::move(Single& move, const Point<Kernel> source_dest, const Point<Kernel> target_dest) {

		EDGEMOVE_OUTPUT("Moving " << *move.edge << (move.left ? " left" : " right"));

		start_move(move);

		graph.move_vertex(move.edge->source(), source_dest);
		graph.move_vertex(move.edge->target(), target_dest);
		//edgeset.add(prev);
		//edgeset.add(move.edge);
		//edgeset.add(next);
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::performStep(Single& contract, Single& compensate) {

		assert(graph.can_perform_operation());

		EDGEMOVE_OUTPUT("Paired move");

		Number<Kernel> area = contract.swept_area();
		const bool contract_compensate = safe_test::leq(compensate.swept_area(), area);

		determineCheckout(contract, true);
		determineCheckout(compensate, contract_compensate);
		checkOut();

		graph.start_operation_group();
		this->contract(contract);
		if (contract_compensate) {
			this->contract(compensate);
		}
		else {
			this->move(compensate, area);
		}
		graph.end_operation_group();

		postProcess();
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::performStep(Single& contract) {

		assert(graph.can_perform_operation());

		EDGEMOVE_OUTPUT("Tiny move");

		determineCheckout(contract, true);
		checkOut();

		graph.start_operation_group();
		this->contract(contract);
		graph.end_operation_group();

		postProcess();
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::performStep(Combo& combo) {

		assert(graph.can_perform_operation());

		EDGEMOVE_OUTPUT("Combo move " << *combo.edge << " for area " << combo.swept_area());

		determineCheckout(*combo.prev_move, !combo.prev_partial);
		determineCheckout(*combo.next_move, !combo.next_partial);
		checkOut();

		// perform

		graph.start_operation_group();
		if (combo.remove_self) {
			// degeneracy handling
			combo.prev_move->remove_next = false;
			combo.next_move->remove_previous = false;
		}
		assert(combo.prev_partial || combo.next_partial || !combo.prev_move->remove_next || !combo.next_move->remove_previous); // cannot contract both ends and both remove the common edge?

		if (combo.prev_partial) {
			move(*combo.prev_move, combo.prev_source_destination(), combo.prev_target_destination());
		}
		else {
			contract(*combo.prev_move);
		}
		if (combo.next_partial) {
			move(*combo.next_move, combo.next_source_destination(), combo.next_target_destination());
		}
		else {
			contract(*combo.next_move);
		}
		assert(graph.edge(combo.edge->graph_index()) == combo.edge); // if violated, something goes wrong with degeneracy handling above?

		if (combo.remove_self) {
			EDGEMOVE_OUTPUT("  combo: remove self");
			// NB: the two points should already be on top of eachother by the two moves
			Vertex_handle src = combo.edge->source();
			Vertex_handle tar = combo.edge->target();
			edgeset.remove(combo.edge);
			graph.merge_vertex(src, true); // removes combo.edge
			if (combo.merge_across) {
				EDGEMOVE_OUTPUT("  combo: merge across");
				edgeset.remove(tar->outgoing());
				graph.merge_vertex(tar, true); // removes tar.outgoing
			}
		}
		graph.end_operation_group();

		postProcess();
	}

	template <detail::EMTraits EMT>
	bool EdgeMoves<EMT>::validate_state() {
		bool correct_state = true;

		for (Move* m : queue.content()) {
			if (graph.edge(m->edge->graph_index()) != m->edge) {
				std::cout << "! move in queue for edge " << *m->edge << "  that is not in graph: ";
				if (m->is_single()) {
					if (static_cast<Single*>(m)->left) {
						std::cout << "left move";
					}
					else {
						std::cout << "right move";
					}
				}
				else {
					std::cout << "combo move";
				}
				std::cout << std::endl;
				correct_state = false;
			}
		}

		Graph_static_edge_map<Graph, int> emap(graph, 0);

		// Edges may lie outside the root box of the tree (boundary nodes extend to infinity), so the query must cover all edges too.
		Rectangle<Kernel> rect = sqt.root_box();
		for (Edge_handle e : graph.edges()) {
			rect = utils::boxOf(rect, EdgeQuadTreeTraits<Graph>::get_bounding_box(e));
		}
		sqt.findOverlapped(rect, [this, &emap](Edge_handle e) {
			if (this->graph.edge(e->graph_index()) != e) {
				std::cout << "! edge in SQT that is not in graph: " << *e << std::endl;
			}
			else {
				emap[e] = emap[e] + 1;
			}
			});
		for (Edge_handle e : graph.edges()) {
			if (emap[e] != 1) {
				std::cout << "! edge not exactly once in SQT: " << *e << " occurs " << emap[e] << " times" << std::endl;
			}
		}


		for (Edge_handle e : graph.edges()) {

			Data& edata = EMT::data(e);
			auto& pdata = EMT::data(e->path());

			// left check
			if (edata.left.contractable()) {
				if (edata.left.is_blocked()) {
					if (queue.contains(&edata.left)) {
						std::cout << "! move in queue that is blocked: left move for " << *e << std::endl;
						correct_state = false;
					}

					for (Edge_handle b : edata.left.blocked_by) {
						Data& bdata = EMT::data(b);
						if (std::find(bdata.blocking.begin(), bdata.blocking.end(), &edata.left) == bdata.blocking.end()) {
							std::cout << "! left move for edge " << *e << " says its blocked by edge" << *b << ", but is not in its blocking list" << std::endl;
							correct_state = false;
						}
					}
				}
				else if (edata.left.waiting_index < 0) {
					if (!queue.contains(&edata.left)) {
						std::cout << "! move not in queue that is contractable, not blocked and not waiting: left move for " << *e << std::endl;
						correct_state = false;
					}
				}
				else {
					if (pdata.left_waiting[edata.left.waiting_index] != &edata.left) {
						std::cout << "! move is waiting, but not in waiting list: left move for " << *e << std::endl;
						correct_state = false;
					}
				}
			}
			else {
				if (queue.contains(&edata.left)) {
					std::cout << "! move in queue that is not contractable: left move for " << *e << std::endl;
					correct_state = false;
				}
			}

			// right check
			if (edata.right.contractable()) {
				if (edata.right.is_blocked()) {
					if (queue.contains(&edata.right)) {
						std::cout << "! move in queue that is blocked: right move for " << *e << std::endl;
						correct_state = false;
					}

					for (Edge_handle b : edata.right.blocked_by) {
						Data& bdata = EMT::data(b);
						if (std::find(bdata.blocking.begin(), bdata.blocking.end(), &edata.right) == bdata.blocking.end()) {
							std::cout << "! right move for edge " << *e << " says its blocked by edge" << *b << ", but is not in its blocking list" << std::endl;
							correct_state = false;
						}
					}
				}
				else if (edata.right.waiting_index < 0) {
					if (!queue.contains(&edata.right)) {
						std::cout << "! move not in queue that is contractable, not blocked and not waiting: right move for " << *e << std::endl;
						correct_state = false;
					}
				}
				else {
					if (pdata.right_waiting[edata.right.waiting_index] != &edata.right) {
						std::cout << "! move is waiting, but not in waiting list: right move for " << *e << std::endl;
						correct_state = false;
					}
				}
			}
			else {
				if (queue.contains(&edata.right)) {
					std::cout << "! move in queue that is not contractable: right move for " << *e << std::endl;
					correct_state = false;
				}
			}

			// combo check
			if (edata.combo.is_executable()) {
				if (edata.combo.is_blocked()) {
					if (queue.contains(&edata.combo)) {
						std::cout << "! move in queue that is blocked: combo move for " << *e << std::endl;
						correct_state = false;
					}

					for (Edge_handle b : edata.combo.blocked_by) {
						Data& bdata = EMT::data(b);
						if (std::find(bdata.blocking.begin(), bdata.blocking.end(), &edata.combo) == bdata.blocking.end()) {
							std::cout << "! combo move for edge " << *e << " says its blocked by edge" << *b << ", but is not in its blocking list" << std::endl;
							correct_state = false;
						}
					}
				}
				else {
					if (!queue.contains(&edata.combo)) {
						std::cout << "! move not in queue that is movable and not blocked: combo move for " << *e << std::endl;
						correct_state = false;
					}
				}
			}
			else {
				if (queue.contains(&edata.combo)) {
					std::cout << "! move in queue that is not movable: combo move for " << *e << std::endl;
					correct_state = false;
				}
			}


			for (Move* m : edata.blocking) {
				if (std::find(m->blocked_by.begin(), m->blocked_by.end(), e) == m->blocked_by.end()) {
					std::cout << "! edge " << *e << " says its blocking a move for " << *m->edge << ", but is not in its blocked_by list" << std::endl;
					correct_state = false;
				}
			}

		}

		return correct_state;
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::determineCheckout(Single& move, const bool contract) {

		// TODO: refine what actually changes around high-degree vertices? (is that possible...?)
		edgeset.addAll(move.edge->source()->incident_edges());
		edgeset.addAll(move.edge->target()->incident_edges());
		if (contract) {
			if (move.merge_previous) {
				edgeset.add(move.edge->prev()->prev());
			}
			else if (move.remove_previous && move.edge->source()->degree() != 2) {
				// vertex of deg 2 will be replaced by vertex of degree != 2
				Edge_handle pe = move.edge->prev(move.left);
				edgeset.add(pe->source() == move.edge->source() ? pe->next() : pe->prev());
			}

			if (move.merge_next) {
				edgeset.add(move.edge->next()->next());
			}
			else if (move.remove_next && move.edge->target()->degree() != 2) {
				// vertex of deg 2 will be replaced by vertex of degree != 2
				Edge_handle ne = move.edge->next(!move.left);
				edgeset.add(ne->source() == move.edge->target() ? ne->next() : ne->prev());
			}
		}
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::checkOut() {

		//std::cout << "checking out..." << std::endl;
		for (Edge_handle e : edgeset) {
			//std::cout << "  " << *e << std::endl;
			sqt.remove(e);

			Data& edata = EMT::data(e);

			if (edata.left.movable()) {
				queue.remove(&edata.left);
				if (edata.left.waiting_index >= 0) {
					EMT::data(e->path()).remove_waiting(&edata.left);
				}

				for (Edge_handle b : edata.left.blocked_by) {
					utils::listRemove<Move>(&edata.left, EMT::data(b).blocking);
				}
				edata.left.blocked_by.clear();
			}
			assert(!queue.contains(&edata.left));

			if (edata.right.movable()) {
				queue.remove(&edata.right);
				if (edata.right.waiting_index >= 0) {
					EMT::data(e->path()).remove_waiting(&edata.right);
				}

				for (Edge_handle b : edata.right.blocked_by) {
					utils::listRemove<Move>(&edata.right, EMT::data(b).blocking);
				}
				edata.right.blocked_by.clear();
			}
			assert(!queue.contains(&edata.right));

			if (edata.combo.is_executable()) {
				queue.remove(&edata.combo);

				for (Edge_handle b : edata.combo.blocked_by) {
					utils::listRemove<Move>(&edata.combo, EMT::data(b).blocking);
				}
				edata.combo.blocked_by.clear();
			}
			assert(!queue.contains(&edata.combo));

			for (Move* move : edata.blocking) {
				assert(!queue.contains(move));

				utils::listRemove(e, move->blocked_by);

				if (!move->is_blocked() && !edgeset.contains(move->edge)) {

					if (move->is_single()) {
						Single* single = static_cast<Single*>(move);
						assert(single->waiting_index < 0);
						// NB: single move can also be blocked because it was checked as a compensating move, but not be actually contractable itself
						if (single->contractable())
							queue.push(move);

						// but it's movable in any case (otherwise, it would not have been blocked)
						// so it may be movable now, pop waiting from the other side
						auto& pdata = EMT::data(single->edge->path());
						std::vector<Single*>& waiting = single->left ? pdata.right_waiting : pdata.left_waiting;
						for (Single* other : waiting) {
							if (!edgeset.contains(other->edge)) {
								queue.push(other);
							}
							other->waiting_index = -1;
						}
						waiting.clear();

					}
					else {
						// Combo move: would only have been blocked if it was indeed executable
						queue.push(move);
					}
				}
			}
			edata.blocking.clear();

		}
		//edgeset.removeAll();
	}

	template <detail::EMTraits EMT>
	void EdgeMoves<EMT>::postProcess() {
		size_t changed = edgeset.size();

		for (size_t i = 0; i < changed; ++i) {

			Edge_handle e = edgeset[i];
			assert(graph.edge(e->graph_index()) == e);
			sqt.insert(e);
			edgeset.addAll(e->source()->incident_edges());
			edgeset.addAll(e->target()->incident_edges());
		}

		size_t movechanged = edgeset.size();
		for (size_t i = 0; i < movechanged; ++i) {
			Edge_handle e = edgeset[i];
			update_singles(e);
			edgeset.addAll(e->source()->incident_edges());
			edgeset.addAll(e->target()->incident_edges());
		}

		for (Edge_handle e : edgeset) {
			update_combo(e);
		}

		edgeset.removeAll();

		EDGEMOVE_OUTPUT("--- reached " << graph.number_of_edges() << " edges -------------------------");
	}

	template <detail::EMTraits EMT>
	bool EdgeMoves<EMT>::run(std::optional<std::function<bool(int, Number<Kernel>)>> stop) {

		assert(graph.can_perform_operation());
		assert(validate_state());

		while (true) {

			std::optional<Operation> next = findNextStep();
			if (!next) {
				return false;
			}

			Operation op = next.value();

			if (std::holds_alternative<Combo*>(op)) {
				Combo* combo = std::get<Combo*>(op);
				if (!stop.has_value() || (*stop)(graph.number_of_edges(), combo->cost)) {
					return true;
				}
				performStep(*combo);
			}
			else if (std::holds_alternative<PairedSingles>(op)) {
				auto [contract, compensate] = std::get<PairedSingles>(op);
				if (!stop.has_value() || (*stop)(graph.number_of_edges(), contract->cost)) {
					return true;
				}
				performStep(*contract, *compensate);
			}
			else {
				Single* contract = std::get<Single*>(op);
				if (!stop.has_value() || (*stop)(graph.number_of_edges(), contract->cost)) {
					return true;
				}
				performStep(*contract);
			}

			assert(validate_state());
		}
	}

	template <detail::EMTraits EMT>
	bool EdgeMoves<EMT>::step() {
		assert(graph.can_perform_operation());
		assert(validate_state());

		std::optional<Operation> next = findNextStep();
		if (!next) {
			return false;
		}

		Operation op = next.value();

		if (std::holds_alternative<Combo*>(op)) {
			Combo* combo = std::get<Combo*>(op);
			performStep(*combo);
		}
		else if (std::holds_alternative<PairedSingles>(op)) {
			auto [contract, compensate] = std::get<PairedSingles>(op);
			performStep(*contract, *compensate);
		}
		else {
			Single* contract = std::get<Single*>(op);
			performStep(*contract);
		}

		assert(validate_state());
		return true;
	}

	template <typename G>
	detail::EMData<typename G::Kernel, G::Graph_traits::historic> data(typename G::Edge_handle e) {
		return e->data();
	}
	template <typename G>
	Number<typename G::Kernel> BuchinEtAlTraits<G>::determineSingleCost(detail::SingleMove<G>& sm) {
		return CGAL::abs(sm.swept.area());
	}

	template <typename G>
	Number<typename G::Kernel> BuchinEtAlTraits<G>::determineComboCost(detail::ComboMove<G>& cm) {
		return CGAL::abs(cm.prev_swept.area());
	}

} // namespace CGAL::internal