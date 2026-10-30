// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------

#include <CGAL/Simplification_2/data_structures/graph_map_2.h>
#include <CGAL/Simplification_2/data_structures/quad_tree.h>

#include "utils.h"
#include "vertex_quad_tree.h"

namespace CGAL::internal {

	namespace detail {

		template<class Graph>
		struct AngleRestriction {

			using Vertex_handle = Graph::Vertex_handle;
			using Edge_handle = Graph::Edge_handle;
			using Kernel = Graph::Kernel;
			using Vec = Vector<Inexact>;
			using Num = Number<Inexact>;

			inline Num approxSquaredLength(Edge_handle e) {
				return approximate((e->source()->point() - e->target()->point()).squared_length());
			}
			inline Vec vector(Vertex_handle v, Edge_handle e) {
				auto start = approximate(v->point());
				auto end = approximate(e->other(v)->point());
				return end - start;
			}

			inline Vec direction(Vertex_handle v, Edge_handle e) {
				auto start = approximate(v->point());
				auto end = approximate(e->other(v)->point());
				return (end - start) / std::sqrt(CGAL::squared_distance(start, end));
			}

			enum class EdgeType {
				ALIGN,
				UNALIGN,
				EVADING,
				DEV_ALIGN,
				DEV_UNALIGN,
				UNDETERMINED
			};

			struct EdgeData {
				Edge_handle edge = nullptr;
				Vertex_handle significant = nullptr;
				std::pair<int, int> associated = std::pair<int, int>(-1, -1);
				int assigned = -1;
				EdgeType type = EdgeType::UNDETERMINED;
				int stepcount = -1;
				Num stepdist = -1;
				Polygon<Inexact> interference;

				inline Rectangle<Inexact> box() {
					return utils::boxOf(interference);
				}

				inline int other_direction() {
					return assigned == associated.first ? associated.second : associated.first;
				}
			};

			struct InterferenceQuadTreeTraits {

				using Element = EdgeData*;
				using Kernel = Inexact;

				static Rectangle<Inexact> get_bounding_box(EdgeData* elt) {
					return elt->box();
				};

				static bool element_overlaps_rectangle(EdgeData* elt, const Rectangle<Inexact>& rect) {
					// TODO: precise overlap check?
					return !utils::disjoint(elt->box(), rect);
				}
			};

			inline Segment<Inexact> directed_segment(Edge_handle e, EdgeData& d) {
				return Segment<Inexact>(approximate(d.significant->point()), approximate(e->other(d.significant)->point()));
			}

			using InterferenceTree = CGAL::internal::QuadTree<InterferenceQuadTreeTraits>;

			Graph& graph;
			std::vector<Vec> directions;
			Graph_static_vertex_map<Graph, bool> significant_vertices;
			Graph_static_edge_map<Graph, EdgeData> edge_data;
			InterferenceTree interference_tree;

			std::optional<std::function<void(std::string, int, int)>> progress;

			AngleRestriction(Graph& graph, std::vector<Vec> dirs, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt) :
				graph(graph), directions(dirs), significant_vertices(graph, false), edge_data(graph), progress(progress), interference_tree(approximate(graph.bounding_rectangle()), 10, 0.05) {
			};

			bool is_approximately(Vec a, Vec b) {
				const Num eps = 0.000001;
				return std::abs(a.x() - b.x()) < eps
					&& std::abs(a.y() - b.y()) < eps;
			}

			std::pair<int, int> associated_directions(Vertex_handle v, Edge_handle e) {

				Vec dir = direction(v, e);

				int i = 0;
				if (is_approximately(directions[i], dir)) {
					return std::pair<int, int>(i, -1);
				}
				int in = 1;
				int ndirs = directions.size();
				while (i < ndirs) {
					if (is_approximately(directions[in], dir)) {
						return std::pair<int, int>(in, -1);
					}
					else if (Direction<Inexact>(dir).counterclockwise_in_between(
						Direction<Inexact>(directions[i]),
						Direction<Inexact>(directions[in])
					)) {
						return std::pair<int, int>(i, in);
					}
					else {
						i++;
						in++;
						if (in >= ndirs) {
							in = 0;
						}
					}
				}

				assert(false);
				return std::pair<int, int>(-1, -1);
			}

			void determine_significant_vertices() {

				if (progress) {
					(*progress)("1/7: Significant vertices", 0, graph.number_of_vertices());
				}

				for (Vertex_handle v : graph.vertices()) {

					if (progress) {
						(*progress)("", v->graph_index(), 0);
					}

					size_t d = v->degree();
					if (d >= 2) {
						std::pair<int, int> assoc_prev = associated_directions(v, v->incident_edge(d - 1));

						for (Edge_handle e : v->incident_edges()) {
							std::pair<int, int> assoc = associated_directions(v, e);

							if (assoc.first == assoc_prev.first
								|| (assoc.second >= 0 && assoc.second == assoc_prev.second)
								|| assoc.second == assoc_prev.first
								|| assoc.first == assoc_prev.second) {
								significant_vertices[v] = true;
								break;
							}

							assoc_prev = assoc;
						}
					}
				}

			};

			void subdivide_edges(Num lambda) {

				if (progress) {
					(*progress)("2/7: Subdividing edges", 0, graph.number_of_edges());
				}

				Num max_sqr_len = 0;
				for (Edge_handle e : graph.edges()) {
					Num sqr_len = approxSquaredLength(e);
					if (max_sqr_len < sqr_len) {
						max_sqr_len = sqr_len;
					}
				}

				max_sqr_len *= lambda * lambda;

				size_t cnt = graph.number_of_edges();
				for (size_t i = 0; i < cnt; i++) {
					Edge_handle e = graph.edge(i);

					if (progress) {
						(*progress)("", i, 0);
					}

					Num sqr_len = approxSquaredLength(e);

					int steps = (int)std::ceil(std::sqrt(sqr_len / max_sqr_len));
					if (steps < 2
						&& significant_vertices[e->source()]
						&& significant_vertices[e->target()]) {
						steps = 2;
					}

					Vector<Kernel> step = (e->target()->point() - e->source()->point()) / steps;
					Point<Kernel> pt = e->source()->point() + step;
					while (steps > 1) {
						e = graph.subdivide_edge(e, pt)->outgoing();

						pt = pt + step;
						steps--;
					}

				}
			};

			void find_best(std::vector<Vec>& out, std::vector<int>& build, Num build_cost, int next, std::vector<int>* best, Num* best_cost) {
				if (build_cost >= *best_cost) {
					// cannot get a better result
					return;
				}
				if (next == out.size()) {
					// handled all, this must outperform best
					*best_cost = build_cost;
					std::copy(build.begin(), build.end(), best->begin());
					return;
				}

				for (int i = 0; i < directions.size(); i++) {
					// test all unused directions for "next"
					int j = 0;
					while (j < next && build[j] != i) {
						j++;
					}
					if (j != next) continue;

					// compute the deviation
					Num dev = std::acos(out[next] * directions[i]);
					dev *= dev;

					// recurse
					build[next] = i;
					find_best(out, build, build_cost + dev, next + 1, best, best_cost);
				}

			}

			bool same_sector(Vertex_handle v, Edge_handle e, std::pair<int, int> assoc) {
				if (assoc.second < 0)
					return false;

				Vec dir = direction(v, e);

				return !is_approximately(directions[assoc.first], dir)
					&& !is_approximately(directions[assoc.second], dir)
					&& Direction<Inexact>(dir).counterclockwise_in_between(
						Direction<Inexact>(directions[assoc.first]),
						Direction<Inexact>(directions[assoc.second])
					);
			}

			void assign_directions() {


				size_t e_cnt = graph.number_of_edges();
				edge_data.resize(e_cnt);
				for (Edge_handle e : graph.edges()) {
					edge_data[e].edge = e;
				}

				// NB: subdivide may have increased vertex count
				// but the new vertices are insignificant by construction
				size_t v_cnt = significant_vertices.size();

				if (progress) {
					(*progress)("3/7: Assigning directions", 0, v_cnt);
				}

				for (int i = 0; i < v_cnt; i++) {

					if (progress) {
						(*progress)("", i, 0);
					}

					Vertex_handle v = graph.vertex(i);

					if (!significant_vertices[v]) continue;

					size_t degree = v->degree();
					std::vector<Vec> out;

					for (Edge_handle e : v->incident_edges()) {
						out.push_back(direction(v, e));
					}

					// TODO: optimize this, brute forced currently
					std::vector<int> best(degree, 0);
					Num best_cost = degree * 10;
					std::vector<int> build(degree);
					find_best(out, build, 0, 0, &best, &best_cost);

					for (size_t i = 0; i < degree; ++i) {

						Edge_handle e = v->incident_edge(i);
						EdgeData& edata = edge_data[e];
						edata.significant = v;
						edata.associated = associated_directions(v, e);
						edata.assigned = best[i];

						if (edata.associated.second < 0) {
							if (best[i] == edata.associated.first) {
								edata.type = EdgeType::ALIGN;
							}
							else {
								edata.type = EdgeType::DEV_ALIGN;
							}
						}
						else if (best[i] == edata.associated.first || best[i] == edata.associated.second) {

							bool has_same = false;
							for (size_t j = 0; j < degree; ++j) {
								if (i != j && same_sector(v, v->incident_edge(j), edata.associated)) {
									has_same = true;
									break;
								}
							}

							if (has_same) {
								edata.type = EdgeType::EVADING;
							}
							else {
								edata.type = EdgeType::UNALIGN;
							}
						}
						else {
							edata.type = EdgeType::DEV_UNALIGN;
						}
					}

				}

			};

			void assign_double_insignificant() {


				if (progress) {
					(*progress)("4/7: Completing directions", 0, graph.number_of_edges());
				}

				for (Edge_handle e : graph.edges()) {

					EdgeData& edata = edge_data[e];
					if (edata.type != EdgeType::UNDETERMINED) {
						continue;
					}

					if (progress) {
						(*progress)("", e->graph_index(), 0);
					}

					edata.significant = e->source();
					edata.associated = associated_directions(e->source(), e);
					if (edata.associated.second < 0) {
						edata.type = EdgeType::ALIGN;
						edata.assigned = edata.associated.first;
					}
					else {
						edata.type = EdgeType::UNALIGN;
						Vec dir = approximate(e->target()->point() - e->source()->point());
						Num det1 = CGAL::determinant(directions[edata.associated.first], dir);
						Num det2 = CGAL::determinant(dir, directions[edata.associated.second]);
						edata.assigned = det1 < det2 ? edata.associated.first : edata.associated.second;
					}
				}
			};

			template<typename K>
			std::pair<Number<K>, Number<K>> solve_vector_addition(Vector<K> d1, Vector<K> d2, Vector<K> s) {
				// find to factors, f and g, such that
				// f * d1 + g * d2 = s

				Number<K> sx = s.x();
				Number<K> sy = s.y();
				Number<K> d1x = d1.x();
				Number<K> d1y = d1.y();
				Number<K> d2x = d2.x();
				Number<K> d2y = d2.y();

				// e.g. find f and g such that:
				// sdx = f * d1x + g * d2x
				// sdy = f * d1y + g * d2y
				Number<K> f, g;

				// g * d2y = sy - f * d1y
				// NB division by 0!
				if (-0.000001 < d2y && d2y < 0.000001) {
					// d2y == 0
					// f * d1y = sy
					// NB: d1y cant be zero, as that would result in two parallel directions

					f = sy / d1y;

					// now, solve for g, using f, knowing that d2x != 0
					// sx = f * d1x + g * d2x hence yields:
					g = (sx - f * d1x) / d2x;
				}
				else {
					// g = (sy - f * d1y) / d2y
					// subsituting in first equation:
					// sx = f * d1x + ((sy - f * d1y) / d2y) * d2x
					// rewriting:
					// f (d1x - d1y * (d2x / d2y)) = sx - (sy / d2y) * d2x

					// can d1x - d1y * (d2x / d2y) = 0 hold?
					// == d2x / d2y = d1x / d1y
					// hence, this will only be zero if the directions d1 and d2 are
					// equal (which cannot occur)
					f = (sx - (sy / d2y) * d2x) / (d1x - d1y * (d2x / d2y));

					// now, solve for g, using f, knowing that d2y != 0
					// sy = f * d1y + g * d2y hence yields:
					g = (sy - f * d1y) / d2y;
				}

				return std::pair(f, g);
			}

			void determine_interference_regions(Num eps) {

				if (progress) {
					(*progress)("5/7: Interference regions", 0, graph.number_of_edges());
				}

				for (Edge_handle e : graph.edges()) {

					EdgeData& edata = edge_data[e];

					if (progress) {
						(*progress)("", e->graph_index(), 0);
					}

					switch (edata.type) {
					default: {
						assert(false); // undetermined should not occur
						break;
					}
					case EdgeType::ALIGN: {
						edata.interference.push_back(approximate(edata.significant->point()));
						edata.interference.push_back(approximate(e->other(edata.significant)->point()));
						break;
					}
					case EdgeType::UNALIGN:
					case EdgeType::EVADING: {
						Vec dir_1 = directions[edata.assigned];
						Vec dir_2 = directions[edata.other_direction()];
						Vec evec = vector(edata.significant, e);

						std::pair<Num, Num> fg = solve_vector_addition(dir_1, dir_2, evec);

						auto v_pt = approximate(edata.significant->point());
						auto o_pt = approximate(e->other(edata.significant)->point());
						edata.interference.push_back(v_pt);
						edata.interference.push_back(v_pt + dir_1 * fg.first);
						edata.interference.push_back(o_pt);
						edata.interference.push_back(v_pt + dir_2 * fg.second);

						if (!edata.interference.is_counterclockwise_oriented()) {
							edata.interference.reverse_orientation();
						}
						break;
					}
					case EdgeType::DEV_UNALIGN: {

						Vec dir_close, dir_far;
						if ((edata.assigned + 1) % directions.size() == edata.associated.first) {
							dir_close = directions[edata.associated.first];
							dir_far = directions[edata.associated.second];
						}
						else {
							dir_close = directions[edata.associated.second];
							dir_far = directions[edata.associated.first];
						}

						Vec step = vector(edata.significant, e) / 3.0;

						std::pair<Num, Num> fg = solve_vector_addition(dir_close, dir_far, step);

						dir_close *= fg.first;
						dir_far *= fg.second;

						Num halfcross = 0.5 * std::abs(CGAL::determinant(dir_close, dir_far));
						Num len = halfcross / std::abs(CGAL::determinant(dir_close, directions[edata.assigned]));
						Vec d = len * directions[edata.assigned];

						auto v_pt = approximate(edata.significant->point());
						auto o_pt = approximate(e->other(edata.significant)->point());

						auto far_pt = v_pt + d + dir_close;
						auto is = CGAL::intersection(Line<Inexact>(o_pt, o_pt + dir_far), Line<Inexact>(far_pt, far_pt + dir_close));
						auto corner = std::get<Point<Inexact>>(*is);

						edata.interference.push_back(v_pt);
						edata.interference.push_back(v_pt + d);
						if (CGAL::squared_distance(v_pt, far_pt) > CGAL::squared_distance(v_pt, corner)) {
							edata.interference.push_back(far_pt);
						}
						else {
							edata.interference.push_back(corner);
						}
						edata.interference.push_back(o_pt);
						edata.interference.push_back(o_pt - 3 * dir_close);

						if (!edata.interference.is_counterclockwise_oriented()) {
							edata.interference.reverse_orientation();
						}
						break;
					}
					case EdgeType::DEV_ALIGN: {

						auto v_pt = approximate(edata.significant->point());
						auto o_pt = approximate(e->other(edata.significant)->point());

						auto len = std::sqrt(CGAL::squared_distance(v_pt, o_pt));

						Vec dir_1 = directions[edata.assigned];
						Vec dir_2 = directions[edata.other_direction()];

						auto step_dist = len * eps;
						auto step_len = (1 - eps) * len / 2.0;

						auto sidestep = step_dist * dir_1;
						auto lengthstep = step_len * dir_2;

						auto a = v_pt + sidestep;
						auto b = a + lengthstep;
						auto c = b - 2 * sidestep;
						auto d = c + lengthstep;

						edata.interference.push_back(v_pt);
						edata.interference.push_back(a);
						edata.interference.push_back(b);
						edata.interference.push_back(o_pt);
						edata.interference.push_back(d); // yes d first, then c (order of polyine vs CH order)		
						edata.interference.push_back(c);

						if (!edata.interference.is_counterclockwise_oriented()) {
							edata.interference.reverse_orientation();
						}
						break;
					}
					}

					interference_tree.insert(&edata);
				}
			}

			void assign_step_counts(Num eps) {

				size_t e_cnt = graph.number_of_edges();

				VertexQuadTree<Graph> deg_zero_vertices(convert_kernel<Kernel>(interference_tree.root_box()), 10);
				for (Vertex_handle v : graph.vertices()) {
					if (v->degree() == 0) {
						deg_zero_vertices.insert(v);
					}
				}

				auto even_rounding = [](Num v) {
					int k = (int)std::ceil(v + 0.000001);
					k += k % 2;
					return k;
					};

				auto uncommon_interference = [](EdgeData& a, EdgeData& b) {
					if (a.type == EdgeType::ALIGN) {
						assert(b.type != EdgeType::ALIGN);


						if (!b.interference.has_on_unbounded_side(a.interference[0])) {
							return true;
						}

						Segment<Inexact> a_ls = Segment<Inexact>(a.interference[0], a.interference[1]);

						for (auto ls : b.interference.edges()) {
							if (CGAL::intersection(ls, a_ls).has_value()) {
								return true;
							}
						}
						return false;
					}
					else if (b.type == EdgeType::ALIGN) {
						assert(a.type != EdgeType::ALIGN);

						if (!a.interference.has_on_unbounded_side(b.interference[0])) {
							return true;
						}

						Segment<Inexact> b_ls = Segment<Inexact>(b.interference[0], b.interference[1]);
						for (auto ls : a.interference.edges()) {
							if (CGAL::intersection(ls, b_ls).has_value()) {
								return true;
							}
						}
						return false;
					}
					else {
						// this only works if the interference polygons are not segments (i.e. for aligned edges...)

						if (!a.interference.has_on_unbounded_side(b.interference[0])) {
							return true;
						}

						if (!b.interference.has_on_unbounded_side(a.interference[0])) {
							return true;
						}

						for (auto ls_a : a.interference.edges()) {
							for (auto ls_b : b.interference.edges()) {
								if (CGAL::intersection(ls_a, ls_b).has_value()) {
									return true;
								}
							}
						}

						return false;
						// this throws a CGAL error -- possibly because of inexact. Above is fine, our polygons are small...
						/*try {
							return CGAL::do_intersect(a.interference, b.interference);
						}
						catch (...) {
							std::cout << "ERROR\n";
							std::cout << "<ipeselection pos = \"0 0\">\n";
							std::cout << "<path layer=\"alpha\" stroke=\"red\" cap=\"1\" join=\"1\">\n";

							std::cout << a.interference[0] << " m\n";
							for (int i = 1; i < a.interference.size(); i++)
								std::cout << a.interference[i] << " l\n";
							std::cout << "h\n";
							std::cout << "</path>\n";
							std::cout << "<path layer=\"alpha\" stroke=\"blue\" cap=\"1\" join=\"1\">\n";

							std::cout << b.interference[0] << " m\n";
							for (int i = 1; i < b.interference.size(); i++)
								std::cout << b.interference[i] << " l\n";
							std::cout << "h\n";
							std::cout << "</path>\n";
							std::cout << "</ipeselection>";

							assert(false);
						}*/
					}
					};

				auto common_interference = [](EdgeData& a, EdgeData& b) {
					if (a.type == EdgeType::ALIGN || b.type == EdgeType::ALIGN) {
						// aligned cannot interfere with anything else
						return false;
					}

					if ((a.type == EdgeType::DEV_ALIGN || a.type == EdgeType::DEV_UNALIGN)
						&& (a.type == EdgeType::DEV_ALIGN || a.type == EdgeType::DEV_UNALIGN)) {
						// deviating cannot interfere with eachother
						return false;
					}

					if ((a.type == EdgeType::UNALIGN && b.type == EdgeType::UNALIGN)
						|| (a.type == EdgeType::UNALIGN && b.type == EdgeType::EVADING)
						|| (a.type == EdgeType::EVADING && b.type == EdgeType::UNALIGN)) {
						// unaligned cannot interfere with other unaligned or evading
						return false;
					}

					if (a.type == EdgeType::EVADING && b.type == EdgeType::EVADING) {
						// two evading edges only interfere if they lie in the same sector
						return a.associated.first == b.associated.first && a.associated.second == b.associated.second;
					}

					if (a.type == EdgeType::UNALIGN || a.type == EdgeType::EVADING) {
						if (b.type == EdgeType::DEV_UNALIGN) {
							// only interfere if they lie in the same sector
							return a.associated.first == b.associated.first && a.associated.second == b.associated.second;
						}
						else {
							assert(b.type == EdgeType::DEV_ALIGN);
							// count the aligned direction for both sectors
							return a.associated.first == b.associated.first || a.associated.second == b.associated.first;
						}
					}

					if (b.type == EdgeType::UNALIGN || b.type == EdgeType::EVADING) {
						if (a.type == EdgeType::DEV_UNALIGN) {
							// only interfere if they lie in the same sector
							return a.associated.first == b.associated.first && a.associated.second == b.associated.second;
						}
						else {
							assert(a.type == EdgeType::DEV_ALIGN);
							// count the aligned direction for both sectors
							return a.associated.first == b.associated.first || a.associated.first == b.associated.second;
						}
					}

					assert(false); // should not occur, the above cases are exhaustive
					return false;
					};

				auto ignore = [](Segment<Inexact> e, Num frac) {
					return Segment<Inexact>(e.start() + (e.end() - e.start()) * frac, e.end());
					};

				int handled = 0;

				if (progress) {
					(*progress)("6/7: Step counts", 0, e_cnt);
				}

				// first we do deviating edges
				for (size_t i = 0; i < e_cnt; i++) {

					Edge_handle e = graph.edge(i);
					EdgeData& edata = edge_data[e];

					switch (edata.type) {
					default: {
						assert(false); // shouldn't happen, undetermined edge somewhere?
						break;
					}
					case EdgeType::ALIGN:
					case EdgeType::UNALIGN:
					case EdgeType::EVADING: {
						break;
					}
					case EdgeType::DEV_UNALIGN: {
						if (progress) {
							(*progress)("", handled++, 0);
						}

						Num min_dist_sqr = std::numeric_limits<Num>::infinity();
						Segment<Inexact> seg = directed_segment(e, edata);
						Segment<Inexact> seg_ignore = ignore(seg, 1.0 / 3.0);

						Rectangle<Inexact> box = edata.box();
						interference_tree.findOverlapped(box, [&](EdgeData* other) {
							if (other->edge == e) {
								return;
							}

							Vertex_handle common = other->edge->common_vertex(e);

							if (common == nullptr) {
								// no shared vertex, treat normally

								if (!uncommon_interference(edata, *other)) {
									return;
								}

								min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, approximate(other->edge->curve())));
							}
							else if (common == edata.significant) {
								// the other edge must be also evading or deviating

								if (!common_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = directed_segment(other->edge, *other);
								min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg_ignore, seg_other));

							} // else: sharing an insignficant vertex, these staircases do not interact

							});

						Rectangle<Kernel> kbox = convert_kernel<Kernel>(box);
						deg_zero_vertices.findContained(kbox, [&](Vertex_handle vv) {
							min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(approximate(vv->point()), seg));
							});

						if (!std::isfinite(min_dist_sqr)) {
							edata.stepcount = 4;
						}
						else {
							Num min_dist = std::sqrt(min_dist_sqr);

							Vec dir_close, dir_far;
							if ((edata.assigned + 1) % directions.size() == edata.associated.first) {
								dir_close = directions[edata.associated.first];
								dir_far = directions[edata.associated.second];
							}
							else {
								dir_close = directions[edata.associated.second];
								dir_far = directions[edata.associated.first];
							}

							Num elen = std::sqrt(approxSquaredLength(e));

							// set dir_close/far such that describe a step of unit length
							Vec step = vector(edata.significant, e) / elen;
							std::pair<Num, Num> fg = solve_vector_addition(dir_close, dir_far, step);
							dir_close *= fg.first;
							dir_far *= fg.second;

							// the area of the unit step is the A = |dir_close X dir_far| / 2
							// we need to find a vector side_step = f D, such that |dir_close X side_step| = A
							//    where D is the assigned direction and some f > 0
							// that is, |dir_close X f D| = A
							// which is the same as f |dir_close X D| = A
							// so this reduces to f = A / |dir_close X D| = 0.5 |dir_close X dir_far| / |dir_close X D|

							Vec dir_assigned = directions[edata.assigned];
							Num f = 0.5 * std::abs(CGAL::determinant(dir_close, dir_far)) / std::abs(CGAL::determinant(dir_close, dir_assigned));
							Vec side_step = f * dir_assigned;

							Vec vp = dir_close + side_step;

							Num d_1 = std::sqrt(CGAL::squared_distance(seg.start() + vp, seg));

							edata.stepcount = std::max(4, even_rounding(2 * d_1 * elen / min_dist + 1));
						}
						break;
					}
					case EdgeType::DEV_ALIGN: {
						if (progress) {
							(*progress)("", handled++, 0);
						}
						Num min_dist_sqr = std::numeric_limits<Num>::infinity();
						Segment<Inexact> seg = directed_segment(e, edata);
						Segment<Inexact> seg_ignore = ignore(seg, (1 - eps) / 2.0);

						Rectangle<Inexact> box = edata.box();
						interference_tree.findOverlapped(box, [&](EdgeData* other) {

							if (other->edge == e) {
								return;
							}

							Vertex_handle common = other->edge->common_vertex(e);

							if (common == nullptr) {
								// no shared vertex, treat normally

								if (!uncommon_interference(edata, *other)) {
									return;
								}

								min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, approximate(other->edge->curve())));
							}
							else if (common == edata.significant) {
								// the other edge must be also evading or deviating

								if (!common_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = directed_segment(other->edge, *other);
								min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg_ignore, seg_other));

							} // else: sharing an insignficant vertex, these staircases do not interact

							});

						Rectangle<Kernel> kbox = convert_kernel<Kernel>(box);
						deg_zero_vertices.findContained(kbox, [&](Vertex_handle vv) {
							min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(approximate(vv->point()), seg));
							});

						Num maxdist = eps * std::sqrt(approxSquaredLength(e));
						if (!std::isfinite(min_dist_sqr)) {
							edata.stepcount = maxdist;
						}
						else {
							Num min_dist = std::sqrt(min_dist_sqr);
							edata.stepdist = std::min(min_dist / 2, maxdist);
						}
						break;
					}

					} // switch
				} // loop

				// then the remaining edges
				for (size_t i = 0; i < e_cnt; i++) {

					Edge_handle e = graph.edge(i);
					EdgeData& edata = edge_data[e];

					switch (edata.type) {
					default: {
						assert(false); // shouldn't happen, undetermined edge somewhere?
						break;
					}
					case EdgeType::ALIGN: {
						// nothing to do
						if (progress) {
							(*progress)("", handled++, 0);
						}
						break;
					}
					case EdgeType::UNALIGN: {
						if (progress) {
							(*progress)("", handled++, 0);
						}
						Num min_dist_sqr = std::numeric_limits<Num>::infinity();
						Segment<Inexact> seg = approximate(e->curve());

						Rectangle<Inexact> box = edata.box();
						interference_tree.findOverlapped(box, [&](EdgeData* other) {

							if (other->edge == e) {
								return;
							}

							Vertex_handle common = other->edge->common_vertex(e);

							if (common == nullptr) {
								// no shared vertex, treat normally

								if (!uncommon_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = approximate(other->edge->curve());
								Num dist_sqr = CGAL::squared_distance(seg, seg_other);
								min_dist_sqr = std::min(min_dist_sqr, dist_sqr);
							}
							else if (common == edata.significant) {
								// the other edge must be deviating, either aligned or unaligned

								if (!common_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = directed_segment(other->edge, *other);

								switch (other->type) {
								case EdgeType::DEV_ALIGN: {
									min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, ignore(seg_other, (1 - eps) / 2.0)));
									break;
								}
								case EdgeType::DEV_UNALIGN: {
									min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, ignore(seg_other, 1.0 / (other->stepcount - 1))));
									break;
								}
								default: {
									assert(false); // unexpected edge type
									break;
								}
								}

							} // else: sharing an insignficant vertex, these staircases do not interact
							});

						Rectangle<Kernel> kbox = convert_kernel<Kernel>(box);
						deg_zero_vertices.findContained(kbox, [&](Vertex_handle vv) {
							min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(approximate(vv->point()), seg));
							});

						if (!std::isfinite(min_dist_sqr)) {
							edata.stepcount = 2;
						}
						else {
							Num min_dist = std::sqrt(min_dist_sqr);
							Vec dir_1 = directions[edata.assigned];
							Vec dir_2 = directions[edata.other_direction()];
							Num elen = std::sqrt(approxSquaredLength(e));
							Vec edir = vector(edata.significant, e) / elen;
							Num alpha_1 = std::abs(std::acos(edir * dir_1));
							Num alpha_2 = std::abs(std::acos(edir * dir_2));
							Num lmax = (1 / std::tan(alpha_1) + 1 / std::tan(alpha_2)) * min_dist / 2;
							edata.stepcount = even_rounding(elen / lmax);
						}
						break;
					}
					case EdgeType::EVADING: {
						if (progress) {
							(*progress)("", handled++, 0);
						}
						Num min_dist_sqr = std::numeric_limits<Num>::infinity();
						Segment<Inexact> seg = directed_segment(e, edata);
						Segment<Inexact> seg_ev = ignore(seg, 0.5);

						Rectangle<Inexact> box = edata.box();
						interference_tree.findOverlapped(box, [&](EdgeData* other) {

							if (other->edge == e) {
								return;
							}

							Vertex_handle common = other->edge->common_vertex(e);

							if (common == nullptr) {
								// no shared vertex, treat normally

								if (!uncommon_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = approximate(other->edge->curve());
								Num dist_sqr = CGAL::squared_distance(seg, seg_other);
								min_dist_sqr = std::min(min_dist_sqr, dist_sqr);
							}
							else if (common == edata.significant) {
								// the other edge must be also evading or deviating

								if (!common_interference(edata, *other)) {
									return;
								}

								Segment<Inexact> seg_other = directed_segment(other->edge, *other);

								switch (other->type) {
								case EdgeType::EVADING: {
									min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg_ev, seg_other));
									break;
								}
								case EdgeType::DEV_ALIGN: {
									min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, ignore(seg_other, (1 - eps) / 2.0)));
									break;
								}
								case EdgeType::DEV_UNALIGN: {
									min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(seg, ignore(seg_other, 1.0 / (other->stepcount - 1))));
									break;
								}
								default: {
									assert(false); // unexpected edge type
									break;
								}
								}

							} // else: sharing an insignficant vertex, these staircases do not interact

							});

						Rectangle<Kernel> kbox = convert_kernel<Kernel>(box);
						deg_zero_vertices.findContained(kbox, [&](Vertex_handle vv) {
							min_dist_sqr = std::min(min_dist_sqr, CGAL::squared_distance(approximate(vv->point()), seg));
							});

						if (!std::isfinite(min_dist_sqr)) {
							edata.stepcount = 2;
						}
						else {
							Num min_dist = std::sqrt(min_dist_sqr);
							Vec dir_1 = directions[edata.assigned];
							Vec dir_2 = directions[edata.other_direction()];
							Num elen = std::sqrt(approxSquaredLength(e));
							Vec edir = vector(edata.significant, e) / elen;
							Num alpha_1 = std::abs(std::acos(edir * dir_1));
							Num alpha_2 = std::abs(std::acos(edir * dir_2));
							Num lmax = (1 / std::tan(alpha_1) + 1 / std::tan(alpha_2)) * min_dist / 2;
							edata.stepcount = even_rounding(elen / lmax);
						}
						break;
					}
					case EdgeType::DEV_UNALIGN:
					case EdgeType::DEV_ALIGN: {
						// already done in previous loop
						break;
					}

					} // switch
				} // loop

			};

			void create_staircases(Num eps) {

				size_t e_cnt = graph.number_of_edges();
				if (progress) {
					(*progress)("7/7: Making staircases", 0, e_cnt);
				}

				for (size_t i = 0; i < e_cnt; i++) {

					if (progress) {
						(*progress)("", i, 0);
					}

					Edge_handle e = graph.edge(i);
					EdgeData& edata = edge_data[e];
					Vertex_handle v = edata.significant;

					std::function<Edge_handle(Edge_handle, Point<Kernel>)> split;
					if (v == e->source()) {
						split = [this](Edge_handle edge, Point<Kernel> point) {
							return graph.subdivide_edge(edge, point)->outgoing();
							};
					}
					else {
						split = [this](Edge_handle edge, Point<Kernel> point) {
							return graph.subdivide_edge(edge, point)->incoming();
							};
					}

					switch (edata.type) {
					case EdgeType::ALIGN: {
						// dont need to do anything
						break;
					}
					case EdgeType::UNALIGN: {
						int k = edata.stepcount;

						assert(k > 1 && k % 2 == 0);

						Vector<Kernel> d1 = convert_kernel<Kernel>(directions[edata.assigned]);
						Vector<Kernel> d2 = convert_kernel<Kernel>(directions[edata.other_direction()]);

						if (e->source() != v) {
							d1 *= -1;
							d2 *= -1;
						}

						const Point<Kernel>& s = e->source()->point();
						const Point<Kernel>& t = e->target()->point();
						Vector<Kernel> step = (t - s) / k;

						std::pair<Number<Kernel>, Number<Kernel>> fg = solve_vector_addition(d1, d2, step);

						d1 = fg.first * d1;
						d2 = fg.second * d2;

						Point<Kernel> pt = s + d1;
						e = graph.subdivide_edge(e, pt)->outgoing();
						pt += 2 * d2;
						e = graph.subdivide_edge(e, pt)->outgoing();
						k -= 2;
						while (k > 0) {
							pt += 2 * d1;
							e = graph.subdivide_edge(e, pt)->outgoing();
							pt += 2 * d2;
							e = graph.subdivide_edge(e, pt)->outgoing();
							k -= 2;
						}

						break;
					}
					case EdgeType::EVADING: {
						int k = edata.stepcount;

						assert(k > 1 && k % 2 == 0);

						Vector<Kernel> d1 = convert_kernel<Kernel>(directions[edata.assigned]);
						Vector<Kernel> d2 = convert_kernel<Kernel>(directions[edata.other_direction()]);

						if (e->source() != v) {
							d1 *= -1;
							d2 *= -1;
						}

						const Point<Kernel>& s = e->source()->point();
						const Point<Kernel>& t = e->target()->point();
						Vector<Kernel> step = (t - s) / k;

						std::pair<Number<Kernel>, Number<Kernel>> fg = solve_vector_addition(d1, d2, step);

						d1 = fg.first * d1;
						d2 = fg.second * d2;

						// outward
						Point<Kernel> pt_s = s + d1;
						e = graph.subdivide_edge(e, pt_s)->outgoing();
						Point<Kernel> pt_t = t - d1;
						e = graph.subdivide_edge(e, pt_t)->incoming();
						k -= 2;
						while (k > 0) {
							// back to central line
							pt_s += d2;
							e = graph.subdivide_edge(e, pt_s)->outgoing();
							pt_t -= d2;
							e = graph.subdivide_edge(e, pt_t)->incoming();

							// step outward
							pt_s += d1;
							e = graph.subdivide_edge(e, pt_s)->outgoing();
							pt_t -= d1;
							e = graph.subdivide_edge(e, pt_t)->incoming();
							k -= 2;
						}
						break;
					}
					case EdgeType::DEV_UNALIGN: {
						int k = edata.stepcount;

						assert(k >= 4 && k % 2 == 0);

						Vec dir_close, dir_far;
						if ((edata.assigned + 1) % directions.size() == edata.associated.first) {
							dir_close = directions[edata.associated.first];
							dir_far = directions[edata.associated.second];
						}
						else {
							dir_close = directions[edata.associated.second];
							dir_far = directions[edata.associated.first];
						}

						// set dir_close/far such that describe a step of unit length
						Vec step = vector(edata.significant, e) / (k - 1);

						std::pair<Num, Num> fg = solve_vector_addition(dir_close, dir_far, step);
						dir_close *= fg.first;
						dir_far *= fg.second;

						// the area of the unit step is the A = |dir_close X dir_far| / 2
						// we need to find a vector side_step = f D, such that |dir_close X side_step| = A
						//    where D is the assigned direction and some f > 0
						// that is, |dir_close X f D| = A
						// which is the same as f |dir_close X D| = A
						// so this reduces to f = A / |dir_close X D| = 0.5 |dir_close X dir_far| / |dir_close X D|

						Vec dir_assigned = directions[edata.assigned];
						Num f = 0.5 * std::abs(CGAL::determinant(dir_close, dir_far)) / std::abs(CGAL::determinant(dir_close, dir_assigned));

						Vector<Kernel> side_step = convert_kernel<Kernel>(f * dir_assigned);
						Vector<Kernel> close_step = convert_kernel<Kernel>(dir_close);
						Vector<Kernel> far_step = convert_kernel<Kernel>(dir_far);

						int hk = k / 2;

						// side step
						auto pt = v->point() + side_step;
						e = split(e, pt);
						pt += close_step;
						e = split(e, pt);
						pt -= side_step;
						if (!is_approximately(-dir_assigned, dir_far / fg.second)) {
							e = split(e, pt);
						}
						hk -= 2;

						// evading steps
						while (hk > 0) {
							hk--;

							pt += far_step;
							e = split(e, pt);
							pt += close_step;
							e = split(e, pt);
						}

						// other steps
						hk = k / 2;
						pt += 2 * far_step;
						e = split(e, pt);
						hk--;
						while (hk > 0) {
							hk--;

							pt += close_step;
							e = split(e, pt);
							pt += far_step;
							e = split(e, pt);
						}
						break;
					}
					case EdgeType::DEV_ALIGN: {

						auto v_pt = approximate(v->point());
						auto o_pt = approximate(e->other(v)->point());

						auto len = std::sqrt(CGAL::squared_distance(v_pt, o_pt));

						Vec dir_1 = directions[edata.assigned];
						Vec dir_2 = directions[edata.associated.first];

						auto step_dist = edata.stepdist;
						auto step_len = (1 - eps) * len / 2.0; // NB: this needs to use eps, not stepdist

						auto sidestep = convert_kernel<Kernel>(step_dist * dir_1);
						auto lengthstep = convert_kernel<Kernel>(step_len * dir_2);

						auto a = v->point() + sidestep;
						auto b = a + lengthstep;
						auto c = b - 2 * sidestep;
						auto d = c + lengthstep;
						auto f = d + sidestep;

						e = split(e, a);
						e = split(e, b);
						e = split(e, c);
						e = split(e, d);
						split(e, f);
						break;
					}
					default:
						assert(false);// undetermined should not occur
						break;
					}

				}
			};
		};

		template<class Graph> requires Graph::Graph_traits::sorted&& Graph::Graph_traits::oriented
			void restrict_directions(Graph& graph, std::vector<Vector<Inexact>> directions, Number<Inexact> lambda, Number<Inexact> eps, std::optional<std::function<void(std::string, int, int)>> progress = std::nullopt) {

			assert(graph.is_initialized());

			// we assume directions are sorted by construction
			// and we assume they are normalized

			detail::AngleRestriction<Graph> ar(graph, directions, progress);

			ar.determine_significant_vertices();
			ar.subdivide_edges(lambda);
			ar.assign_directions();
			ar.assign_double_insignificant();
			ar.determine_interference_regions(eps);
			ar.assign_step_counts(eps);
			ar.create_staircases(eps);

			assert(graph.is_initialized());
		}

	} // namespace detail

	template<class Graph> requires Graph::Graph_traits::sorted&& Graph::Graph_traits::oriented
		void restrict_orientations(Graph& graph, std::vector<Vector<Inexact>> orientations, Number<Inexact> lambda, Number<Inexact> eps, std::optional<std::function<void(std::string, int, int)>> progress) {

		// convert orientations ("mod 180 degrees")
		// to directions ("mod 360 degrees")
		std::vector<Vector<Inexact>> dirs;

		for (Vector<Inexact> v : orientations) {
			assert(0.999999 <= v.squared_length() && v.squared_length() <= 1.000001);
			dirs.push_back(v);
		}
		for (Vector<Inexact> v : orientations) {
			dirs.push_back(-1 * v);
		}

		detail::restrict_directions(graph, dirs, lambda, eps, progress);
	}

	template<class Graph> requires Graph::Graph_traits::sorted&& Graph::Graph_traits::oriented
		void restrict_orientations(Graph& graph, int count, Number<Inexact> initial_angle, Number<Inexact> lambda, Number<Inexact> eps, std::optional<std::function<void(std::string, int, int)>> progress) {

		using Vec = Vector<Inexact>;
		using Transform = CGAL::Aff_transformation_2<Inexact>;

		assert(count >= 2);

		std::vector<Vec> dirs;

		// first direction
		Vec dir(1, 0);
		dir = dir.transform(Transform(CGAL::ROTATION, std::sin(initial_angle), std::cos(initial_angle)));
		dirs.push_back(dir);

		// the other directions
		Number<Inexact> a = std::numbers::pi / count;
		int i = 1;
		while (i < 2 * count) {

			Transform t(CGAL::ROTATION, std::sin(i * a), std::cos(i * a));
			dirs.push_back(dir.transform(t));
			i++;
		}

		detail::restrict_directions(graph, dirs, lambda, eps, progress);
	}

	template<class Graph> requires Graph::Graph_traits::sorted&& Graph::Graph_traits::oriented
		void restrict_orientations(Graph& graph, std::vector<Number<Inexact>> angles, Number<Inexact> lambda, Number<Inexact> eps, std::optional<std::function<void(std::string, int, int)>> progress) {

		using Vec = Vector<Inexact>;
		using Transform = CGAL::Aff_transformation_2<Inexact>;

		// right
		Vec dir(1, 0);

		std::vector<Vec> dirs;

		for (Number<Inexact> a : angles) {
			Transform t(CGAL::ROTATION, std::sin(a), std::cos(a));
			dirs.push_back(dir.transform(t));
		}

		for (int i = 0; i < angles.size(); i++) {
			dirs.push_back(-1 * dirs[i]);
		}

		detail::restrict_directions(graph, dirs, lambda, eps, progress);
	}

	template<class Graph> requires Graph::Graph_traits::sorted&& Graph::Graph_traits::oriented
		void restrict_orientations(Graph& graph, std::initializer_list<Number<Inexact>> angles, Number<Inexact> lambda, Number<Inexact> eps, std::optional<std::function<void(std::string, int, int)>> progress) {

		using Vec = Vector<Inexact>;
		using Transform = CGAL::Aff_transformation_2<Inexact>;

		// right
		Vec dir(1, 0);

		std::vector<Vec> dirs;

		for (Number<Inexact> a : angles) {
			Transform t(CGAL::ROTATION, std::sin(a), std::cos(a));
			dirs.push_back(dir.transform(t));
		}

		for (int i = 0; i < angles.size(); i++) {
			dirs.push_back(-1 * dirs[i]);
		}

		detail::restrict_directions(graph, dirs, lambda, eps, progress);
	}
}