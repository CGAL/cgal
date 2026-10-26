#pragma once

namespace CGAL::internal {

	template<class Elt, typename Kernel>
	struct GraphQueueTraits {

		using Element_handle = Elt;

		static void setIndex(Element_handle elt, int index) {
			elt->data().queue_index = index;
		}

		static int getIndex(Element_handle elt) {
			return elt->data().queue_index;
		}

		static int compare(Element_handle a, Element_handle b) {
			Number<Kernel> ac = a->data().cost;
			Number<Kernel> bc = b->data().cost;
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
}