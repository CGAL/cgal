/*
Copyright (C) 2026  TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------
 
#include <cassert>

namespace CGAL::internal {

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::siftUp(int k, Element_handle elt) {
		while (k > 0) {
			int parent = (k - 1) >> 1;
			Element_handle e = queue[parent];
			if (traits.compare(elt, e) >= 0) {
				break;
			}
			queue[k] = e;
		    traits.setIndex(e, k);
			k = parent;
		}
		queue[k] = elt;
	    traits.setIndex(elt, k);
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::siftDown(int k, Element_handle elt) {
		int half = queue.size() >> 1;
		while (k < half) {
			int child = (k << 1) + 1;
			Element_handle c = queue[child];
			int right = child + 1;
		    if (right < queue.size() && traits.compare(c, queue[right]) > 0) {
				c = queue[child = right];
			}
		    if (traits.compare(elt, c) <= 0) {
				break;
			}
			queue[k] = c;
		    traits.setIndex(c, k);
			k = child;
		}
		queue[k] = elt;
	    traits.setIndex(elt, k);
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::empty() const {
		return queue.empty();
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::push(Element_handle elt) {
		assert(!contains(elt));
		assert(std::find(queue.begin(), queue.end(), elt) == queue.end());

		queue.push_back(elt);
		siftUp(queue.size() - 1, elt);
	}

	template <QueueTraits QT>
	QT::Element_handle IndexedPriorityQueue<QT>::pop() {
		Element_handle result = queue[0];
	    traits.setIndex(result, -1);

		Element_handle last = queue[queue.size() - 1];
		queue.pop_back();
		if (!queue.empty()) {
			siftDown(0, last);
		}

		assert(!contains(result));
		assert(std::find(queue.begin(), queue.end(), result) == queue.end());

		return result;
	}

	template <QueueTraits QT>
	QT::Element_handle IndexedPriorityQueue<QT>::peek() const {
		if (queue.empty()) {
			return nullptr;
		}
		return queue[0];
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::remove(Element_handle elt) {
	    int id = traits.getIndex(elt);
		if (id < 0 || id >= queue.size() || queue[id] != elt) {
			assert(!contains(elt));
			assert(std::find(queue.begin(), queue.end(), elt) == queue.end());
			return false;
		}
		else {
		    traits.setIndex(elt, -1);
			if (id == queue.size() - 1) {
				queue.pop_back();
			}
			else {
				Element_handle moved = queue[queue.size() - 1];
				queue.pop_back();
				siftDown(id, moved);
				if (queue[id] == moved) {
					siftUp(id, moved);
				}
			}

			assert(!contains(elt));
			assert(std::find(queue.begin(), queue.end(), elt) == queue.end());
			return true;
		}
	}

	template <QueueTraits QT>
	bool IndexedPriorityQueue<QT>::contains(Element_handle elt) const {
	    int id = traits.getIndex(elt);
		if (id < 0 || id >= queue.size()) {
			return false;
		}
		else {
			return queue[id] == elt;
		}
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::update(Element_handle elt) {
		assert(contains(elt));
		assert(std::find(queue.begin(), queue.end(), elt) != queue.end());

		siftUp(traits.getIndex(elt), elt);
	    siftDown(traits.getIndex(elt), elt);
	}

	template <QueueTraits QT>
	void IndexedPriorityQueue<QT>::clear() {
		for (Element_handle elt : queue) {
		    traits.setIndex(elt, -1);
		}
		queue.clear();
	}

	template <QueueTraits QT> 
	const std::vector<typename QT::Element_handle>& IndexedPriorityQueue<QT>::content() const {
	    return queue;
    }
    } // namespace CGAL::internal