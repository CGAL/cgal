#include "DeleteCurveCallback.h"
#include "Utils.h"

/*! Constructor */
template <typename Arr_>
DeleteCurveCallback<Arr_>::DeleteCurveCallback(
	Arrangement* arr_, QObject* parent_) :
	CGAL::Qt::Callback(parent_),
	scene(NULL), highlightedCurve(new CGAL::Qt::CurveGraphicsItem<Traits>()),
	arr(arr_), deleteOriginatingCurve(false)
{
	QObject::connect(
		this, SIGNAL(modelChanged()), this->highlightedCurve,
		SLOT(modelChanged()));
}

//! Setter Function.
/*!
  sets the current scene of the viewport
*/
template <typename Arr_>
void DeleteCurveCallback<Arr_>::setScene(QGraphicsScene* scene_)
{
	this->scene = scene_;
	this->highlightedCurve->setScene(scene_);
	if (this->scene)
	{
		this->scene->addItem(this->highlightedCurve);
	}
}

//! Getter function.
/*!
  returns the current scene of the viewport
*/
template <typename Arr_>
QGraphicsScene* DeleteCurveCallback<Arr_>::getScene() const
{
	return this->scene;
}

template <typename Arr_>
void DeleteCurveCallback<Arr_>::reset()
{
	this->highlightedCurve->clear();
	this->removableHalfedge = Halfedge_handle();
	this->deleteOriginatingCurve = false;
	Q_EMIT modelChanged();

	QGraphicsView* view = this->scene->views().first();
	view->scale(1.01, 1.01);
	view->scale(1 / 1.01, 1 / 1.01);
}

template <typename Arr_>
void DeleteCurveCallback<Arr_>::partialReset()
{
	this->highlightedCurve->clear();
	this->removableHalfedge = Halfedge_handle();
	Q_EMIT modelChanged();
}

template <typename Arr_>
void DeleteCurveCallback<Arr_>::mousePressEvent(
	QGraphicsSceneMouseEvent* /* event */)
{
	if (this->removableHalfedge == Halfedge_handle())
	{
		return;
	}

	if (this->deleteOriginatingCurve)
	{
		Originating_curve_iterator it =
			this->arr->originating_curves_begin(this->removableHalfedge);
		Originating_curve_iterator it_end =
			this->arr->originating_curves_end(this->removableHalfedge);
		while (it != it_end)
		{
			Originating_curve_iterator temp = it;
			++temp;
			CGAL::remove_curve(*(this->arr), it);
			it = temp;
		}
	}
	else
	{
		// CGAL::remove_edge( *(this->arr), this->removableHalfedge->curve( ) );
		this->arr->remove_edge(this->removableHalfedge);
	}

	this->partialReset();

	QGraphicsView* view = this->scene->views().first();
	view->scale(1.01, 1.01);
	view->scale(1 / 1.01, 1 / 1.01);
}

template <typename Arr_>
void DeleteCurveCallback<Arr_>::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
	this->highlightNearestCurve(event);
}

template <typename Arr_>
void DeleteCurveCallback<Arr_>::highlightNearestCurve(
	QGraphicsSceneMouseEvent* event)
{
	// find the nearest curve to the cursor to be the new highlighted curve
	Point p = this->convert(event->scenePos());

	Find_nearest_edge<Arr_> findNearestEdge(this->arr);
	findNearestEdge.setScene(this->scene);
	Halfedge_const_handle nearestEdge = findNearestEdge(p);
	this->removableHalfedge = this->arr->non_const_handle(nearestEdge);

	// now 'removableHalfedge' holds the closest halfedge to the point of the
	// mouse
	// this->removableHalfedge = nearestHei;
	// if ( isFirst )
	if (this->removableHalfedge == Halfedge_handle())
	{
		return;
	}

	// create a curve graphics item and add it to the scene
	this->highlightedCurve->clear();
	if (this->deleteOriginatingCurve)
	{ // highlight the originating curve
		Originating_curve_iterator ocit, temp;
		ocit = this->arr->originating_curves_begin(this->removableHalfedge);
		while (ocit !=
			   this->arr->originating_curves_end(this->removableHalfedge))
		{
			temp = ocit;
			++temp;

			Curve_handle ch = ocit;
			Induced_edge_iterator itr;
			for (itr = this->arr->induced_edges_begin(ch);
				 itr != this->arr->induced_edges_end(ch); ++itr)
			{
				X_monotone_curve_2 curve = (*itr)->curve();
				this->highlightedCurve->insert(curve);
			}
			ocit = temp;
		}
	}
	else
	{ // highlight just the edge
		this->highlightedCurve->insert(this->removableHalfedge->curve());
	}

	Q_EMIT modelChanged();
}

template class DeleteCurveCallback<Seg_arr>;
template class DeleteCurveCallback<Pol_arr>;
template class DeleteCurveCallback<Conic_arr>;
template class DeleteCurveCallback<Lin_arr>;
template class DeleteCurveCallback<Alg_seg_arr>;
