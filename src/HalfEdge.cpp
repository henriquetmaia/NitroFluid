#include "HalfEdge.h"
#include "Mesh.h"

namespace DDG
{

	double HalfEdge::cotan() const
	{
		Vector e1 = vertex->position - next->next->vertex->position; // *vertex->position - next->next->vertex->position; 
		Vector e2 = next->vertex->position - next->next->vertex->position;

		double cosTheta = dot( e1, e2 ) / ( e1.norm() * e2.norm() );
		double sinTheta = sin( acos( cosTheta ) );

		return cosTheta/sinTheta;
	}

}

