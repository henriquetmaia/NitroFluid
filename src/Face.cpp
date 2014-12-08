#include "Face.h"
#include "Mesh.h"
#include "Vector.h"

namespace DDG
{
   Vector Face::normal( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;

      return cross( p1-p0, p2-p0 ).unit();
   }

   bool Face::isBoundary( void ) const
   {
      return he->onBoundary;
   }

   double Face::area( void ) const
   {
      Vector e1 = he->vertex->position - he->next->next->vertex->position;
      Vector e2 = he->next->vertex->position - he->next->next->vertex->position;
      double theta = acos ( dot( e1, e2 ) / ( e1.norm() * e2.norm() ) ); 

      return 0.5 * e1.norm() * e2.norm() * sin( theta );
   }
}

