#include <vector>
using namespace std;

#include "Vertex.h"
#include "Mesh.h"
#include "HalfEdge.h"

namespace DDG
{
   Vector Vertex::normal( void ) const
   // returns the vertex normal
   {
      Vector N;

      HalfEdgeCIter h = he;
      do
      {
         N += h->face->normal();
         h = h->flip->next;
      }
      while( h != he );

      return N.unit();
   }

   vector<HalfEdge> isolated; // all isolated vertices point to isolated.begin()

   bool Vertex::isIsolated( void ) const
   // returns true if the vertex is not contained in any face or edge; false otherwise
   {
      return he == isolated.begin();
   }

   int Vertex::valence( void ) const
   // returns the number of incident faces
   {
      int n = 0;

      HalfEdgeCIter h = he;
      do
      {
         n++;
         h = h->flip->next;
      }
      while( h != he );

      return n;
   }

   void Vertex::setID( const int& new_id )
   {
      v_id = new_id;
   }

   double Vertex::dualArea( void ) const
   {
      double incidentAreaSum = 0;
      HalfEdgeCIter h = he;
      do
      {
         incidentAreaSum += h->face->area();
         h = h->flip->next;
      }
      while( h != he );

      return incidentAreaSum / 3;
   }

}

