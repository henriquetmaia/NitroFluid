// -----------------------------------------------------------------------------
// libDDG -- Vertex.h
// -----------------------------------------------------------------------------
//
// Vertex stores attributes associated with a mesh edge.  The iterator he
// points to its "outgoing" halfedge.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_VERTEX_H
#define DDG_VERTEX_H

#include "Vector.h"
#include "Types.h"

namespace DDG
{
   class Vertex
   {
      public:
         Vector normal( void ) const;
         // returns the vertex normal

         bool isIsolated( void ) const;
         // returns true if the vertex is not contained in any face or edge; false otherwise

         int valence( void ) const;
         // returns the number of incident faces / edges

         void setID( const int& new_id ) { index = new_id; };
         // assigns a new ID to the vertex

         int getID( void ) const { return index; }

         double dualArea( void ) const;
         // one third the are of incident faces, area of dual polygon surrounding vertex

         HalfEdgeIter he;
         // points to the "outgoing" halfedge

         Vector position;
         // location of vertex in Euclidean 3-space

      private:
         int index;// = -1;
         // unique ID with reference to a mesh in range 0, ... , NumVertices - 1

         // double rho;
         // Marker value associated with the fluid at this position in space
         // Could represent color, density, temperature, or any other scalar.

   };
}

#endif

