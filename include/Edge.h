// -----------------------------------------------------------------------------
// libDDG -- Edge.h
// -----------------------------------------------------------------------------
//
// Edge stores attributes associated with a mesh edge.  The iterator he points
// to one of its two associated halfedges.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_EDGE_H
#define DDG_EDGE_H

#include "Types.h"

namespace DDG
{
   class Edge
   {

      double getCoef( void ) const { return ref_coef; }

      void setCoef( double& newCoef ) { mod_coef = newCoef; }

      void updateRefCoef( void ) { ref_coef = mod_coef; }

      public:
         HalfEdgeIter he;
         // points to one of the two halfedges associated with this edge

      private:
   		double ref_coef;
   		// Integrated fluid velocity weight associated with flow across this edge
   		// READ ONLY Value from the start of the timestep used as reference for computing flow
            //TODO: Should this be (this potentially could be...) an array of weights for multiple velocity samples along edge?

   		double mod_coef;
   		 // WRITE ONLY Modified value for next timestep
   };
}

#endif
