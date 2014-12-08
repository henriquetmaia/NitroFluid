// -----------------------------------------------------------------------------
// libDDG -- Fluid.h
// -----------------------------------------------------------------------------
//
// Fluid handles fluid flow tools and operations.
//

#ifndef DDG_FLUID_H
#define DDG_FLUID_H

#include "Face.h"
#include "HalfEdge.h"
#include "Vector.h"
#include "Mesh.h"

namespace DDG
{

   class FLUID
   {
      public:

         enum AdvectionScheme
         {
            SEMI_LAGRANGIAN
         };

         enum Interpolation
         {
            WHITNEY
         };

         enum ProjectionComponent
         {
            CURL,
            DIV,
            HARMONIC
         };

         Fluid( Mesh*& surface_ptr );

         ~Fluid( void );

         static void flow( const float& dt,
                           const AdvectionScheme& advectType,
                           const ProjectionComponent& projectType,
                           const Interpolation& interpType 
                         );
         // step fluid velocity and markers foward by dt with choice of advection, projection, and interpolation

      protected:
         // static void prescribeField( ? );
         // // initates the fluid field somehow

         // static void interact( ? );
         // // updates/forces the field according to user interaction        

      private:

         static void advectMarkers( const float& dt ) = 0;

         static void advectSemiLagrangian( const float& dt );

         static void projectCurl( void ) = 0;

         static void projectDivergence( void );

         static void projectHarmonic( void ) = 0;

         static double rayEdgeIntersectDistance( const Vector& coordinate, const Vector& direction, const HalfEdge& half_edge );

         static Vector rotateAcrossBy( const Vector& direction, const Vector& axis, const double& angle );

         static void BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k );

         static void updateEdgeWeights( void );

         static Vector whitneyInterpolate( const Vector& coordinate, const Edge& edge );

         static Vector whitneyInterpolate( const Vector& coordinate, const HalfEdge& he );

         // viscosity/density/number of samples/other parameters?

         Mesh* fluid_ptr;
   };
}

#endif

