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
#include "DiscreteExteriorCalculus.h"

namespace DDG
{

   class Fluid
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

         Fluid( Mesh*& surface_ptr, const ProjectionComponent& projectType = DIV );

         ~Fluid( void );

         // static void flow( const float& dt,
         //                   const AdvectionScheme& advectType,
         //                   const ProjectionComponent& projectType,
         //                   const Interpolation& interpType 
         //                 );
         // step fluid velocity and markers foward by dt with choice of advection, projection, and interpolation

      protected:
         // static void prescribeField( ? );
         // // initates the fluid field somehow

         // static void interact( ? );
         // // updates/forces the field according to user interaction        

      private:
         
         void buildOperators( );
         
         void advectMarkers( const float& dt );

         void advectSemiLagrangian( const float& dt );

         void projectCurl( void );

         void projectDivergence( void );

         void projectHarmonic( void );

         void updateEdgeWeights( void );

         // Static helper functions:
         static double intersectRay( const Vector& coordinate, const Vector& direction, const HalfEdgeIter& half_edge, const double tmax )

         static Vector rotateAcrossBy( const Vector& direction, const Vector& axis, const double& angle );

         static void BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k );

         static Vector whitneyInterpolate( const Vector& coordinate, const EdgeIter& edge );

         static Vector whitneyInterpolate( const Vector& coordinate, const HalfEdgeIter& he );

         // viscosity/density/number of samples/other parameters?

         Mesh* fluid_ptr;

         // Discrete Operators
         SparseMatrix<Real> star0;
         SparseMatrix<Real> star1;
         SparseMatrix<Real> star2;
         SparseMatrix<Real> d0;
         SparseMatrix<Real> d1;

   };
}

#endif

