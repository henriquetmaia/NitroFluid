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

         Fluid( Mesh* surface_ptr, const ProjectionComponent& projectType = DIV );

         ~Fluid( void );

         void advectColorAlongField( const float& dt );

         void advectVelocitySemiLagrangian( const float& dt );

         void projectCurl( void );

         void projectDivergence( void );

         void projectHarmonic( void );

         void updateEdgeWeights( void );

      protected:

         void prescribeVelocityField( int vf );
 
         void prescribeDensity( int d );

         // static void interact( ? );
         // // updates/forces the field according to user interaction        

      private:

         void buildOperators( );

         void backtrackAlongField( const float& dt, const Vector start_pt, const Vector& start_vel, HalfEdgeIter& current_he, Vector& final_coord, Quaternion& accumulated_angle );

         // Static helper functions:
         static double intersectRay( Vector& coordinate, const Vector& direction, const HalfEdgeIter& half_edge, const double tmax, HalfEdgeIter& crossing_half_edge );

         static Vector rotateAcrossBy( const Vector& direction, const Vector& axis, const double& angle );

         static void BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k );

         static Vector whitneyInterpolateVelocity( const Vector& coordinate, const EdgeIter& edge, HalfEdgeIter& chosen_he );

         static Vector whitneyInterpolateVelocity( const Vector& coordinate, const HalfEdgeIter& he );

         static Vector barycentricInterpolateColors( const Vector& coordinate, const HalfEdgeIter& half_edge );

         // viscosity/density/number of samples/other parameters?

         Mesh* fluid_ptr;

         // Discrete Operators
         SparseMatrix<Real> star1;
         SparseMatrix<Real> d0;

   };
}

#endif

