
#include "Fluid.h"
#include "HalfEdge.h"
#include "Quaternion.h"
#include "Vector.h"

namespace DDG
{
   Fluid::Fluid( Mesh*& surface_ptr, const double ProjectionComponent& projectType = DIV )
   :fluid_ptr( &surface )
   {
      // Initiate anything else?
     
      //project pressure to begin with under some criteria:
      if( projectType == CURL ){
         projectCurl( );
      }
      else if( projectType == DIV ){
         projectDivergence( );
      }
      else if ( projectType == HARMONIC ){
         projectHarmonic( );
      }
      else{
         std::cerr << "Projection Component not implemented, exiting" << std::endl;
         return;
      }
   }

   Fluid::~Fluid( void )
   {}

   void Fluid :: flow( const float& dt,
                       const AdvectionScheme& advectType,
                       const ProjectionComponent& projectType,
                       const Interpolation& interpType
                     )
   {
      assert( dt > 0 );

      //TODO: check for interaction/input/Forces

      //advect velocity field
      if( advectType == SEMI_LAGRANGIAN )
      {
         advectSemiLagrangian( dt );
      }
      else{
         std::cerr << "Advection Scheme not implemeted, exiting" << std::endl;
         return;
      }
      updateEdgeWeights( );

      //project pressure under some criteria:
      if( projectType == CURL ){
         projectCurl( );
      }
      else if( projectType == DIV ){
         projectDivergence( );
      }
      else if ( projectType == HARMONIC ){
         projectHarmonic( );
      }
      else{
         std::cerr << "Projection Component not implemented, exiting" << std::endl;
         return;
      }
      updateEdgeWeights( );

      //advectMarker along the flow... update visualization (here? or in view?)

      advectMarkers( dt );
   }

   void Fluid :: advectField( const float& dt )
   {
      for( EdgeIter e = fluid->edges.begin(); e != edges.end(); ++e )
      {
         // Get center of edge
         Vector edge_vec = e->he->flip->vertex - e->he->vertex;
         Vector edge_crossing_coord = e->he->vertex + ( edge_vec / 2 );

         // Integrate the velocity along the edge (get velocity at edge midpoint)
         // by interpolating one of the two bordering triangle faces
         Vector original_velocity = whitneyInterpolate( edge_crossing_coord, e );

         // Determine which halfEdge to use (face associated with velocity direction)
         HalfEdge current_he = cross( edge_vec, original_velocity ) > 0.0 ? e->he : e->he->flip;
 
         // Walk along this direction for distance [ h = v*dt ]
         Quaternion accumulated_angle;
         Vector final_coordinate;

         // Keep track of how much distance is left and the direction you are going (which must be tangent to the surface)
         // This assumes for small timesteps velocity is constant and so we step straight in one direction
         double h_remains = original_velocity.norm() * dt;
         Vector direction = original_velocity.normalize();
         while( h_remains > 0 ) ////what happens at mesh boundaries???
         {
            // determine which of two other edges would be intersected first
            double counter = rayEdgeIntersectDistance( edge_crossing_coord, direction, current_he->next );
            double clockwise = rayEdgeIntersectDistance( edge_crossing_coord, direction, current_he->next->next );
            HalfEdge crossing_he = counter < clockwise ? current_he->next : current_he->next->next;
            double distance = counter < clockwise ? counter : clockwise;

            if( h_remains > distance )
            // Flip to the next face, continue to travel remaider of h
            {
               edge_crossing_coord += distance * direction;
               double theta = acos( dot( crossing_he->flip->face.normal(), current_he->face.normal() ) );
               Quaternion angleChange = Quaternion( cos(theta/2), sin(theta/2)*direction);
               accumulated_angle = accumulated_angle * angleChange;

               direction = ( angleChange * Quaternion( direction ) * angleChange.conj() ).im();
               //TODO: At every crossing should I also change direction to follow new field direction at crossing (?)
                  // If we do this then we are curving along the field path as we travel a distance h, removing the small timestep assumption (!)
                  // Although this is more advanced/accurate, if dt is small enough this will introduce unnecessary computation
                  // This can be done simply by recomputing the velocity (WhitneyInterp) at edge crossings everytime instead of turning/curving the direction vec 

               current_he = crossing_he->flip;
            }
            else{
               final_coordinate = edge_crossing_coord + h_remains * direction;
            }
            h_remains -= distance;
         }

         // Compute interpolated velocity at this point
         Vector acquired_velocity = whitneyInterpolate( final_coordinate, current_he->face );

         // Parallel transport the velocity back to the original face
            // (rotate direction vector back by accumulated angle)
         acquired_velocity = ( accumulated_angle.conj() * Quaternion( acquired_velocity ) * accumulated_angle ).im();

         //TODO:
         assert( acquired_velocity lives on plane of original velocity/face );

         // Update the velocity at initial edge to be the dot product (TODO: why?) of aquired velocity with original velocity
         e->mod_coef = dot( acquired_velocity, edge_vec ); //TODO: should things be normalized here? (original is)
      }
   }

   void Fluid :: projectCurl( void )
   {
      std::cerr << "Projection of Curl is not yet implemented, exiting" << std::endl;
      return;
   }

   void Fluid :: projectDivergence( void )
   // first solve [ d * dp = d * u ] then solve [ u_new = u - dp ] 
   {

      // retrive vector of edge weights ( omega )

      // ROHANS MAGIC CODE

      // Iterate and update
         // Edge.mod_coeff = Edge.ref_coeff - d0Alpha

   }          

   void Fluid :: projectHarmonic( void )
   {
      std::cerr << "Projection of Harmonic is not yet implemented, exiting" << std::endl;
      return;
   }   

   double Fluid :: rayLineIntersectDistance( Vector& coordinate, Vector& direction, HalfEdge& edge )
   // ray should never negatively intersect with halfEdge
   // since rays are emitted internally to a triangle from the border
   {

//      assert that the direction vector is on the plane of the triangle (tangent to the plane)

      Vector edge_vec = edge->flip->vertex - edge->vertex;
      Vector edge_center = edge->vertex + ( edge_vec / 2 );
      Vector plane_normal = cross( edge_vec, edge->face.normal() );
      plane_normal.normalize();
      direction.normalize();

      double denom = dot( plane_normal, direction );
      if( denom > 1e-6 ){
      // if really small or negative I should ignore/skip?
         Vector coord_to_plane_vec = edge_center - coordinate;
         return dot( coord_to_plane_vec, plane_normal );
      }
      std::cerr << "Cause for concern in rayEdgeIntersectDistance?" << std::endl;
      return -1;
   }

   void Fluid :: updateEdgeWeights( )
   {
      for( EdgeIter e = fluid_ptr->edges.begin(); e != edges.end(); ++e )
      {
         // Do we actually care about updating mod_coeff?
         // // double temp = e->ref_coef;
         // e->ref_coef = e->mod_coef;
         // // e->mod_coef = temp;

         e->updateRefCoef();
      }
   }

   Vector whitneyInterpolate( const Vector& coordinate, const Edge& edge )
   {
      assert( coordinate lies on plane created by the edges neighboring triangle faces );

      // Calls whiteneyInterpolate( coordinate, HALF_EDGE ) twice for neighboring triangles of edge

      assert( edge_wise componeny of vectors are the same )

      // averages the flux/normal components to the edge and keeps the same edgewise component
         // TODO: if we find we get a lot of viscosity/damping on our flow, we might find taking the MAX or Min side here is better (?)

      // need to take the averaged component and project it down onto the face that dominates (otherwise the edge_normal component is off the surface)

      //returns this hybrid vector

   }

   void BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k)
   {
       Vector v0 = v_j - v_i;
       Vector v1 = v_k - v_i;
       Vector v2 = coordinate - v_i;

       float d00 = dot( v0, v0 );
       float d01 = dot( v0, v1 );
       float d11 = dot( v1, v1 );
       float d20 = dot( v2, v0 );
       float d21 = dot( v2, v1 );
       float denom = d00 * d11 - d01 * d01;
       a_j = (d11 * d20 - d01 * d21) / denom;
       a_k = (d00 * d21 - d01 * d20) / denom;
       a_i = 1.0f - a_j - a_k;

       assert( 0. <= a_i && a_i <= 1. );
       assert( 0. <= a_j && a_j <= 1. );
       assert( 0. <= a_k && a_k <= 1. );
   }

   /*
      For Details, see Design of Tangent Vector Fields [Fisher et alia 07]

             k
            / \
      C_ki /   \ C_jk
          /     \
        i ------- j 
            C_ij

      Barycentric Coords:
         a_i, a_j, a_k

      Triangle Fact Area:
         Area = face.area()

      Perpendicular Edge Vectors:
         P_ij, P_jk, P_ki
         90 degrees rotation about plane of triangle:
            theta = PI/2;
            Q90 = Quaternion( cos(theta/2), sin(theta/2)* Face_normal );
            P_ij =  ( Q90 ( Pj - Pi ) Q90.conj() ).im()

      Interpolated Vector:
      u( a_i, a_j, a_k ) = ( 
                             ( C_ki*a_k - C_ij*a_j ) * P_jk +
                             ( C_ij*a_i - C_jk*a_k ) * P_ki +
                             ( C_jk*a_j - C_ki*a_i ) * P_ij
                           ) / ( 2 * Area );
   */
   Vector whitneyInterpolate( const Vector& coordinate, const HalfEdge& he )
   {
      assert( coordinate lies on plane created by face of halfEdge );

      // compute whitney interpolation for point
         // Find barycentric weights of coordinate
         // [ rotate 90 deg in plane of triangle ]

         // divide by [ 2 * Area_of_Triangle ]

      assert( vector lies on the plane of the triangle face );

      // return vector
   }

}

