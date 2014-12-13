
#include "Fluid.h"
#include "HalfEdge.h"
#include "Quaternion.h"
#include "Vector.h"
#include <Math.h> /* pi */
#include <algorithm>
#include <stdlib.h>

#define EPSILON 1e-6

namespace DDG
{
  Fluid :: Fluid( Mesh* surface_ptr, const ProjectionComponent& projectType )
  :fluid_ptr( surface_ptr )
  {
    prescribeVelocityField( 1 );
    prescribeDensity( 1 );
    buildOperators( );  

    // SparseMatrix<Real> d = d1 * d0;
    // std::cout << "d: " << d1 * d0 << std::endl;
    // std::cout << "*: " << star1 <<std::endl;

    // Project pressure to guarantee initially divergence free field:
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
    updateEdgeWeights();
  }

  Fluid :: ~Fluid( void )
  {}
   
  void Fluid :: prescribeVelocityField( int vf )
  {	
    switch ( vf )
    {
      case 0:
        for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e )
        {
          e->setCoef( 8.0 );
        }		
        break; 
      case 1:
        Vector x( 8, 0, 0 );
        for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e )
        {
      	  Vector v = e->he->flip->vertex->position - e->he->vertex->position;
          double weight = dot( v, x );
          e->setCoef( weight );
          e->updateRefCoef();
        }
	break;
    }    
  }

  void Fluid:: prescribeDensity( int d )
  {
    switch ( d )
    {
      case 0:
      	for( VertexIter v = fluid_ptr->vertices.begin(); v != fluid_ptr->vertices.end(); ++v )
        {
          v->color = Vector( (double) rand() / RAND_MAX, (double) rand() / RAND_MAX, (double) rand() / RAND_MAX );
        }
      	break;
      case 1:
        for( VertexIter v = fluid_ptr->vertices.begin(); v != fluid_ptr->vertices.end(); ++v )
        {
          v->color = (v->position + Vector(1.0,1.0,1.0)) / 2.0;
        }
      	// texture map???
      	break;
        // default:
        // break;
    }
  }

  void Fluid :: buildOperators( )
  {
    HodgeStar1Form<Real>::build(*fluid_ptr, star1);
    ExteriorDerivative0Form<Real>::build(*fluid_ptr, d0);
    ExteriorDerivative1Form<Real>::build(*fluid_ptr, d1);
  } 

  void Fluid :: advectVelocitySemiLagrangian( const float& dt )
  {
    for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e )
    {
      // Get center of edge
      Vector edge_vec = e->he->flip->vertex->position - e->he->vertex->position;
      Vector edge_midpoint = e->he->vertex->position + ( edge_vec / 2 );

      // Integrate the velocity along the edge (get velocity at edge midpoint)
      // by interpolating one of the two bordering triangle faces
      HalfEdgeIter current_he;
      Vector original_velocity = whitneyInterpolateVelocity( edge_midpoint, e, current_he );

      Vector final_coord;
      Quaternion accumulated_angle;
      backtrackAlongField( dt, edge_midpoint, original_velocity, current_he, final_coord, accumulated_angle );

      // Compute interpolated velocity at this point
      Vector acquired_velocity = whitneyInterpolateVelocity( final_coord, current_he ); // current_he is now inside triangle face we care about

      // Parallel transport the velocity back to the original face (rotate direction vector back by accumulated angle)
      acquired_velocity = ( accumulated_angle.conj() * Quaternion( acquired_velocity ) * accumulated_angle ).im();

      assert( std::abs( dot( acquired_velocity, starting_he->face->normal() )  ) <= EPSILON ); // acquired_velocity lives on plane of original face 

      // Update the velocity at initial edge to be the dot product of aquired velocity with original velocity
      double advectedWeight = dot( acquired_velocity, edge_vec );


      if( e->getID() == 4 ){
        std::cout << "before: " << e->getCoef() <<" advectedWeight: " << advectedWeight << std::endl;
      }
      e->setCoef( advectedWeight );
    }
  }

  void Fluid :: advectColorAlongField( const float& dt )
  {
    // Update all the face colors based on previous vertex colors
    for( FaceIter f = fluid_ptr->faces.begin(); f != fluid_ptr->faces.end(); ++f )
    {
      // Get center of edge
      Vector face_midpoint = ( f->he->vertex->position + f->he->next->vertex->position + f->he->next->next->vertex->position ) / 3;

      Vector original_velocity = whitneyInterpolateVelocity( face_midpoint, f->he ); // returns a velocity

      HalfEdgeIter final_he = f->he;
      Vector final_coord;
      Quaternion acc_angl;

      backtrackAlongField( dt, face_midpoint, original_velocity, final_he, final_coord, acc_angl );

      Vector acquired_color = barycentricInterpolateColors( final_coord, final_he );

      // Update the velocity at initial edge to be the dot product of aquired velocity with original velocity
      f->color = acquired_color;
    }

    // Update vertex colors by weighting the face colors on the ring around it and area/angle of contributing faces
    for( VertexIter v = fluid_ptr->vertices.begin(); v != fluid_ptr->vertices.end(); ++v )
    {
      HalfEdgeCIter h = v->he;
      Vector accumulated_color;
      do
      {
         accumulated_color +=  h->face->color * ( h->face->area() / 3 );
         h = h->flip->next;
      }
      while( h != v->he );
      v->color = accumulated_color / v->dualArea();
    }
  }

  void Fluid :: backtrackAlongField( const float& dt, const Vector start_pt, const Vector& start_vel, HalfEdgeIter& current_he, Vector& final_coord, Quaternion& accumulated_angle )
  {
      // Walk along this direction for distance [ h = v*dt ]
      final_coord = start_pt;

      // Keep track of how much distance is left and the direction you are going (which must be tangent to the surface)
      // This assumes for small timesteps velocity is constant and so we step straight in one direction
      double h_remains = start_vel.norm() * dt;
      Vector direction = -start_vel; direction.normalize();
      while( h_remains > 0 ) ////what happens at mesh boundaries???
      {
          // determine which of two other edges would be intersected first

          HalfEdgeIter crossing_he;
          double distance = intersectRay( final_coord, direction, current_he, h_remains, crossing_he );

          if( h_remains > distance )
          // Flip to the next face, continue to travel remaider of h
          {
            /* edge_crossing_coord += distance * direction; */

            double theta = acos( dot( crossing_he->flip->face->normal(), current_he->face->normal() ) );
            Quaternion angleChange = Quaternion( cos(theta/2), sin(theta/2)*direction);
            accumulated_angle = accumulated_angle * angleChange;

            direction = ( angleChange * Quaternion( direction ) * angleChange.conj() ).im();
            //TODO: At every crossing should I also change direction to follow new field direction at crossing (?)
              // If we do this then we are curving along the field path as we travel a distance h, removing the small timestep assumption (!)
              // Although this is more advanced/accurate, if dt is small enough this will introduce unnecessary computation
              // This can be done simply by recomputing the velocity (WhitneyInterp) at edge crossings everytime instead of turning/curving the direction vec 

            current_he = crossing_he->flip;
          }

          h_remains -= distance;
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
    unsigned V = fluid_ptr->vertices.size();
    unsigned E = fluid_ptr->edges.size();

    // retrive vector of edge weights ( u )      
    DenseMatrix<Real> u = DenseMatrix<Real>(E, 1);
    for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e ){
      assert( 0 <= e->index && e->index < E );
      u( e->getID() ) = fluid_ptr->edges[ e->getID() ].getCoef();
    }

    SparseMatrix<Real> A( V, V );
    A = ( d0.transpose() ) * star1 * d0;
    A.shift( 1e-15 );

    SparseFactor<Real> L;
    L.build( A );
    DenseMatrix<Real> divU( V, 1 );
    divU = ( d0.transpose() ) * star1 * u;

    DenseMatrix<Real> p( V, 1 );
    backsolvePositiveDefinite( L, p, divU );

    DenseMatrix<Real> gradP( E, 1 );
    gradP = d0 * p;

    for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e ){
      double u_new = e->getCoef() - gradP( e->getID() );
            if( e->getID() == 4 ){
      std::cout << "getCoef before u_new: " << e->getCoef() << std::endl;
            }
      e->setCoef( u_new );
        if( e->getID() == 4 ){
            std::cout << "getCoef after u_new: " << e->getCoef() << std::endl;
    }
    } 

    // std::cout << "A: " << A << std::endl;
    // std::cout << "divU: " << divU << std::endl;
    // std::cout << "gradP: " << gradP << std::endl;

  }          

  void Fluid :: projectHarmonic( void )
  {
    std::cerr << "Projection of Harmonic is not yet implemented, exiting" << std::endl;
    return;
  }   

  double Fluid :: intersectRay( Vector& coordinate, const Vector& direction, const HalfEdgeIter& half_edge, const double tmax, HalfEdgeIter& crossing_half_edge )
  // ray should never negatively intersect with halfEdge
  // since rays are emitted internally to a triangle from the border
  {

    //  !!! cannot assume ray is shot from center of edge or any edge at all,
    // arbitrary ray, that lives on the face of triangle, intersect with edges of triangle
    // use 'half_edge' as handle to this triangle,
    // coordinate and direction give you the ray,
    // return t, update coordinate, and update crossing_he 

    Vector e1 = half_edge->flip->vertex->position - half_edge->vertex->position;
    Vector v1 = half_edge->vertex->position - coordinate;
    double t1 = cross( v1, e1 ).norm() / cross( direction, e1 ).norm(); 

    Vector e2 = half_edge->next->flip->vertex->position - half_edge->next->vertex->position;
    Vector v2 = half_edge->next->vertex->position - coordinate;
    double t2 = cross( v2, e2 ).norm() / cross( direction, e2 ).norm();

    double t = 0.0;
    if ( t1 < t2 ) {
      t = t1;
      crossing_half_edge = half_edge->next;
    } else {
    	t = t2;
    	crossing_half_edge = half_edge->next->next;
    }

    t = std::min( t, tmax );  
 
    coordinate = coordinate + direction * t; 
    return t;
  }

  void Fluid :: updateEdgeWeights( )
  {
    for( EdgeIter e = fluid_ptr->edges.begin(); e != fluid_ptr->edges.end(); ++e )
    {
      e->updateRefCoef();
    }
  }

  Vector Fluid :: whitneyInterpolateVelocity( const Vector& coordinate, const EdgeIter& edge, HalfEdgeIter& chosen_he )
  {
    assert( coordinate lies on plane created by the edges neighboring triangle faces );

    // Calls whitneyInterpolate( coordinate, HALF_EDGE ) twice for neighboring triangles of edge
    Vector left = whitneyInterpolateVelocity( coordinate, edge->he );
    Vector right = whitneyInterpolateVelocity( coordinate, edge->he->flip );

    // check if edge_wise components of vectors are the same
    assert( dot( left, ( e->he->flip->vertex->position - e->he->vertex->position ) ) == dot( right, ( e->he->flip->vertex->position - e->he->vertex->position ) ) );

    // averages the flux/normal components to the edge and keeps the same edgewise component?
    // need to take the averaged component and project it down onto the face that dominates (otherwise the edge_normal component is off the surface)
       // TODO: if we find we get a lot of viscosity/damping on our flow, we might find taking the MAX or Min side here is better (?)
    // ...For now we just take the largest.

    Vector interp_vec;
    if( left.norm() > right.norm() ){
      chosen_he = edge->he;
      interp_vec = left;
    }
    else{
      chosen_he = edge->he->flip;
      interp_vec = right;
    }

    return interp_vec;
  }

  Vector Fluid :: barycentricInterpolateColors( const Vector& coordinate, const HalfEdgeIter& he )
  {
    float a_i, a_j, a_k;
    Vector i = he->vertex->position;
    Vector j = he->next->vertex->position;
    Vector k = he->next->next->vertex->position;
    BarycentricWeights( coordinate, i, j, k, a_i, a_j, a_k ); // Already asserts coordinate lies on triangle face

    Vector color;
    color += a_i * he->vertex->color;
    color += a_j * he->next->vertex->color;
    color += a_k * he->next->next->vertex->color;
    return color;
  }

  void Fluid :: BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k)
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

    assert( 0. <= a_i && a_i <= 1. && "Coordinate lies outside of triangle!" );
    assert( 0. <= a_j && a_j <= 1. && "Coordinate lies outside of triangle!" );
    assert( 0. <= a_k && a_k <= 1. && "Coordinate lies outside of triangle!" );
  }

  /*
    For details, see Design of Tangent Vector Fields [Fisher et alia 07]

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
  Vector Fluid :: whitneyInterpolateVelocity( const Vector& coordinate, const HalfEdgeIter& he )
  { // compute whitney interpolation vector quantity for a point inside a triangle

    // Find barycentric weights of coordinate:
    float a_i, a_j, a_k;
    Vector i = he->vertex->position;
    Vector j = he->next->vertex->position;
    Vector k = he->next->next->vertex->position;
    BarycentricWeights( coordinate, i, j, k, a_i, a_j, a_k ); // Already asserts coordinate lies on triangle face

    // Compute Perpendicular Edge Vectors:
    double theta = M_PI / 2;
    Quaternion Q_perp = Quaternion( cos(theta/2), sin(theta/2) * he->face->normal() );
    Vector P_ij = ( Q_perp * ( j - i ) * Q_perp.conj() ).im();
    Vector P_jk = ( Q_perp * ( k - j ) * Q_perp.conj() ).im();
    Vector P_ki = ( Q_perp * ( i - k ) * Q_perp.conj() ).im();

    assert( std::abs( dot( P_ij, ( j - i ) ) ) <= EPSILON );
    assert( std::abs( dot( P_jk, ( k - j ) ) ) <= EPSILON );
    assert( std::abs( dot( P_ki, ( i - k ) ) ) <= EPSILON );

    //Retrieve edge weights
    double C_ij = he->edge->getCoef();
    double C_jk = he->next->edge->getCoef();
    double C_ki = he->next->next->edge->getCoef();

    // Return Interpolated Vector 
    Vector interp_vec = ( 
                          ( C_ki * a_k  -  C_ij * a_j ) * P_jk +
                          ( C_ij * a_i  -  C_jk * a_k ) * P_ki +
                          ( C_jk * a_j  -  C_ki * a_i ) * P_ij
                        ) / ( 2 * he->face->area() );
    assert( std::abs( dot( he->face->normal(), interp_vec ) ) <= EPSILON ); // vector lies on the plane of the triangle face
    return interp_vec;
  }

}

