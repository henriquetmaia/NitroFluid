
#include "Fluid.h"
#include "HalfEdge.h"
#include "Quaternion.h"
#include "Vector.h"
#include <Math.h> /* pi */
#include <algorithm>
#include <stdlib.h>

#define EPSILON 1e-6

using namespace std;
namespace DDG
{
  Fluid :: Fluid( Mesh& fluid_ptr, const ProjectionComponent& projectType )
  {
    prescribeVelocityField( fluid_ptr, 0 );
    prescribeDensity( fluid_ptr, 3 );
    buildOperators( fluid_ptr );  

    // Project pressure to guarantee initially divergence free field:
    if( projectType == CURL ){
      projectCurl( fluid_ptr );
    }
    else if( projectType == DIV ){
      projectDivergence( fluid_ptr );
    }
    else if ( projectType == HARMONIC ){
      projectHarmonic( fluid_ptr );
    }
    else{
      cerr << "Projection Component not implemented, exiting" << endl;
      return;
    }
    updateEdgeWeights( fluid_ptr );

  }

  Fluid :: ~Fluid( void )
  {}
   
  void Fluid :: prescribeVelocityField( Mesh& fluid_ptr, int vf )
  {	
    switch ( vf )
    {
      case 0:{
        for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e )
        {
          double constant = 1.0;
          e->setCoef( constant );
          e->updateRefCoef();
        }		
        break; 
      }
      case 1:{
        Vector x( 1, 0, 0 );
        for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e )
        {
          double weight = dot( e->vector(), x );
          e->setCoef( weight );
          e->updateRefCoef();
        }
      	break;
      }
      case 2:{
        for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e )
        {
          double random = (double) rand() / RAND_MAX ;
          e->setCoef( random );
          e->updateRefCoef();
        }   
        break; 
      }
    }    
  }

  void Fluid:: prescribeDensity( Mesh& fluid_ptr, int d )
  {
    switch ( d )
    { // White: [1,1,1] Black: [0,0,0]
      case 0:
      	for( VertexIter v = fluid_ptr.vertices.begin(); v != fluid_ptr.vertices.end(); ++v )
        {
          v->color = Vector( (double) rand() / RAND_MAX, (double) rand() / RAND_MAX, (double) rand() / RAND_MAX );
        }
      	break;
      case 1:
        for( VertexIter v = fluid_ptr.vertices.begin(); v != fluid_ptr.vertices.end(); ++v )
        {
          v->color = (v->position + Vector(1.0,1.0,1.0)) / 2.0;
        }
      	break;
      case 2:
        for( VertexIter v = fluid_ptr.vertices.begin(); v != fluid_ptr.vertices.end(); ++v )
        {
          if( (int) ( v->position.x * 30 ) % 2 == 0 ){
            v->color = Vector( 0.549, 0.839, 1. );// Vector( 1., 1., 1. ); // ( 0.549, 0.839, 1. );         
          }
          else{
            v->color = Vector( 1., .4, 0 );// Vector( 0., 0., 0. ); // ( 1., .4, 0 ); 
          }
        }  
        break;
      case 3:
        for( VertexIter v = fluid_ptr.vertices.begin(); v != fluid_ptr.vertices.end(); ++v )
        {
          if( v->position.y > 0.7 ){
            v->color = Vector( 1.0, 0.0, 0.0 );
          }
          else{
            v->color = Vector( 1.0, 1.0, 1.0 );           
          }
        }
        break;
        // default:
        // break;
    }
  }

  void Fluid :: buildOperators( Mesh& fluid_ptr )
  {
    HodgeStar1Form<Real>::build(fluid_ptr, star1);
    ExteriorDerivative0Form<Real>::build(fluid_ptr, d0);
    ExteriorDerivative1Form<Real>::build(fluid_ptr, d1);
  } 

  void Fluid :: advectVelocitySemiLagrangian( Mesh& fluid_ptr, const float& dt )
  {
    for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e )
    {
      // Get center of edge
      Vector edge_midpoint = e->he->vertex->position + ( e->vector() / 2 );

      // Integrate the velocity along the edge (get velocity at edge midpoint)
      // by interpolating one of the two bordering triangle faces
      HalfEdgeIter starting_he;
      Vector original_velocity = whitneyInterpolateVelocity( edge_midpoint, e, starting_he );
      HalfEdgeIter current_he = starting_he;
      Vector final_coord = edge_midpoint;
      Quaternion accumulated_angle(1);
      assert( insideTriangle( final_coord, current_he ) );
      backtrackAlongField( dt, original_velocity, current_he, final_coord, accumulated_angle );
      assert( insideTriangle( final_coord, current_he ) );

      // Compute interpolated velocity at this point
      Vector acquired_velocity = whitneyInterpolateVelocity( final_coord, current_he ); // current_he is now inside triangle face we care about

      // Parallel transport the velocity back to the original face (rotate direction vector back by accumulated angle)

      acquired_velocity = ( accumulated_angle.conj() * Quaternion( acquired_velocity ) * accumulated_angle ).im();

cout << "accumulated_angle: " << accumulated_angle <<  " dot " << dot( acquired_velocity / acquired_velocity.norm(), starting_he->face->normal() ) << endl;
// cout << "starting_he-face " << starting_he->face->getID() << " current_he-face " << current_he->face->getID() << endl;

      assert( abs( dot( acquired_velocity / acquired_velocity.norm(), starting_he->face->normal() )  ) <= 0.2/*EPSILON*/ ); // acquired_velocity lives on plane of original face 

      // Update the velocity at initial edge to be the dot product of aquired velocity with original velocity
      double advectedWeight = dot( acquired_velocity, e->vector() );
      e->setCoef( advectedWeight );
    }
  }

  void Fluid :: advectColorAlongField( Mesh& fluid_ptr, const float& dt )
  {
    // Update all the face colors based on previous vertex colors
    for( FaceIter f = fluid_ptr.faces.begin(); f != fluid_ptr.faces.end(); ++f )
    {
cout << "face: " << f->getID() << endl;
      // Get center of edge
      Vector face_midpoint = ( f->he->vertex->position + f->he->next->vertex->position + f->he->next->next->vertex->position ) / 3;

      Vector original_velocity = whitneyInterpolateVelocity( face_midpoint, f->he ); // returns a velocity

      HalfEdgeIter final_he = f->he;
      Vector final_coord = face_midpoint;
      Quaternion acc_angl(1);

      assert( insideTriangle( final_coord, final_he ) );
      backtrackAlongField( dt, original_velocity, final_he, final_coord, acc_angl );
      assert( insideTriangle( final_coord, final_he ) && "On way out" );

      Vector acquired_color = barycentricInterpolateColors( final_coord, final_he );

      // Update the velocity at initial edge to be the dot product of aquired velocity with original velocity
      f->color = acquired_color;
    }

    // Update vertex colors by weighting the face colors on the ring around it and area/angle of contributing faces
    for( VertexIter v = fluid_ptr.vertices.begin(); v != fluid_ptr.vertices.end(); ++v )
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

  void Fluid :: backtrackAlongField( const float& dt, const Vector& start_vel, HalfEdgeIter& current_he, Vector& final_coord, Quaternion& accumulated_angle )
  {
      // Walk along this direction for distance [ h = v*dt ]
      assert( accumulated_angle.re() == 1. && accumulated_angle.im().norm() <= EPSILON );

      // Keep track of how much distance is left and the direction you are going (which must be tangent to the surface)
      // This assumes for small timesteps velocity is constant and so we step straight in one direction
      double h_remains = start_vel.norm() * dt;
      Vector direction = -start_vel; direction.normalize();
      assert( abs( dot( current_he->face->normal(), direction ) ) <= EPSILON ); // vector lies on the plane of the triangle face

      unsigned count = 0;
      HalfEdgeIter crossing_he;
      while( h_remains > 0 ) ////what happens at mesh boundaries???
      {
          // determine which of two other edges would be intersected first
cout << endl <<  "count: " << count <<endl;
          count++;
          double distance = intersectRay( final_coord, direction, current_he, h_remains, crossing_he );

          if( h_remains > distance )
          // Flip to the next face, continue to travel remaider of h
          {
            double theta = acos( dot( crossing_he->flip->face->normal(), current_he->face->normal() ) );
            Vector rot_axis = cross( current_he->face->normal(), crossing_he->flip->face->normal() );
            rot_axis.normalize();

            Quaternion angleChange = Quaternion( cos(theta/2), sin(theta/2) * rot_axis );
// cout << "orig_accum_angle: " << accumulated_angle << " theta: " << theta << " angleChange: " << angleChange << endl;
            accumulated_angle = accumulated_angle * angleChange;
// cout << "accumulated_angle_mult: " << accumulated_angle << endl;

            direction = ( angleChange * Quaternion( direction ) * angleChange.conj() ).im();
            direction.normalize();
// cout << "new_dir: " << direction << endl;
            //TODO: At every crossing should I also change direction to follow new field direction at crossing (?)
              // If we do this then we are curving along the field path as we travel a distance h, removing the small timestep assumption (!)
              // Although this is more advanced/accurate, if dt is small enough this will introduce unnecessary computation
              // This can be done simply by recomputing the velocity (WhitneyInterp) at edge crossings everytime instead of turning/curving the direction vec 

            current_he = crossing_he->flip;
            assert( abs( dot( current_he->face->normal(), direction ) ) <= EPSILON ); // vector lies on the plane of the triangle face
            assert( insideTriangle( final_coord, current_he ) );

          }

          h_remains -= distance;
      }
  }

  void Fluid :: projectCurl( Mesh& fluid_ptr )
  {
    cerr << "Projection of Curl is not yet implemented, exiting" << endl;
    return;
  }

  void Fluid :: projectDivergence( Mesh& fluid_ptr )
  // first solve [ d * dp = d * u ] then solve [ u_new = u - dp ] 
  {
    unsigned V = fluid_ptr.vertices.size();
    unsigned E = fluid_ptr.edges.size();

    // retrive vector of edge weights ( u )      
    DenseMatrix<Real> u = DenseMatrix<Real>(E, 1);
    for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e ){

      assert( 0 <= e->getID() && e->getID() < E );
      u( e->getID() ) = fluid_ptr.edges[ e->getID() ].getCoef();
    }

    // SparseMatrix<Real> A( V, V );
    // A = ( d0.transpose() ) * star1 * d0;
    // A.shift( 1e-15 );

    SparseFactor<Real> L;
    fluid_ptr.buildLaplacian();
    fluid_ptr.laplacian.shift( 1e-15 );
    L.build( fluid_ptr.laplacian );

    DenseMatrix<Real> divU( V, 1 );
    divU = ( d0.transpose() ) * star1 * u;

    DenseMatrix<Real> p( V, 1 );
    backsolvePositiveDefinite( L, p, divU );

    DenseMatrix<Real> gradP( E, 1 );
    gradP = d0 * p;

    for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e ){

      double u_new = e->getCoef() - gradP( e->getID() );

// if( e->getID() == 0 ){
//   cout << "u_new: " << u_new << endl;
//   cout << "u_before: " << e->getCoef() << endl;
//   cout << "gradP: " << gradP( e->getID() )<< endl;
// }
      
      e->setCoef( u_new );
    }

  }          

  void Fluid :: projectHarmonic( Mesh& fluid_ptr )
  {
    cerr << "Projection of Harmonic is not yet implemented, exiting" << endl;
    return;
  }   

  double Fluid :: intersectRay( Vector& coordinate, const Vector& direction, const HalfEdgeIter& half_edge, const double tmax, HalfEdgeIter& crossing_half_edge )
  // ray should never negatively intersect with halfEdge
  // since rays are emitted internally to a triangle from the border
  {

    // R = orig + dir * t,  [0 <= t]
    // L = v0 + (vf - v0) * s, [0 <= s <= 1]
/*
     o
    | \
    |  \
  orig  \
    |    \
    |     \
    p ---- q
*/

    assert( insideTriangle( coordinate, half_edge ) );
    assert( abs( dot( half_edge->face->normal(), direction ) ) <= EPSILON ); // vector lies on the plane of the triangle face

    Vector e1 = half_edge->next->vector();
    Vector v1 = half_edge->next->vertex->position - coordinate;
    double t1 = cross( v1, e1 ).norm() / cross( direction, e1 ).norm(); 

    Vector e2 = half_edge->next->next->vector();
    Vector v2 = half_edge->next->next->vertex->position - coordinate;
    double t2 = cross( v2, e2 ).norm() / cross( direction, e2 ).norm();


    double t = 0.0;
    double s = 0.0;

    if ( t1 <= t2 ) {
      t = t1;
      crossing_half_edge = half_edge->next;
    } else if( t2 < t1 ){
    	t = t2;
    	crossing_half_edge = half_edge->next->next;
    }
    else{
      cerr << "_strange_" << " t1 "<< t1 << " t2 " <<t2 << endl;
      exit(EXIT_FAILURE);
    }

    s = cross( coordinate - crossing_half_edge->vertex->position, direction ).norm() / cross( crossing_half_edge->vector(), direction ).norm();


    t = min( t, tmax );  
// cout << "t: " << t << " tmax: " << tmax << " coordinate: " << coordinate << endl;
    assert( insideTriangle( coordinate, crossing_half_edge ) );
    assert( abs( dot( crossing_half_edge->face->normal(), direction ) ) <= EPSILON ); // vector lies on the plane of the triangle face

    if( !insideTriangle( coordinate + direction * t, crossing_half_edge ) ){
      coordinate = coordinate - direction * t;
      assert( insideTriangle( coordinate, crossing_half_edge ) );
    }
    else{
      coordinate = coordinate + direction * t;
      assert( insideTriangle( coordinate, crossing_half_edge ) );
    }
    if( t != tmax && !insideTriangle( coordinate, crossing_half_edge->flip ) ){
      
      // s = ( coordinate - crossing_half_edge->vertex->position ).norm() / crossing_half_edge->vector().norm();
      // s = dot( coordinate - crossing_half_edge->vertex->position , crossing_half_edge->vector() ) /crossing_half_edge->vector().norm() ;

      // s = ((coordinate - crossing_half_edge->vertex->position ) + direction * t).norm() / crossing_half_edge->vector().norm();
cout << "s: " << s << " t1 "<< t1 << " t2 " << t2 <<  endl;
// cout << "direction: " << direction << " coordinate: " << coordinate << " crossing_half_edgevec: " << crossing_half_edge->vector() << " crossing_half_edge pos " << crossing_half_edge->vertex->position <<  endl;
// cout << "v2: " << crossing_half_edge->next->vertex->position << " e2: " << crossing_half_edge->next->vector() << endl;
// cout << "v3: " << crossing_half_edge->next->next->vertex->position << " e3: " << crossing_half_edge->next->next->vector() << endl;
      assert( 0. <= s && s <= 1. );
      coordinate = crossing_half_edge->vertex->position + s * crossing_half_edge->vector();  
          assert( insideTriangle( coordinate, crossing_half_edge->flip ) );    
    }
// cout << "new_Coord: " << coordinate << " direction: " << direction << endl;
    assert( insideTriangle( coordinate, crossing_half_edge ) );

    return t;
  }

  void Fluid :: updateEdgeWeights( Mesh& fluid_ptr )
  {
    for( EdgeIter e = fluid_ptr.edges.begin(); e != fluid_ptr.edges.end(); ++e )
    {
      e->updateRefCoef();
    }
  }

  Vector Fluid :: whitneyInterpolateVelocity( const Vector& coordinate, const EdgeIter& edge, HalfEdgeIter& chosen_he )
  {
    // assert( coordinate lies on plane created by the edges neighboring triangle faces );

    assert( insideTriangle( coordinate, edge->he ) );

    // Calls whitneyInterpolate( coordinate, HALF_EDGE ) twice for neighboring triangles of edge
    Vector left = whitneyInterpolateVelocity( coordinate, edge->he );
    Vector right = whitneyInterpolateVelocity( coordinate, edge->he->flip );

    // check if edge_wise components of vectors are the same
// cout << "dot with edge " << dot( left, edge->vector() ) << " vs " << dot( right, edge->vector() ) << endl;
    assert( abs( dot( left, edge->vector() ) - dot( right, edge->vector() ) )  <= EPSILON );

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

    assert( insideTriangle( coordinate, chosen_he ) );

    return interp_vec;
  }

  Vector Fluid :: barycentricInterpolateColors( const Vector& coordinate, const HalfEdgeIter& he )
  {
    float a_i, a_j, a_k;
    Vector i = he->vertex->position;
    Vector j = he->next->vertex->position;
    Vector k = he->next->next->vertex->position;

    BarycentricWeights( coordinate, i, j, k, a_i, a_j, a_k );

    Vector color;
    color += a_i * he->vertex->color;
    color += a_j * he->next->vertex->color;
    color += a_k * he->next->next->vertex->color;

// Debug
    double a = he->vertex->color.x;
    double b = he->next->vertex->color.x;
    double c = he->next->next->vertex->color.x;
    double max_c = max( max( a, b ), c );
    double min_c = min( min( a, b ), c );

    if( color.x - min_c  < -EPSILON || max_c - color.x < -EPSILON ){
      cout << "a_i: " << a_i << " a_j " << a_j << " a_k " << a_k << " sum: " << a_i + a_j + a_k << endl;
      cerr << "Failed convex color assertion: " << color << " min_c: " << min_c << " max_c: " << max_c  << endl;
      exit(EXIT_FAILURE);
    }
// Remove later

    return color;
  }

  bool Fluid :: insideTriangle( const Vector coordinate, const HalfEdgeIter he )
  {
    Vector v_i = he->vertex->position;
    Vector v_j = he->next->vertex->position;
    Vector v_k = he->next->next->vertex->position;

    double area = cross( v_i - v_j, v_i - v_k ).norm();
    double a_i = cross( coordinate - v_j, coordinate - v_k ).norm() / area;
    double a_j = cross( coordinate - v_k, coordinate - v_i ).norm() / area;
    double a_k = cross( coordinate - v_i, coordinate - v_j ).norm() / area;

    if( 0. <= a_i && a_i <= 1. && 
        0. <= a_j && a_j <= 1. &&
        0. <= a_k && a_k <= 1. &&
        a_i + a_j + a_k - 1. <= EPSILON )
    {
      return true;
    }
    // cerr << "Failing Inside Triangle:: a_i: " << a_i << " a_j: " <<a_j << " a_k: " << a_k << " sum: " << a_i + a_j + a_k  << endl;
    return false;
  }


  void Fluid :: BarycentricWeights( const Vector coordinate, const Vector v_i, const Vector v_j, const Vector v_k, float &a_i, float &a_j, float &a_k)
  {

    double area = cross( v_i - v_j, v_i - v_k ).norm();
    a_i = cross( coordinate - v_j, coordinate - v_k ).norm() / area;
    a_j = cross( coordinate - v_i, coordinate - v_k ).norm() / area;
    a_k = cross( coordinate - v_j, coordinate - v_i ).norm() / area;

// cout << "a_i: " << a_i << endl;
// cout << "cross: "<< cross( coordinate - v_j, coordinate - v_k ).norm()  <<endl;
// cout << "area:" << area << endl;
// cerr << (0. <= a_i &&  a_i <= 1.) << endl;

    assert( 0. <= a_i && a_i <= 1. );
    assert( 0. <= a_j && a_j <= 1. );
    assert( 0. <= a_k && a_k <= 1. );
    assert(  a_i + a_j + a_k - 1. <= EPSILON );

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

    Vector P_ij = cross( he->face->normal(), j - i );
    Vector P_jk = cross( he->face->normal(), k - j );
    Vector P_ki = cross( he->face->normal(), i - k );

    assert( abs( dot( P_ij, ( j - i ) ) ) <= EPSILON );
    assert( abs( dot( P_jk, ( k - j ) ) ) <= EPSILON );
    assert( abs( dot( P_ki, ( i - k ) ) ) <= EPSILON );

    //Retrieve edge weights
    double C_ij = he->weight();
    double C_jk = he->next->weight();
    double C_ki = he->next->next->weight();

    // Return Interpolated Vector 
    Vector interp_vec = ( 
                          ( C_ki * a_k  -  C_ij * a_j ) * P_jk +
                          ( C_ij * a_i  -  C_jk * a_k ) * P_ki +
                          ( C_jk * a_j  -  C_ki * a_i ) * P_ij
                        ) / ( 2 * he->face->area() );

    assert( abs( dot( he->face->normal(), interp_vec ) ) <= EPSILON ); // vector lies on the plane of the triangle face

    return interp_vec;
  }

}

