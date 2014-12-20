#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

#include "Viewer.h"
#include "Image.h"
#include "Fluid.h"

namespace DDG
{
   // declare static member variables
   Mesh Viewer::mesh;
   Viewer::RenderMode Viewer::mode = renderShaded;
   GLuint Viewer::surfaceDL = 0;
   int Viewer::windowSize[2] = { 512, 512 };
   Camera Viewer::camera;
   Shader Viewer::shader;
   Fluid* fluid;
   bool color_markers = false;

   void Viewer :: init( void )
   {
      // using Fluid::ProjectionComponent;
      // ProjectionComponent d = 
    for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e )
    {
      if( e->getID() == 4 ){
        std::cout << "before init " << e->getCoef() << std::endl;
      }
    }

      fluid = new Fluid( mesh, Fluid::DIV );
      restoreViewerState();
      initGLUT();
      initGL();
      initGLSL();
    for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e )
    {
      if( e->getID() == 4 ){
        std::cout << "after init: " << e->getCoef() << std::endl;
      }
    }
      updateDisplayList();
   
      glutMainLoop();
   }

   void Viewer :: initGL( void )
   {
      glClearColor( .5, .5, .5, 1. );
   }
   
   void Viewer :: initGLUT( void )
   {
      int argc = 0;
      vector< vector<char> > argv(1);
   
      // initialize window
      glutInitWindowSize( windowSize[0], windowSize[1] );
      glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
      glutInit( &argc, (char**)&argv );
      glutCreateWindow( "DDG" );
   
      // specify callbacks
      glutDisplayFunc  ( Viewer::display  );
      glutIdleFunc     ( Viewer::idle     );
      glutKeyboardFunc ( Viewer::keyboard );
      glutSpecialFunc  ( Viewer::special  );
      glutMouseFunc    ( Viewer::mouse    );
      glutMotionFunc   ( Viewer::motion   );
   
      // initialize menus
      int viewMenu = glutCreateMenu( Viewer::view );
      glutSetMenu( viewMenu );
      glutAddMenuEntry( "[s] Smooth Shaded",  menuSmoothShaded );
      glutAddMenuEntry( "[f] Wireframe",      menuWireframe    );
      glutAddMenuEntry( "[v] Vector Field",   menuVectorField  );
      glutAddMenuEntry( "[↑] Zoom In",    menuZoomIn       );
      glutAddMenuEntry( "[↓] Zoom Out",   menuZoomOut      );

      int mainMenu = glutCreateMenu( Viewer::menu );
      glutSetMenu( mainMenu );
      glutAddMenuEntry( "[space] Process Mesh", menuProcess    );
      glutAddMenuEntry( "[r] Reset Mesh",       menuResetMesh  );
      glutAddMenuEntry( "[w] Write Mesh",       menuWriteMesh  );
      glutAddMenuEntry( "[\\] Screenshot",      menuScreenshot );
      glutAddMenuEntry( "[esc] Exit",           menuExit       );
      glutAddSubMenu( "View", viewMenu );
      glutAttachMenu( GLUT_RIGHT_BUTTON );
   }

   void Viewer :: initGLSL( void )
   {
      shader.loadVertex( "shaders/vertex.glsl" );
      shader.loadFragment( "shaders/fragment.glsl" );
   }
   
   void Viewer :: menu( int value )
   {
      switch( value )
      {
         case( menuProcess ):
            mProcess();
            break;
         case( menuResetMesh ):
            mResetMesh();
            break;
         case( menuWriteMesh ):
            mWriteMesh();
            break;
         case( menuScreenshot ):
            mScreenshot();
            break;
         case( menuExit ):
            mExit();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: view( int value )
   {
      switch( value )
      {
         case( menuSmoothShaded ):
            mSmoothShaded();
            break;
         case( menuWireframe ):
            mWireframe();
            break;
         case( menuVectorField ):
            mVectorField();
            break;
         case( menuZoomIn ):
            mZoomIn();
            break;
         case( menuZoomOut ):
            mZoomOut();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: mProcess( void )
   {
      // TODO process geometry here!
      static double sim_time = 0.;
      if( fluid == NULL ){
         std::cerr << "[Viewer::mProcess] Fluid not initialized" << std::endl;
         exit( EXIT_FAILURE );
      }

      const float dt = 0.01;
      // const Fluid::AdvectionScheme advectType = Fluid::SEMI_LAGRANGIAN;
      // const Fluid::ProjectionComponent projectType = Fluid::DIV;

for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); ++e ){
if( e->getID() == 4 ){
std::cout << "Before advect: " << e->getCoef() << std::endl;
}
}

      {
         assert( dt > 0. );
         //TODO: check for interaction/input/Forces
/*
         //advect velocity field
         if( advectType == Fluid::SEMI_LAGRANGIAN )
         {
            fluid->advectVelocitySemiLagrangian( mesh, dt );
         }
         else{
            std::cerr << "Advection Scheme not implemeted, exiting" << std::endl;
            return;
         }
         std::cout << "AFTER ADVECT: " << std::endl;
         fluid->updateEdgeWeights( mesh );


         //project pressure under some criteria:
         if( projectType == Fluid::CURL ){
            fluid->projectCurl( mesh );
         }
         else if( projectType == Fluid::DIV ){
            fluid->projectDivergence( mesh );
         }
         else if ( projectType == Fluid::HARMONIC ){
            fluid->projectHarmonic( mesh );
         }
         else{
            std::cerr << "Projection Component not implemented, exiting" << std::endl;
            return;
         }
         std::cout << "AFTER PROJECT: " << std::endl;
         fluid->updateEdgeWeights( mesh );
*/
         fluid->advectColorAlongField( mesh, dt );
      }

      updateDisplayList();
      sim_time += dt;
      std::cout << "t: " << sim_time << std::endl;

   }
   
   void Viewer :: mResetMesh( void )
   {
      mesh.reload();
      updateDisplayList();
   }
   
   void Viewer :: mWriteMesh( void )
   {
      mesh.write( "out.obj" );
   }
   
   void Viewer :: mExit( void )
   {
      // delete fluid;
      storeViewerState();
      exit( 0 );
   }
   
   void Viewer :: mSmoothShaded( void )
   {
      mode = renderShaded;
      updateDisplayList();
   }

   void Viewer :: mVectorField( void )
   {
      mode = renderVectorField;
      updateDisplayList();
   }
   
   void Viewer :: mWireframe( void )
   {
      mode = renderWireframe;
      updateDisplayList();
   }

   void Viewer :: mZoomIn( void )
   {
      camera.zoomIn();
   }

   void Viewer :: mZoomOut( void )
   {
      camera.zoomOut();
   }

   void Viewer :: mScreenshot( void )
   {
      static int index = 0;
   
      // get window width and height
      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      int w = view[2];
      int h = view[3];
   
      // get pixels
      Image image( w, h );
      glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );
   
      stringstream filename;
      filename << "frames/viewer" << setw(8) << setfill( '0' ) << index << ".tga";
      image.write( filename.str().c_str() );
   
      index++;
   }
   
   void Viewer :: keyboard( unsigned char c, int x, int y )
   {
      switch( c )
      {
         case 's':
            mSmoothShaded();
            break;
         case 'f':
            mWireframe();
            break;
         case 'v':
            mVectorField();
            break;
         case 'w':
            mWriteMesh();
            break;
         case 'c':
            color_markers = !color_markers;
            updateDisplayList();
            break;            
         case 'r':
            mResetMesh();
            break;
         case '\\':
            mScreenshot();
            break;
         case ' ':
            mProcess();
            break;
         case 'd':
// fluid->prescribeVelocityField( 1 ); 
//   fluid->prescribeDensity( 1 );	
            break;
         case 'q':
         case 27:
            mExit();
            break;
         default:
            break;
      }
   }

   void Viewer :: special( int i, int x, int y )
   {
      switch( i )
      {
         case GLUT_KEY_UP:
            camera.zoomIn();
            break;
         case GLUT_KEY_DOWN:
            camera.zoomOut();
            break;
         case 27:
            mExit();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: display( void )
   {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

      shader.enable();
   
      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
   
      GLint viewport[4];
      glGetIntegerv( GL_VIEWPORT, viewport );
      double aspect = (double) viewport[2] / (double) viewport[3];
      const double fovy = 50.;
      const double clipNear = .01;
      const double clipFar = 1000.;
      gluPerspective( fovy, aspect, clipNear, clipFar );
   
      glMatrixMode( GL_MODELVIEW );
      glLoadIdentity();

      Quaternion    eye = Vector( 0., 0., -2.5*camera.zoom );
      Quaternion center = Vector( 0., 0., 0. );
      Quaternion     up = Vector( 0., 1., 0. );

      gluLookAt(    eye[1],    eye[2],    eye[3],
                 center[1], center[2], center[3],
                     up[1],     up[2],     up[3] );


      Quaternion r = camera.currentRotation();

      eye = r.conj() * eye * r;
      GLint uniformEye = glGetUniformLocation( shader, "eye" );
      glUniform3f( uniformEye, eye[1], eye[2], eye[3] );

      Quaternion light = Vector( -1., 1., -2. );
      light = r.conj() * light * r;
      GLint uniformLight = glGetUniformLocation( shader, "light" );
      glUniform3f( uniformLight, light[1], light[2], light[3] );

      camera.setView();
   
      drawSurface();

      shader.disable();
   
      glutSwapBuffers();
   }

   void Viewer :: drawSurface( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_DEPTH_TEST );
      glEnable( GL_LIGHTING );
   
      glCallList( surfaceDL );
   
      glPopAttrib();
   }
   
   void Viewer :: drawMesh( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_POLYGON_OFFSET_FILL );
      glPolygonOffset( 1., 1. );
   
      drawPolygons();
   
      glDisable( GL_POLYGON_OFFSET_FILL );
   
      if( mode == renderWireframe )
      {
         shader.disable();
         drawWireframe();
      }
      else if( mode == renderVectorField )
      {
         shader.disable();
         drawVectorField();
      }

      drawIsolatedVertices();

      glPopAttrib();
   }

   void Viewer :: drawPolygons( void )
   {
      // Default color
      glColor3f( 0.549, 0.839, 1. );

      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         if( f->isBoundary() ) continue;

         if( color_markers )
         {
            //INTERP_COLOR
            Vector color;
            HalfEdgeCIter h = f->he;
            do
            {
               color += h->vertex->color;
               h = h->next;
            }
            while( h != f->he );
            color = color / 3;
            // color.normalize();
            glColor3f( color.x, color.y, color.z );
            if( f->getID() == 4 ){
               std::cout << "F4c: " << f->color <<std::endl;
               std::cout << "final_color: " << color <<std::endl;            
            }
         }

         glBegin( GL_POLYGON );
         if( mode == renderWireframe || mode == renderVectorField )
         {
            Vector N = f->normal();
            glNormal3dv( &N[0] );
         }

         HalfEdgeCIter he = f->he;
         do
         {
            if( mode != renderWireframe && mode != renderVectorField ) // Comment back in for triangles in vectorfield
            {
               Vector N = he->vertex->normal();
               glNormal3dv( &N[0] );
            }

            glVertex3dv( &he->vertex->position[0] );

            he = he->next;
         }
         while( he != f->he );

         glEnd();
      }
   }

   void Viewer :: drawVectorField( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glDisable( GL_LIGHTING );
      glColor4f( 0., 0., 0., .5 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

      glBegin( GL_LINES );
      for( FaceCIter f  = mesh.faces.begin(); f != mesh.faces.end(); ++f )
      {
         Vector face_midpoint = ( f->he->vertex->position + f->he->next->vertex->position + f->he->next->next->vertex->position ) / 3;
         Vector field_vel = fluid->whitneyInterpolateVelocity( face_midpoint, f->he );
         field_vel *= 100 * f->area();

         Vector endpoint = face_midpoint + field_vel;
         glVertex3dv( &face_midpoint[0] );
         glVertex3dv( &endpoint[0] );


      face_midpoint = ( f->he->vertex->position * 0.8 + f->he->next->vertex->position * 0.1 + f->he->next->next->vertex->position * 0.1 );
      field_vel = fluid->whitneyInterpolateVelocity( face_midpoint, f->he );
      field_vel *= 100 * f->area();

      endpoint = face_midpoint + field_vel;
      glVertex3dv( &face_midpoint[0] );
      glVertex3dv( &endpoint[0] );

      face_midpoint = ( f->he->vertex->position * 0.1 + f->he->next->vertex->position * 0.8 + f->he->next->next->vertex->position * 0.1 );
      field_vel = fluid->whitneyInterpolateVelocity( face_midpoint, f->he );
      field_vel *= 100 * f->area();

      endpoint = face_midpoint + field_vel;
      glVertex3dv( &face_midpoint[0] );
      glVertex3dv( &endpoint[0] );

      face_midpoint = ( f->he->vertex->position * 0.1 + f->he->next->vertex->position * 0.1 + f->he->next->next->vertex->position * 0.8 );
      field_vel = fluid->whitneyInterpolateVelocity( face_midpoint, f->he );
      field_vel *= 100 * f->area();

      endpoint = face_midpoint + field_vel;
      glVertex3dv( &face_midpoint[0] );
      glVertex3dv( &endpoint[0] );


      }

      
      glEnd();

      glPopAttrib();
   }


   void Viewer :: drawWireframe( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glDisable( GL_LIGHTING );
      glColor4f( 0., 0., 0., .5 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

      glBegin( GL_LINES );
      for( EdgeCIter e  = mesh.edges.begin();
            e != mesh.edges.end();
            e ++ )
      {
         glVertex3dv( &e->he->vertex->position[0] );
         glVertex3dv( &e->he->flip->vertex->position[0] );
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: drawIsolatedVertices( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      // draw with big, round, red dots
      glPointSize( 5 );
      glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
      glEnable( GL_POINT_SMOOTH );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glColor4f( 1., 0., 0., 1. ); // red

      glBegin( GL_POINTS );
      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         if( v->isIsolated() )
         {
            glVertex3dv( &v->position[0] );
         }
      }
      glEnd();

      glPopAttrib();
   }
   
   void Viewer :: updateDisplayList( void )
   {
      if( surfaceDL )
      {
         glDeleteLists( surfaceDL, 1 );
         surfaceDL = 0;
      }
   
      surfaceDL = glGenLists( 1 );
   
      glNewList( surfaceDL, GL_COMPILE );
      drawMesh();
      glEndList();
   }
   
   void Viewer :: mouse( int button, int state, int x, int y )
   {
      camera.mouse( button, state, x, y );
   }

   void Viewer :: motion( int x, int y )
   {
      camera.motion( x, y );
   }
   
   void Viewer :: idle( void )
   {
      camera.idle();
      glutPostRedisplay();
   }

   void Viewer :: storeViewerState( void )
   {
      ofstream out( ".viewer_state.txt" );

      out << camera.rLast[0] << endl;
      out << camera.rLast[1] << endl;
      out << camera.rLast[2] << endl;
      out << camera.rLast[3] << endl;

      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      out << view[2] << endl;
      out << view[3] << endl;

      out << (int) mode << endl;
   }

   void Viewer :: restoreViewerState( void )
   {
      ifstream in( ".viewer_state.txt" );

      if( !in.is_open() ) return;

      in >> camera.rLast[0];
      in >> camera.rLast[1];
      in >> camera.rLast[2];
      in >> camera.rLast[3];

      in >> windowSize[0];
      in >> windowSize[1];

      int m;
      in >> m;
      mode = (RenderMode) m;
   }
}

