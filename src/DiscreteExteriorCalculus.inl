#include "DiscreteExteriorCalculus.h"

namespace DDG
{
   template <class T>
   void HodgeStar0Form<T> :: build(const Mesh& mesh, SparseMatrix<T>& star0) {
	
	int V = mesh.vertices.size();
	star0.resize(V,V);
	for(vector<Vertex>::const_iterator it = mesh.vertices.begin(); it != mesh.vertices.end(); it++) {
		star0(it->getID(), it->getID()) = it->dualArea(); 
	}
   }

   template <class T>
   void HodgeStar1Form<T> :: build(const Mesh& mesh, SparseMatrix<T>& star1) {
	
	int E = mesh.edges.size();
	star1.resize(E,E);
	for(vector<Edge>::const_iterator it = mesh.edges.begin(); it != mesh.edges.end(); it++) {
		HalfEdgeIter he = it->he;
		double cotAlpha = he->cotan();
		
		he = he->flip;
		double cotBeta = he->cotan();
				
		star1(it->getID(), it->getID()) = 0.5 * (cotAlpha + cotBeta); 
	}	
   }

   template <class T>
   void HodgeStar2Form<T> :: build(const Mesh& mesh, SparseMatrix<T>& star2) {
	
	int F = mesh.faces.size();
	star2.resize(F,F);
	for(vector<Face>::const_iterator it = mesh.faces.begin(); it != mesh.faces.end(); it++) {
		star2(it->getID(), it->getID()) = 1.0 / it->area(); 
	}
   }

   template< class T >
   void ExteriorDerivative0Form<T> :: build(const Mesh& mesh, SparseMatrix<T>& d0) {
	
	int E = mesh.edges.size();
	int V = mesh.vertices.size();
	d0.resize(E,V);	
   	
	for(vector<Edge>::const_iterator it = mesh.edges.begin(); it != mesh.edges.end(); it++) {
		d0(it->getID(), it->he->vertex->getID()) = -1;
		d0(it->getID(), it->he->flip->vertex->getID()) = 1;	
	}
   }

   template< class T >
   void ExteriorDerivative1Form<T> :: build(const Mesh& mesh, SparseMatrix<T>& d1) {
	
	int F = mesh.faces.size();
	int E = mesh.edges.size();
 	d1.resize(F,E);
	
   	for(vector<Face>::const_iterator it = mesh.faces.begin(); it != mesh.faces.end(); it++) {
		
		HalfEdgeIter he = it->he;
				 
		do {	
			int val = 1;
			if(he->flip == he->edge->he) val = -1;
			d1(it->getID(), he->edge->getID()) = val;
	
			he = he->next;			
		} while(he != it->he);
	}
    }
}
