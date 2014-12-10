#ifndef DDG_DISCRETEEXTERIORCALCULUS_H
#define DDG_DISCRETEEXTERIORCALCULUS_H

#include "Mesh.h"
#include "SparseMatrix.h"

namespace DDG
{
   template< class T > struct HodgeStar0Form { static void build( const Mesh& mesh, SparseMatrix<T>& star0 ); };
   template< class T > struct HodgeStar1Form { static void build( const Mesh& mesh, SparseMatrix<T>& star1 ); };
   template< class T > struct HodgeStar2Form { static void build( const Mesh& mesh, SparseMatrix<T>& star2 ); };
   template< class T > struct ExteriorDerivative0Form { static void build( const Mesh& mesh, SparseMatrix<T>& d0 ); };
   template< class T > struct ExteriorDerivative1Form { static void build( const Mesh& mesh, SparseMatrix<T>& d1 ); };
}

#include "DiscreteExteriorCalculus.inl"

#endif
