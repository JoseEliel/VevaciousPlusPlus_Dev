/*
 * RealMassesSquaredMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/OldRealMassesSquaredMatrix.hpp"

namespace VevaciousPlusPlus
{

  OldRealMassesSquaredMatrix::OldRealMassesSquaredMatrix( size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    OldMassesSquaredFromMatrix< double >( numberOfRows,
                                          attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    PolynomialSum() )
  {
    // This constructor is just an initialization list.
  }

  OldRealMassesSquaredMatrix::OldRealMassesSquaredMatrix(
                                  OldRealMassesSquaredMatrix const& copySource ) :
    OldMassesSquaredFromMatrix< double >( copySource ),
    matrixElements( copySource.matrixElements )
  {
    // This constructor is just an initialization list.
  }

  OldRealMassesSquaredMatrix::OldRealMassesSquaredMatrix() :
    OldMassesSquaredFromMatrix< double >(),
    matrixElements()
  {
    // This constructor is just an initialization list.
  }

  OldRealMassesSquaredMatrix::~OldRealMassesSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the last scale which was used to update them.
  Eigen::MatrixXd OldRealMassesSquaredMatrix::CurrentValues(
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( fieldConfiguration );
      for( size_t columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex )
        = matrixElements[ rowsTimesLength + columnIndex ](
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the natural exponent of logarithmOfScale.
  Eigen::MatrixXd OldRealMassesSquaredMatrix::CurrentValues(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( fieldConfiguration,
                                                      logarithmOfScale );
      for( size_t columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex )
        = matrixElements[ rowsTimesLength + columnIndex ]( fieldConfiguration,
                                                           logarithmOfScale );
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

} /* namespace VevaciousPlusPlus */