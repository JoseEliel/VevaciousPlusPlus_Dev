/*
 * MinimizingPotentialOnParallelPlanes.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnParallelPlanes.hpp"

namespace VevaciousPlusPlus
{
  MinimizingPotentialOnParallelPlanes::MinimizingPotentialOnParallelPlanes(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        numberOfPathSegments,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    notYetProvidedPath( true ),
    planeDifferenceFraction( 1.0
                             / static_cast< double > ( numberOfPathSegments ) )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnParallelPlanes::~MinimizingPotentialOnParallelPlanes()
  {
    // This does nothing.
  }



  // This minimizes the potential on a series of hyperplanes, all
  // perpendicular to the vector difference of the vacua. It tries to avoid
  // "zig-zag-iness" from Minuit2 coming to rest on a jittery line reasonably
  // far from the starting path by setting the starting point on each plane
  // to be the previous node plus the vector difference of previous node from
  // the node before it, with special cases for the first and last varying
  // nodes, of course. It ignores both arguments, and also sets
  // notYetProvidedPath to false.
  TunnelPath const* MinimizingPotentialOnParallelPlanes::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    SetParallelVector( pathNodes.front(),
                       pathNodes.back() );
    SetUpHouseholderReflection();
    ROOT::Minuit2::MnMigrad mnMigrad( *this,
                                      nodeZeroParameterization,
                                      minuitInitialSteps,
                                      minuitStrategy );

    currentHyperplaneOrigin = pathNodes.front();
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentHyperplaneOrigin[ fieldIndex ] += ( planeDifferenceFraction
                                    * currentParallelComponent[ fieldIndex ] );
    }
    RunMigradAndPutTransformedResultIn( mnMigrad,
                                        pathNodes[ 1 ] );

  }

} /* namespace VevaciousPlusPlus */
