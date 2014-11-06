/*
 * FixedScaleOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FIXEDSCALEONELOOPPOTENTIAL_HPP_
#define FIXEDSCALEONELOOPPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"
#include "SlhaManagement/RunningParameterManager.hpp"
#include "SlhaManagement/SlhaManager.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class FixedScaleOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    FixedScaleOneLoopPotential( std::string const& modelFilename,
                                double const scaleRangeMinimumFactor,
            bool const treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                               double const assumedPositiveOrNegativeTolerance,
                           RunningParameterManager& runningParameterManager );
    FixedScaleOneLoopPotential(
          PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses );
    virtual ~FixedScaleOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double operator()( std::vector< double > const& fieldConfiguration,
                               double const temperatureValue = 0.0 ) const;

    // This returns the tree-level potential energy density evaluated at the
    // correct scale.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 );

    // This returns the square of the renormalization scale.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
    { return ( currentMinimumRenormalizationScale
               * currentMinimumRenormalizationScale ); }

    // This performs all relevant updates for the new SLHA data except for
    // propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void UpdateSelfForNewSlha( SlhaManager const& slhaManager );

    virtual PolynomialGradientTargetSystem* HomotopyContinuationTargetSystem()
    { return &homotopyContinuationTargetSystem; }


  protected:
    double inverseRenormalizationScaleSquared;
    PolynomialGradientTargetSystem homotopyContinuationTargetSystem;
  };




  // This returns the tree-level potential energy density evaluated at the
  // correct scale.
  inline double FixedScaleOneLoopPotential::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    return treeLevelPotential( fieldConfiguration );
  }

  // This sets dsbFieldValueInputs based on the SLHA file just read in.
  inline void FixedScaleOneLoopPotential::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    currentMinimumRenormalizationScale = runningParameters.LowestBlockScale();
    squareOfMinimumRenormalizationScale = ( currentMinimumRenormalizationScale
                                        * currentMinimumRenormalizationScale );
    inverseRenormalizationScaleSquared
    = ( 1.0 / squareOfMinimumRenormalizationScale );
    currentMaximumRenormalizationScale = runningParameters.HighestBlockScale();
    if( currentMaximumRenormalizationScale
        < ( scaleRangeMinimumFactor * currentMinimumRenormalizationScale ) )
    {
      currentMaximumRenormalizationScale
      = ( scaleRangeMinimumFactor * currentMinimumRenormalizationScale );
    }
    squareOfMaximumRenormalizationScale = ( currentMaximumRenormalizationScale
                                        * currentMaximumRenormalizationScale );
    std::vector< double > fieldOrigin( numberOfFields,
                                       0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      dsbFieldValueInputs[ fieldIndex ]
      = dsbFieldValuePolynomials[ fieldIndex ]( fieldOrigin );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* FIXEDSCALEONELOOPPOTENTIAL_HPP_ */
