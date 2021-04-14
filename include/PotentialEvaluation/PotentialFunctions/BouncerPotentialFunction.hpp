/*
 * PotentialFunction.hpp
 *
 *  Created on: Feb 25, 2021
 *      Author: Eliel Camargo-Molina
 */

#ifndef BOUNCERPOTENTIALFUNCTION_HPP_
#define BOUNCERPOTENTIALFUNCTION_HPP_
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include "PotentialMinimization/PotentialMinimum.hpp"


namespace VevaciousPlusPlus
{

  class BouncerPotentialFunction : public PotentialFunction 
  {
  public:
    BouncerPotentialFunction(double (*potential)(std::vector< double >, double), double renormalizationScale, 
      LagrangianParameterManager& lagrangianParameterManager);
 

    virtual ~BouncerPotentialFunction() {}

    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
    { return ( renormalizationScale * renormalizationScale ); }

    // This should return the energy density in GeV^4 of the potential for a
    // state strongly peaked around expectation values (in GeV) for the fields
    // given by the values of fieldConfiguration and temperature in GeV given
    // by temperatureValue.
    double operator()( std::vector< double > const& fieldConfiguration,
                               double const temperatureValue = 0.0 ) const 
    {
      return potential(fieldConfiguration, temperatureValue);
    }

    

  protected:
    double (*potential)(std::vector< double >, double);
    double renormalizationScale;

  };





} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFUNCTION_HPP_ */
