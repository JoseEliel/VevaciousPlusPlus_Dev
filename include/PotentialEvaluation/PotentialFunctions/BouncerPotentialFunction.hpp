/*
 * PotentialFunction.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFUNCTION_HPP_
#define POTENTIALFUNCTION_HPP_

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
    BouncerPotentialFunction(const size_t numberOfFields, 
                              std::vector<std::string> fieldNames, 
                              std::vector< double > dsbFieldValueInputs,
                              double (*potential)(std::vector< double >, double), 
                              LagrangianParameterManager& lagrangianParameterManager);
 

    virtual ~BouncerPotentialFunction() {}

    // This should return the energy density in GeV^4 of the potential for a
    // state strongly peaked around expectation values (in GeV) for the fields
    // given by the values of fieldConfiguration and temperature in GeV given
    // by temperatureValue.
    virtual double operator()( std::vector< double > const& fieldConfiguration,
                               double const temperatureValue = 0.0 ) const 
    {
      return potential(fieldConfiguration, temperatureValue);
    }


  protected:
    double (*potential)(std::vector< double >);

  };





} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFUNCTION_HPP_ */
