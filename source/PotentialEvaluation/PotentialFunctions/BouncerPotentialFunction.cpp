/*
 * TreeLevelPotential.cpp
 *
 *  Created on: Oct 12, 2020
 *      Author: José Eliel Camargo-Molina (Eliel@camargo-molina.com) 
 */

#include "PotentialEvaluation/PotentialFunctions/BouncerPotentialFunction.hpp"

namespace VevaciousPlusPlus
{

BouncerPotentialFunction::BouncerPotentialFunction(double (*potential)(std::vector< double >, double), double renormalizationScale, 
      LagrangianParameterManager& lagrangianParameterManager):
   PotentialFunction(lagrangianParameterManager),
   renormalizationScale(renormalizationScale),
   potential(potential){}
} /* namespace VevaciousPlusPlus */
