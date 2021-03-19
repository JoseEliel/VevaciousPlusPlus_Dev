/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 22 2021
 *      Authors: 
 *               José Eliel Camargo-Molina (elielcamargomolina@gmail.com)
 */

#include "VevaciousPlusPlus.hpp"
#include "LHPC/Utilities/RestrictedXmlParser.hpp"
#include "Utilities/FilePlaceholderManager.hpp"
#include "PotentialEvaluation/PotentialFunctions/BouncerPotentialFunction.hpp"

int main( int argumentCount,
          char** argumentCharArrays )
{
  std::string tunnelingStrategy = "QuantumThenThermal";
  double survivalProbabilityThreshold = 0.01;
  int thermalIntegrationResolution = 5;
  int temperatureAccuracy = 7;
  int resolutionOfPathPotential = 100;
  int pathFindingTimeout = 10000000;
  double vacuumSeparationFraction = 0.2;
  
  std::string bouncePotentialFitClass = "BubbleShootingOnPathInFieldSpace";
  std::string bouncePotentialFitArguments(
  "<NumberShootAttemptsAllowed>"
            "32"
          "</NumberShootAttemptsAllowed>"
          "<RadialResolution>"
          "  0.05"
          "</RadialResolution>");


  std::string tunnelPathFinders(
"    <TunnelPathFinders>"
"        <PathFinder>"
"          <ClassType>"
"            MinuitOnPotentialOnParallelPlanes"
"          </ClassType>"
"          <ConstructorArguments>"
"            <NumberOfPathSegments>"
"              50"
"            </NumberOfPathSegments>"
"            <MinuitStrategy>"
"              1"
"            </MinuitStrategy>"
"            <MinuitTolerance>"
"              0.5"
"            </MinuitTolerance>"
"          </ConstructorArguments>"
"        </PathFinder>"
"        <PathFinder>"
"          <ClassType>"
"            MinuitOnPotentialPerpendicularToPath"
"          </ClassType>"
"          <ConstructorArguments>"
"            <NumberOfPathSegments>"
"              100"
"            </NumberOfPathSegments>"
"            <NumberOfAllowedWorsenings>"
"              1"
"            </NumberOfAllowedWorsenings>"
"            <ConvergenceThresholdFraction>"
"              0.05"
"            </ConvergenceThresholdFraction>"
"            <MinuitDampingFraction>"
"              0.75"
"            </MinuitDampingFraction>"
"            <NeighborDisplacementWeights>"
"              0.5"
"              0.25"
"            </NeighborDisplacementWeights>"
"            <MinuitStrategy>"
"              1"
"            </MinuitStrategy>"
"            <MinuitTolerance>"
"              0.5"
"            </MinuitTolerance>"
"          </ConstructorArguments>"
"        </PathFinder>"
"      </TunnelPathFinders>"
      );

  std::vector< std::unique_ptr<VevaciousPlusPlus::BouncePathFinder>> 
  pathFinders = std::move(VevaciousPlusPlus::VevaciousPlusPlus::CreateBouncePathFinders( tunnelPathFinders ));


  std::unique_ptr<VevaciousPlusPlus::BounceActionCalculator> bounceActionCalculator(  std::move(
                       VevaciousPlusPlus::VevaciousPlusPlus::CreateBounceActionCalculator( bouncePotentialFitClass,
                                               bouncePotentialFitArguments ) ) );

  std::unique_ptr<VevaciousPlusPlus::BounceAlongPathWithThreshold> bouncerAlongPathWithThreshold = VevaciousPlusPlus::Utils::make_unique<VevaciousPlusPlus::BounceAlongPathWithThreshold>( std::move(pathFinders),
                                             std::move(bounceActionCalculator),
                               VevaciousPlusPlus::VevaciousPlusPlus::InterpretTunnelingStrategy( tunnelingStrategy ),
                                             survivalProbabilityThreshold,
                                             thermalIntegrationResolution,
                                             temperatureAccuracy,
                                             resolutionOfPathPotential,
                                             pathFindingTimeout,
                                             vacuumSeparationFraction );

// Here I create the LagrangianParameterManager I need to construct the potential

 std::unique_ptr<VevaciousPlusPlus::LesHouchesAccordBlockEntryManager> LagrangianParameterManager
    = std::move(VevaciousPlusPlus::VevaciousPlusPlus::CreateLagrangianParameterManager( "LesHouchesAccordBlockEntryManager", "Empty" ));

// Potential function 




  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}
