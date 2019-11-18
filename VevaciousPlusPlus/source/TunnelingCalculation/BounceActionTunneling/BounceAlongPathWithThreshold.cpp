/*
 * BounceAlongPathWithThreshold.cpp
 *
 *  Created on: Jul 1, 2014
 *      Authors: Ben O'Leary (benjamin.oleary@gmail.com)
 *               José Eliel Camargo-Molina (eliel@camargo-molina.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"

namespace VevaciousPlusPlus
{

    BounceAlongPathWithThreshold::BounceAlongPathWithThreshold(
            std::vector< std::unique_ptr<BouncePathFinder> > pathFinders,
            std::unique_ptr<BounceActionCalculator> actionCalculator,
            TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
            bool const pathDeformation,
            double const survivalProbabilityThreshold,
            unsigned int const thermalIntegrationResolution,
            unsigned int const temperatureAccuracy,
            unsigned int const pathPotentialResolution,
            unsigned int const pathFindingTimeout,
            double const vacuumSeparationFraction ) :
            BounceActionTunneler( tunnelingStrategy,
                                  survivalProbabilityThreshold,
                                  temperatureAccuracy,
                                  vacuumSeparationFraction ),
            pathFinders( std::move(pathFinders) ),
            actionCalculator( std::move(actionCalculator) ),
            pathDeformation( pathDeformation ),
            thermalIntegrationResolution( thermalIntegrationResolution ),
            pathPotentialResolution( pathPotentialResolution ),
            pathFindingTimeout( pathFindingTimeout )
    {
      // This constructor is just an initialization list.
    }

    BounceAlongPathWithThreshold::~BounceAlongPathWithThreshold()
    {
    }


    // This sets thermalSurvivalProbability by numerically integrating up to the
    // critical temperature for tunneling to be possible from T = 0 unless the
    // integral already passes a threshold, and sets
    // dominantTemperatureInGigaElectronVolts to be the temperature with the
    // lowest survival probability.
    void BounceAlongPathWithThreshold::ContinueThermalTunneling(
            PotentialFunction const& potentialFunction,
            PotentialMinimum const& falseVacuum,
            PotentialMinimum const& trueVacuum,
            double const potentialAtOriginAtZeroTemperature )
    {
      // First we set up the (square of the) threshold distance that we demand
      // between the vacua at every temperature to trust the tunneling
      // calculation.
      double const thresholdSeparationSquared( vacuumSeparationFractionSquared
      * falseVacuum.SquareDistanceTo( trueVacuum ));

      // We sum up decay widths over increasing temperatures.
      double partialDecayWidth( 0.0 );
      // The partial decay width scaled by the volume of the observable Universe
      // is recorded in partialDecayWidth so that the bounce action threshold for
      // each temperature can be calculated taking into account the contributions
      // from higher temperatures.
      double const temperatureStep( rangeOfMaxTemperatureForOriginToTrue.first
                                    / static_cast< double >( thermalIntegrationResolution + 1 ) );
      double currentTemperature( 0.0 );
      MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
      thermalPotentialMinimizer.SetTemperature( currentTemperature );
      PotentialMinimum thermalFalseVacuum( falseVacuum );
      PotentialMinimum thermalTrueVacuum( trueVacuum );
      double const thresholdDecayWidth( -log( survivalProbabilityThreshold )
                                        / ( temperatureStep * exp( lnOfThermalIntegrationFactor ) ) );

      double smallestExponent( maximumPowerOfNaturalExponent );
      dominantTemperatureInGigaElectronVolts = 0.0;

      for( unsigned int whichStep( 0 );
           whichStep < thermalIntegrationResolution;
           ++whichStep )
      {
        currentTemperature += temperatureStep;
        thermalPotentialMinimizer.SetTemperature( currentTemperature );
        // We update the positions of the thermal vacua based on their positions
        // at the last temperature step.
        thermalFalseVacuum
                = thermalPotentialMinimizer( thermalFalseVacuum.FieldConfiguration() );
        // We have to keep checking to see if the field origin should be the
        // thermal false vacuum. The result of thermalPotentialMinimizer already
        // has the value of the potential at the field origin subtracted, so we
        // just compare with zero.
        if( thermalFalseVacuum.PotentialValue() > 0.0 )
        {
          thermalFalseVacuum
                  = PotentialMinimum( potentialFunction.FieldValuesOrigin(),
                                      0.0 );
        }
        thermalTrueVacuum
                = thermalPotentialMinimizer( thermalTrueVacuum.FieldConfiguration() );

        if( !( thermalTrueVacuum.FunctionValue()
               < thermalFalseVacuum.FunctionValue() ) )
        {
          // If the thermal vacua are the wrong way around in depth order
          // (possibly due to the thermal minimizer failing to converge properly)
          // then we break without adding in any decay width for this
          // temperature.
          std::stringstream warningBuilder;
          warningBuilder << "At temperature " << currentTemperature
                         << " GeV, minimizer rolled from panic vacuum from a lower temperature"
                         << " to configuration which is not deeper than the configuration to"
                         << " which it rolled from the DSB vacuum from that lower temperature."
                         << " Skipping the contribution of this temperature.";
          WarningLogger::LogWarning( warningBuilder.str() );
          continue;
        }
        else if( thermalTrueVacuum.SquareDistanceTo( thermalFalseVacuum )
                 < thresholdSeparationSquared )
        {
          // If the thermal vacua have gotten so close that a tunneling
          // calculation is suspect, we break without adding in any decay width
          // for this temperature.
          std::stringstream warningBuilder;
          warningBuilder << "At temperature " << currentTemperature
                         << " GeV, minimizer found DSB vacuum and panic vacuum to be so close"
                         << " that a tunneling calculation is not trustworthy. Skipping the"
                         << " contribution of this temperature and higher temperatures.";
          WarningLogger::LogWarning( warningBuilder.str() );
          break;
        }

        double const actionThreshold( -currentTemperature
                                      * log( currentTemperature
                                             * currentTemperature
                                             * ( thresholdDecayWidth
                                                 - partialDecayWidth ) ) );
        double const
                bounceOverTemperature( BoundedBounceAction( potentialFunction,
                                                            thermalFalseVacuum,
                                                            thermalTrueVacuum,
                                                            currentTemperature,
                                                            actionThreshold,
                                                            thresholdSeparationSquared )
                                       / currentTemperature );

        if( bounceOverTemperature < maximumPowerOfNaturalExponent )
        {
          partialDecayWidth += ( exp( -bounceOverTemperature )
                                 / ( currentTemperature * currentTemperature ) );
        }
        if( bounceOverTemperature < smallestExponent )
        {
          smallestExponent = bounceOverTemperature;
          dominantTemperatureInGigaElectronVolts = currentTemperature;
        }

        if( partialDecayWidth > thresholdDecayWidth )
        {
          // We don't bother calculating the rest of the contributions to the
          // integral of the decay width if it is already large enough that the
          // survival probability is below the threshold.
          break;
        }

      }
      if( partialDecayWidth > 0.0 )
      {
        logOfMinusLogOfThermalProbability = ( lnOfThermalIntegrationFactor
                                              + log( partialDecayWidth * temperatureStep ) );
      }
      else
      {
        logOfMinusLogOfThermalProbability
                = -exp( maximumPowerOfNaturalExponent );
        std::stringstream warningBuilder;
        warningBuilder
                << "The calculated integrated thermal decay width was so close to zero"
                << " that taking its logarithm would be problematic, so setting the"
                << " logarithm of the negative of the logarithm of the thermal survival"
                << " probability to " << logOfMinusLogOfThermalProbability << ".";
        WarningLogger::LogWarning( warningBuilder.str() );
      }
      SetThermalSurvivalProbability();
    }

    // This returns either the dimensionless bounce action integrated over four
    // dimensions (for zero temperature) or the dimensionful bounce action
    // integrated over three dimensions (for non-zero temperature) for tunneling
    // from falseVacuum to trueVacuum at temperature tunnelingTemperature, or an
    // upper bound if the upper bound drops below actionThreshold during the
    // course of the calculation. The vacua are assumed to already be the minima
    // at tunnelingTemperature.
    double BounceAlongPathWithThreshold::BoundedBounceAction(
            PotentialFunction const& potentialFunction,
            PotentialMinimum const& falseVacuum,
            PotentialMinimum const& trueVacuum,
            double const tunnelingTemperature,
            double const actionThreshold,
            double const requiredVacuumSeparationSquared )
    {
      std::vector< std::vector< double > > straightPath( 2,
                                                         falseVacuum.FieldConfiguration() );
      straightPath.back() = trueVacuum.FieldConfiguration();

      LinearSplineThroughNodes bestPathInstance( straightPath,
                                std::vector< double >( 0 ),
                                tunnelingTemperature );

      std::cout << "(EC)              pointer of bestPathInstance after creation is " << &bestPathInstance << std::endl;


      std::shared_ptr< const TunnelPath> bestPath = std::make_shared< const LinearSplineThroughNodes>(bestPathInstance);



      //std::shared_ptr< const TunnelPath>  bestPath( new LinearSplineThroughNodes( straightPath,
      //                                                                            std::vector< double >( 0 ),
      //                                                                            tunnelingTemperature ) );

      std::cout << "(EC)              raw pointer of bestPath after creation is " << bestPath.get()
                << std::endl;

      actionCalculator->ResetVacua( potentialFunction,
                                    falseVacuum,
                                    trueVacuum,
                                    tunnelingTemperature );
      //std::cout<<"(JR) after ResetVacua " << std::endl;
//    std::cout<<"(JR) pathPotentialResolution "<< pathPotentialResolution << std::endl;

      SplinePotential pathPotential( potentialFunction,
                                     *bestPath,
                                     pathPotentialResolution,
                                     requiredVacuumSeparationSquared );
      //std::cout<<"(JR) after SplinePotential " << std::endl;

      if( !(pathPotential.EnergyBarrierWasResolved()) )
      {
//      std::cout<<"(JR) in if  " << std::endl;
        std::stringstream warningBuilder;
        warningBuilder << "Unable to resolve an energy barrier between false"
                       << " vacuum and true vacuum: returning bounce action of zero (which"
                       << " should be sufficient to exclude the parameter point).";
        WarningLogger::LogWarning( warningBuilder.str() );
        bestPath.reset();
        return 0.0;
      }

      //UndershootOvershootBubble* bestBubbleInstance = dynamic_cast<UndershootOvershootBubble*>((*actionCalculator)( *bestPath, pathPotential ));

      //std::shared_ptr< const BubbleProfile> bestBubble = std::make_shared<const UndershootOvershootBubble>(*bestBubbleInstance);

      std::shared_ptr< const BubbleProfile> bestBubble = (*actionCalculator)( *bestPath, pathPotential );

      std::cout << "(EC)              raw pointer of bestBubble after creation is " << bestBubble.get()
                << std::endl;

      std::cout << std::endl
                << "Initial path bounce action = " << bestBubble->BounceAction();
      if( bestPath->NonZeroTemperature() )
      {
        std::cout << " GeV";

        thermalThresholdAndActions.push_back(actionThreshold);
        thermalThresholdAndActions.push_back(bestBubble->BounceAction());
      }
      else
      {
        thresholdAndActions.push_back(actionThreshold);
        thresholdAndActions.push_back(bestBubble->BounceAction());
      }
      std::cout << ", threshold is " << actionThreshold;
      if( bestPath->NonZeroTemperature() )
      {
        std::cout << " GeV";
      }
      std::cout << ".";
      std::cout << std::endl;

      // Checking if initial path already has a very low action

      if( bestBubble->BounceAction() < actionThreshold )
      {
        std::cout
                << std::endl
                << "Bounce action dropped below threshold, breaking off from looking"
                << " for further path improvements.";
        std::cout << std::endl;
        double const bounceAction( bestBubble->BounceAction() );
        bestBubble.reset();
        bestPath.reset();
        return bounceAction;
      }

      if( !pathDeformation )
      {
        std::cout
                << std::endl
                << "No path deformation. Taking straight path action.";
        std::cout << std::endl;
        std::cout
                << std::endl
                << pathDeformation ;
        std::cout << std::endl;
        double const bounceAction( bestBubble->BounceAction() );
        bestBubble.reset();
        bestPath.reset();
        return bounceAction;
      }

      // Declaring variables for timing
      time_t pathFindingStartTime;
      time_t currentTime;
      // Setting the starting time of the path finding
      time( &pathFindingStartTime );




      // (JR) moved outside of for loop to be able to set to null after , attempt to fix mem leak
      std::cout<<" ====== (JR) this is a memory leak fix attempt # 7 ==== " << std::endl;
      std::shared_ptr< const TunnelPath> currentPath(bestPath);
      std::shared_ptr< const BubbleProfile> currentBubble(bestBubble);

      for( auto pathFinder( pathFinders.begin() );
           pathFinder < pathFinders.end();
           ++pathFinder )
      {
        time(&currentTime);

        if (difftime( currentTime, pathFindingStartTime ) > pathFindingTimeout )
        {
          std::stringstream errorBuilder;
          errorBuilder << "Path finder has been running longer than the specified timeout of " << pathFindingTimeout << " seconds.";
          bestBubble.reset();
          bestPath.reset();
          currentBubble.reset();
          currentPath.reset();
          throw std::runtime_error( errorBuilder.str() );
          break;
        };

        std::cout << std::endl
                  << "Passing best path so far to next path finder.";
        std::cout << std::endl;

        (*pathFinder)->SetPotentialAndVacuaAndTemperature( potentialFunction,
                                                           falseVacuum,
                                                           trueVacuum,
                                                           tunnelingTemperature );
        std::cout << std::endl
                  << "(EC)              Before assignment of currentPath and currentBubble.";
        std::cout << std::endl;


        std::cout << std::endl
                  << "(EC)              After assignment of currentPath and currentBubble.";
        std::cout << std::endl;

        std::cout << "(EC)              raw pointer of currentPath after assignment is " << currentPath.get()
                  << std::endl;
        std::cout << "(EC)              raw pointer of currentBubble after assignment is " << currentBubble.get()
                  << std::endl;


        // The paths produced in sequence by pathFinder are kept separate from
        // bestPath to give more freedom to pathFinder internally (though I
        // cannot right now think of any way in which it would actually be
        // useful, apart from maybe some crazy MCMC which might try to get out of
        // a local minimum, which wouldn't work if it was sent back to the local
        // minimum at each step).


        // This loop will get a path from pathFinder and then repeat if
        // pathFinder decides that the path can be improved once the bubble
        // profile is obtained, as long as the bounce action has not dropped
        // below the threshold.

        do
        {
          // The nextPath and nextBubble pointers are not strictly necessary,
          // but they make the logic of the code clearer and will probably be
          // optimized away by the compiler anyway.
          time(&currentTime);

          if (difftime( currentTime, pathFindingStartTime ) > pathFindingTimeout ) // HERE THE CUTOFF IS SET, MAKE A VARIABLE IN INPUT
          {
            std::stringstream errorBuilder;
            errorBuilder << "Path finder has been running longer than the specified timeout of " << pathFindingTimeout << "seconds.";
            throw std::runtime_error( errorBuilder.str() );
            break;
          };

          std::cout << std::endl
                    << "(EC)              Before creating of nextPath and nextBubble.";
          std::cout << std::endl;

          std::shared_ptr<const TunnelPath> nextPath = (*pathFinder)->TryToImprovePath( *currentPath, *currentBubble ) ;

          std::cout << std::endl
                    << "(EC)             Adress of nextPath after creation is "<< nextPath.get();
          std::cout << std::endl;

          std::cout << std::endl
                    << "(EC)             Before creating Spline Potential.";
          std::cout << std::endl;

          SplinePotential potentialApproximation( potentialFunction,
                                                  *nextPath,
                                                  pathPotentialResolution,
                                                  requiredVacuumSeparationSquared );
          std::cout << std::endl
                    << "(EC)             After creating Spline Potential.";
          std::cout << std::endl;

          //UndershootOvershootBubble* nextBubbleInstance = dynamic_cast<UndershootOvershootBubble*>((*actionCalculator)( *nextPath, potentialApproximation ));

          std::cout << std::endl
                    << "(EC)             bextBubbleInstance is created";
          std::cout << std::endl;
          //std::cout << std::endl << "(EC)             address after creation is " << nextBubbleInstance <<std::endl;

          //std::shared_ptr< const BubbleProfile> nextBubble = std::make_shared<const UndershootOvershootBubble>(*nextBubbleInstance);

          std::shared_ptr< const BubbleProfile> nextBubble = (*actionCalculator)( *nextPath, potentialApproximation );

          std::cout << std::endl
                    << "(EC)              After creating of nextPath and nextBubble.";
          std::cout << std::endl;

          std::cout << "(EC)              raw pointer of nextPath is " << nextPath.get()
                    << std::endl;
          std::cout << "(EC)              raw pointer of nextBubble is " << nextBubble.get()
                    << std::endl;


          if( nextBubble->BounceAction() < bestBubble->BounceAction() )
          {
            // If nextBubble was an improvement on bestBubble, then bestPath
            // is set to point at nextPath, with the corresponding operations for
            // the bubble pointers.

            bestPath = nextPath;
            bestBubble = nextBubble;

            std::cout << "(EC)              bestPath now points to " << bestPath.get()
                      << std::endl;
            std::cout << "(EC)              bestBubble now points to " << bestBubble.get()
                      << std::endl;
          }


          currentBubble = nextBubble;
          currentPath = nextPath;


          std::cout << std::endl
                    << "bounce action for new path = " << currentBubble->BounceAction();
          if( currentPath->NonZeroTemperature() )
          {
            std::cout << " GeV";
          }
          std::cout << ", lowest bounce action so far = "
                    << bestBubble->BounceAction();
          if( currentPath->NonZeroTemperature() )
          {
            std::cout << " GeV";
          }
          std::cout << ", threshold is " << actionThreshold;
          if( currentPath->NonZeroTemperature() )
          {
            std::cout << " GeV";
          }
          std::cout << ".";
          std::cout << std::endl;

          std::cout<<"       (JR) nextPath managing "<< nextPath.get_cout() << "" << std::endl;
          std::cout<<"       (JR) nextBubble managing "<< nextBubble.get_cout() << "" << std::endl;
          nextPath.reset();
          nextBubble.reset();
          std::cout<<"       (JR) nextPath managing "<< nextPath.get_cout() << " after reset" << std::endl;
          std::cout<<"       (JR) nextBubble managing "<< nextBubble.get_cout() << " after reset" << std::endl;
        
        } 

        while( ( bestBubble->BounceAction() > actionThreshold )
                 &&
                 (*pathFinder)->PathCanBeImproved( *currentBubble ) );

        // Recodring the best action for each pathfinder
        if( bestPath->NonZeroTemperature() )
        {
          thermalThresholdAndActions.push_back(bestBubble->BounceAction());
        }
        else
        {
          thresholdAndActions.push_back(bestBubble->BounceAction());
        }

        // We don't bother with the rest of the path finders if the action has
        // already dropped below the threshold.
        if( bestBubble->BounceAction() < actionThreshold )
        {
          std::cout
                  << std::endl
                  << "Bounce action dropped below threshold, breaking off from looking"
                  << " for further path improvements.";
          std::cout << std::endl;
          break;
        }
      }

      std::cout<<"       (JR) currentPath managing "<< currentPath.get_cout() << "" << std::endl;
      std::cout<<"       (JR) currentBubble managing "<< currentBubble.get_cout() << "" << std::endl;
      currentPath.reset();
      currentBubble.reset();
      std::cout<<"       (JR) currentPath managing "<< currentPath.get_cout() << " after reset" << std::endl;
      std::cout<<"       (JR) currentBubble managing "<< currentBubble.get_cout() << " after reset" << std::endl;

      std::cout << std::endl
                << "Lowest path bounce action at " << tunnelingTemperature << " GeV was "
                << bestBubble->BounceAction();
      if( bestPath->NonZeroTemperature() )
      {
        std::cout << " GeV";
      }
      std::cout << ", threshold is " << actionThreshold;
      if( bestPath->NonZeroTemperature() )
      {
        std::cout << " GeV";
      }
      std::cout << ".";
      std::cout << std::endl;

      double const bounceAction( bestBubble->BounceAction() );

      bestPath.reset();
      bestBubble.reset();
      std::cout<<"       (JR) END OF FUNCTION bestBubble managing "<< bestBubble.use_count() << " after reset" << std::endl;
      std::cout<<"       (JR) END OF FUNCTION bestPath managing "<< bestPath.use_count() << " after reset" << std::endl;
      return bounceAction;
    }

} /* namespace VevaciousPlusPlus */

