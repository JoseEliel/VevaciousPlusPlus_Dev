/*
 * NodesOnPlanes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef NODESONPLANES_HPP_
#define NODESONPLANES_HPP_

#include "CommonIncludes.hpp"
#include "NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class NodesOnPlanes : public NodesFromParameterization
  {
  public:
    NodesOnPlanes( std::vector< double > const& falseVacuum,
                   std::vector< double > const& trueVacuum,
                   size_t const numberOfIntermediateNodes );
    virtual
    ~NodesOnPlanes();


    // This returns all the nodes based on the numbers given in
    // pathParameterization, ordered in the sequence that they are visited in
    // the path from the false vacuum to the true vacuum, with front() being
    // the false vacuum and back() being the true vacuum.
    std::vector< std::vector< double > >
    PathNodeSet( std::vector< double > const& pathParameterization ) const;

    // This returns the vector which would represent the node at
    // pathNodes[ adjustmentOrderIndex ] if it were to be set based on
    // nodeParameterization and the rest of the nodes in pathNodes, unless it
    // would mean moving the true vacuum or false vacuum, in which case an
    // exception is thrown. The node is set to be on a plane perpendicular to
    // two other nodes (decided by the derived class), by
    // nodeParameterization being taken to be a vector on a plane with the
    // field with index referenceField being zero, through
    // ProjectPerpendicularToAndShift.
    virtual std::vector< double >
    NodeProposal( size_t const adjustmentOrderIndex,
                  std::vector< double > const& nodeParameterization ) const;


  protected:
    size_t referenceField;
    size_t numberOfParametersPerNode;


    // This should add the perpendicular component from the parameterization
    // given by nodeParameterization along with startNode and endNode to
    // nodeVector.
    virtual void
    AddTransformedNode( std::vector< double > const& nodeVector,
                        std::vector< double > const& startNode,
                        std::vector< double > const& endNode,
                 std::vector< double > const& nodeParameterization ) const = 0;

    // This sets nodeVector to be the vector sum of startNode plus
    // shiftFraction times the difference between startNode and endNode, plus
    // the node given by nodeParameterization transformed to be perpendicular
    // to the difference between startNode and endNode.
    void TransformPerpendicularToAndShift( std::vector< double >& nodeVector,
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                             std::vector< double > const& nodeParameterization,
                                           double const shiftFraction ) const;

    // This should return the false-vacuum-side node of the pair of nodes
    // from which the node at adjustmentOrderIndex should be set.
    std::vector< double > const&
    FalseSideNode( size_t const adjustmentOrderIndex,
                   std::vector< std::vector< double > > const& nodeSet
                                                       = pathNodes ) const = 0;

    // This should return the true-vacuum-side node of the pair of nodes
    // from which the node at nodeIndex should be set.
    std::vector< double > const&
    TrueSideNode( size_t const adjustmentOrderIndex,
                  std::vector< std::vector< double > > const& nodeSet
                                                       = pathNodes ) const = 0;

    // This should return the fraction along the node difference vector that
    // the rotated plane should be shifted appropriate for
    // pathNodes[ nodeIndex ].
    double ShiftFraction( size_t const nodeIndex ) const = 0;
  };




  // This returns all the nodes based on the numbers given in
  // pathParameterization, ordered in the sequence that they are visited in
  // the path from the false vacuum to the true vacuum, with front() being
  // the false vacuum and back() being the true vacuum.
  std::vector< std::vector< double > > NodesOnPlanes::PathNodeSet(
                      std::vector< double > const& pathParameterization ) const
  {
    std::vector< std::vector< double > > nodeSet( pathNodes );
    std::vector< double > nodeParameterization( numberOfParametersPerNode );
    std::vector< double >::const_iterator
    currentParameterizationStart( pathParameterization.begin() );
    for( size_t nodeIndex( 1 );
         nodeIndex <= numberOfIntermediateNodes;
         ++nodeIndex )
    {
      nodeParameterization.assign( currentParameterizationStart,
                ( currentParameterizationStart + numberOfParametersPerNode ) );
      TransformPerpendicularToAndShift( nodeSet[ nodeIndex ],
                                        FalseSideNode( nodeIndex ),
                                        TrueSideNode( nodeIndex ),
                                        nodeParameterization,
                                        ShiftFraction( nodeIndex ) );
      currentParameterizationStart += numberOfParametersPerNode;
    }
    return nodeSet;
  }

  // This returns the vector which would represent the node at
  // pathNodes[ adjustmentOrderIndex ] if it were to be set based on
  // nodeParameterization and the rest of the nodes in pathNodes, unless it
  // would mean moving the true vacuum or false vacuum, in which case an
  // exception is thrown. The node is set to be on a plane perpendicular to
  // two other nodes (decided by the derived class), by
  // nodeParameterization being taken to be a vector on a plane with the
  // field with index referenceField being zero, through
  // ProjectPerpendicularToAndShift.
  inline std::vector< double >
  NodesOnPlanes::NodeProposal( size_t const adjustmentOrderIndex,
                      std::vector< double > const& nodeParameterization ) const
  {
    std::vector< double > returnVector( numberOfFields );
    TransformPerpendicularToAndShift( returnVector,
                                      FalseSideNode( adjustmentOrderIndex ),
                                      TrueSideNode( adjustmentOrderIndex ),
                                      nodeParameterization,
                                      ShiftFraction( adjustmentOrderIndex ) );
    return returnVector;
  }


  // This sets nodeVector to be the vector sum of startNode plus
  // shiftFraction times the difference between startNode and endNode, plus
  // the node given by nodeParameterization transformed to be perpendicular
  // to the difference between startNode and endNode.
  inline void NodesOnPlanes::TransformPerpendicularToAndShift(
                                             std::vector< double >& nodeVector,
                                        std::vector< double > const& startNode,
                                          std::vector< double > const& endNode,
                             std::vector< double > const& nodeParameterization,
                                             double const shiftFraction ) const
  {
    // We actually take the fraction along the straight line between nodes (the
    // "shift") before adding the perpendicular component.
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      nodeVector[ fieldIndex ]
      = ( ( 1.0 - shiftFraction ) * startNode[ fieldIndex ]
          + ( shiftFraction * endNode[ fieldIndex ] ) );
    }
    if( nodeParameterization != zeroParameterization )
    {
      AddTransformedNode( nodeVector,
                          startNode,
                          endNode,
                          nodeParameterization );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* NODESONPLANES_HPP_ */
