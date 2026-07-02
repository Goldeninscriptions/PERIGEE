#ifndef PDNSOLUTION_NS_HPP
#define PDNSOLUTION_NS_HPP
// ==================================================================
// PDNSolution_NS.hpp
//
// This is the solution class for NS Solver. 
// 
// The solution vector has 4 dofs, pressure, u, v, w. 
//
// Author: Ju Liu
// Date Created: Jan 7 2020
// ==================================================================
#include "PDNSolution.hpp"
#include "FEANode.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_FarFieldInflowBC.hpp"

class PDNSolution_NS : public PDNSolution
{
  public:
    PDNSolution_NS( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_FarFieldInflowBC * const &infbc,
        const int &type, const double &uniform_speed = 0.0,
        bool isprint = true );

    PDNSolution_NS( const APart_Node * const &pNode,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc,
        const int &type, const double &uniform_speed = 0.0,
        bool isprint = true );

    PDNSolution_NS( const APart_Node * const &pNode, 
        const int &type, bool isprint = true );

    virtual ~PDNSolution_NS() = default;

  private:
    // case 0 : generate full zero solution
    void Init_zero( const APart_Node * const &pNode_ptr, bool isprint );

    // case 1: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_FarFieldInflowBC * const &infbc,
        bool isprint );

    void Init_flow_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr,
        const ALocal_InflowBC * const &infbc,
        bool isprint );

    // case 2: generate flow parabolic for an arbitrary inlet face
    //         with the unit flow rate.
    void Init_pipe_parabolic( const APart_Node * const &pNode_ptr,
        const FEANode * const &fNode_ptr, bool isprint );

    // case 3: generate a uniform inflow velocity on all inlet nodes.
    void Init_uniform_inflow( const APart_Node * const &pNode_ptr,
        const ALocal_FarFieldInflowBC * const &infbc,
        const double &uniform_speed, bool isprint );

    void Init_uniform_inflow( const APart_Node * const &pNode_ptr,
        const ALocal_InflowBC * const &infbc,
        const double &uniform_speed, bool isprint );
};

#endif
