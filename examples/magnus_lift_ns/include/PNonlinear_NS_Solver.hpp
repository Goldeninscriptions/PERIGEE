#ifndef PNONLINEAR_NS_SOLVER_HPP
#define PNONLINEAR_NS_SOLVER_HPP
// ==================================================================
// PNonlinear_NS_Solver.hpp
// 
// Parallel nonlinear solver for Navier-Stokes equations. 
//
// Author: Ju Liu
// Date: Feb 11 2020
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
#include "PDNSolution_NS.hpp"
#include "ALocal_FarFieldInflowBC.hpp"
#include "ALocal_RotatedBC.hpp"

class PNonlinear_NS_Solver
{
  public:
    PNonlinear_NS_Solver( 
        std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
        std::unique_ptr<Matrix_PETSc> in_bc_mat,
        std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
        const double &input_freestream_speed,
        const double &input_freestream_thd_time,
        const double &input_nrtol, const double &input_natol, 
        const double &input_ndtol, const int &input_max_iteration, 
        const int &input_renew_freq, 
        const double &input_angular_velo,
        const double &input_angular_thd_time,
        const int &input_renew_threshold = 4 );

    ~PNonlinear_NS_Solver() = default;

    int get_non_max_its() const {return nmaxits;}

    void print_info() const;

    void print_lsolver_info() const {lsolver->print_info();}

    // --------------------------------------------------------------
    // GenAlpha_Solve_NS:
    // This is a solver for fluid dynamics.
    //
    // This solver solves the Navier-Stokes using 2nd-order Generalized
    // alpha method.
    // --------------------------------------------------------------
    int GenAlpha_Solve_NS(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_dot_sol,
        const PDNSolution * const &pre_sol,
        PDNSolution * const &dot_sol,
        PDNSolution * const &sol,
        const ALocal_FarFieldInflowBC * const &infnbc_part,
        const ALocal_RotatedBC * const &rotbc_part,
        const IGenBC * const &gbc,
        IPGAssem * const &gassem_ptr ) const;

  private:
    const double freestream_speed, nr_tol, na_tol, nd_tol;
    const double freestream_thd_time, angular_velo, angular_thd_time;
    const int nmaxits, nrenew_freq, nrenew_threshold;

    const std::unique_ptr<PLinear_Solver_PETSc> lsolver;
    const std::unique_ptr<Matrix_PETSc> bc_mat;
    const std::unique_ptr<TimeMethod_GenAlpha> tmga;

#ifdef PETSC_USE_LOG
    PetscLogEvent mat_assem_0_event, mat_assem_1_event;
    PetscLogEvent vec_assem_0_event, vec_assem_1_event;
    PetscClassId classid_assembly;
#endif

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {
      SYS_T::commPrint("  === NR ite: %d, r_error: %e, a_error: %e \n",
          count, rel_err, abs_err);
    }

    void update_rotating_wall_value(
        const double &stime,
        const ALocal_RotatedBC * const &rotbc,
        PDNSolution * const &sol ) const;

    void update_rotating_wall_dot_value(
        const double &stime,
        const ALocal_RotatedBC * const &rotbc,
        PDNSolution * const &dot_sol ) const;

    void update_uniform_inflow_value(
        const double &stime,
        const ALocal_FarFieldInflowBC * const &infbc,
        PDNSolution * const &sol ) const;

    void update_uniform_inflow_dot_value(
        const double &stime,
        const ALocal_FarFieldInflowBC * const &infbc,
        PDNSolution * const &dot_sol ) const;

    double get_ramped_freestream_speed( const double &stime ) const;

    double get_ramped_freestream_accel( const double &stime ) const;

    double get_ramped_angular_velocity( const double &stime ) const;

    double get_ramped_angular_accel( const double &stime ) const;

};

#endif
