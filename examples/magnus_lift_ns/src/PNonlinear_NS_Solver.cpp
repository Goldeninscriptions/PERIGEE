#include "PNonlinear_NS_Solver.hpp"
#include "LoadData.hpp"
#include <cmath>

PNonlinear_NS_Solver::PNonlinear_NS_Solver(
    std::unique_ptr<PLinear_Solver_PETSc> in_lsolver,
    std::unique_ptr<Matrix_PETSc> in_bc_mat,
    std::unique_ptr<TimeMethod_GenAlpha> in_tmga,
    const double &input_freestream_speed,
    const double &input_freestream_thd_time,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, 
    const int &input_renew_freq,
    const double &input_angular_velo,
    const double &input_angular_thd_time,
    const int &input_renew_threshold )
: freestream_speed(input_freestream_speed), nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  freestream_thd_time(input_freestream_thd_time), angular_velo(input_angular_velo),
  angular_thd_time(input_angular_thd_time),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold),
  lsolver(std::move(in_lsolver)),
  bc_mat(std::move(in_bc_mat)),
  tmga(std::move(in_tmga))
{
#ifdef PETSC_USE_LOG
  PetscClassIdRegister("mat_vec_assembly", &classid_assembly);
  PetscLogEventRegister("assembly mat 0", classid_assembly, &mat_assem_0_event);
  PetscLogEventRegister("assembly mat 1", classid_assembly, &mat_assem_1_event);
  PetscLogEventRegister("assembly vec 0", classid_assembly, &vec_assem_0_event);
  PetscLogEventRegister("assembly vec 1", classid_assembly, &vec_assem_1_event);
#endif
}

void PNonlinear_NS_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Nonlinear solver setted up:\n");
  SYS_T::commPrint("  relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("  absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("  divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("  maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("  tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::commPrint("  tangent matrix renew threshold: %d \n", nrenew_threshold);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

int PNonlinear_NS_Solver::GenAlpha_Solve_NS(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    const ALocal_InflowBC * const &infnbc_part,
    const ALocal_RotatedBC * const &rotbc_part,
    const IGenBC * const &gbc,
    IPGAssem * const &gassem_ptr ) const
{
  // Initialize the counter and error
  int nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga->get_gamma();
  const double alpha_m = tmga->get_alpha_m();
  const double alpha_f = tmga->get_alpha_f();

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );

  // ------------------------------------------------- 
  // Update the inflow boundary values
  update_uniform_inflow_value(curr_time + dt, infnbc_part, sol);
  update_uniform_inflow_value(curr_time + alpha_f * dt, infnbc_part, &sol_alpha);
  // ------------------------------------------------- 

  // ------------------------------------------------- 
  // Update rotating wall velocity values
  update_rotating_wall_value(curr_time + dt, rotbc_part, sol);
  update_rotating_wall_value(curr_time + alpha_f * dt, rotbc_part, &sol_alpha);
  // -------------------------------------------------

  // ------------------------------------------------- 
  // Update the dot_inflow boundary values
  update_uniform_inflow_dot_value(curr_time + dt, infnbc_part, dot_sol);
  update_uniform_inflow_dot_value(curr_time + alpha_m * dt, infnbc_part, &dot_sol_alpha);
  update_rotating_wall_dot_value(curr_time + dt, rotbc_part, dot_sol);
  update_rotating_wall_dot_value(curr_time + alpha_m * dt, rotbc_part, &dot_sol_alpha);
  // ------------------------------------------------- 

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol, 
        curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    SYS_T::commPrint("  --- M updated");
    
    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver->SetOperator( gassem_ptr->K );
  }
  else
  {
    gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  auto dot_step = SYS_T::make_unique<PDNSolution>( pre_sol );

  // Now do consistent Newton-Raphson iteration
  do
  {
    // solve the equation K dot_step = G
    lsolver->Solve( gassem_ptr->G, dot_step.get() );

    bc_mat->MatMultSol( dot_step.get() );

    nl_counter += 1;

    dot_sol->PlusAX( dot_step.get(), -1.0 );
    sol->PlusAX( dot_step.get(), (-1.0) * gamma * dt );

    dot_sol_alpha.PlusAX( dot_step.get(), (-1.0) * alpha_m );
    sol_alpha.PlusAX( dot_step.get(), (-1.0) * alpha_f * gamma * dt );

    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif

      SYS_T::commPrint("  --- M updated");
      lsolver->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    
    SYS_T::print_fatal_if( residual_norm != residual_norm, "Error: nonlinear solver residual norm is NaN. Job killed.\n" );
    
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  return nl_counter;
}

double PNonlinear_NS_Solver::get_ramped_freestream_speed( const double &stime ) const
{
  if( freestream_thd_time <= 0.0 || stime >= freestream_thd_time ) return freestream_speed;
  if( stime <= 0.0 ) return 0.0;

  return 0.5 * freestream_speed
    * ( 1.0 - std::cos( MATH_T::PI * stime / freestream_thd_time ) );
}

double PNonlinear_NS_Solver::get_ramped_freestream_accel( const double &stime ) const
{
  if( freestream_thd_time <= 0.0 || stime <= 0.0 || stime >= freestream_thd_time ) return 0.0;

  return 0.5 * freestream_speed * MATH_T::PI / freestream_thd_time
    * std::sin( MATH_T::PI * stime / freestream_thd_time );
}

double PNonlinear_NS_Solver::get_ramped_angular_velocity( const double &stime ) const
{
  if( angular_thd_time <= 0.0 || stime >= angular_thd_time ) return angular_velo;
  if( stime <= 0.0 ) return 0.0;

  return 0.5 * angular_velo
    * ( 1.0 - std::cos( MATH_T::PI * stime / angular_thd_time ) );
}

double PNonlinear_NS_Solver::get_ramped_angular_accel( const double &stime ) const
{
  if( angular_thd_time <= 0.0 || stime <= 0.0 || stime >= angular_thd_time ) return 0.0;

  return 0.5 * angular_velo * MATH_T::PI / angular_thd_time
    * std::sin( MATH_T::PI * stime / angular_thd_time );
}

void PNonlinear_NS_Solver::update_rotating_wall_value(
    const double &stime,
    const ALocal_RotatedBC * const &rotbc,
    PDNSolution * const &sol ) const
{
  constexpr double cylinder_radius = 0.5;
  const double ramped_angular_velocity = get_ramped_angular_velocity(stime);
  const int numnode = rotbc->get_Num_LD();

  for( int ii=0; ii<numnode; ++ii )
  {
    const int node_index = rotbc->get_LDN(ii);
    const Vector_3 pt = rotbc->get_LDN_pt_xyz(ii);
    const double radial_distance = std::sqrt( pt.x() * pt.x() + pt.y() * pt.y() );

    SYS_T::print_fatal_if( radial_distance <= 1.0e-12,
        "Error: rotating wall node is too close to the z-axis.\n" );

    const double wall_speed = ramped_angular_velocity * cylinder_radius;
    const double vals[3] = {
      -wall_speed * pt.y() / radial_distance,
       wall_speed * pt.x() / radial_distance,
       0.0
    };
    const int sol_idx[3] = { node_index * 4 + 1, node_index * 4 + 2, node_index * 4 + 3 };

    VecSetValues(sol->solution, 3, sol_idx, vals, INSERT_VALUES);
  }

  sol->Assembly_GhostUpdate();
}

void PNonlinear_NS_Solver::update_rotating_wall_dot_value(
    const double &stime,
    const ALocal_RotatedBC * const &rotbc,
    PDNSolution * const &dot_sol ) const
{
  constexpr double cylinder_radius = 0.5;
  const double ramped_angular_accel = get_ramped_angular_accel(stime);
  const int numnode = rotbc->get_Num_LD();

  for( int ii=0; ii<numnode; ++ii )
  {
    const int node_index = rotbc->get_LDN(ii);
    const Vector_3 pt = rotbc->get_LDN_pt_xyz(ii);
    const double radial_distance = std::sqrt( pt.x() * pt.x() + pt.y() * pt.y() );

    SYS_T::print_fatal_if( radial_distance <= 1.0e-12,
        "Error: rotating wall node is too close to the z-axis.\n" );

    const double wall_accel = ramped_angular_accel * cylinder_radius;
    const double vals[3] = {
      -wall_accel * pt.y() / radial_distance,
       wall_accel * pt.x() / radial_distance,
       0.0
    };
    const int sol_idx[3] = { node_index * 4 + 1, node_index * 4 + 2, node_index * 4 + 3 };

    VecSetValues(dot_sol->solution, 3, sol_idx, vals, INSERT_VALUES);
  }

  dot_sol->Assembly_GhostUpdate();
}

void PNonlinear_NS_Solver::update_uniform_inflow_value(
    const double &stime,
    const ALocal_InflowBC * const &infbc,
    PDNSolution * const &sol ) const
{
  const double ramped_freestream_speed = get_ramped_freestream_speed(stime);
  const int num_nbc = infbc->get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const Vector_3 outvec = infbc->get_outvec(nbc_id);
    const double vals[3] = {
      -ramped_freestream_speed * outvec.x(),
      -ramped_freestream_speed * outvec.y(),
      -ramped_freestream_speed * outvec.z()
    };

    const int numnode = infbc->get_Num_LD(nbc_id);
    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc->get_LDN(nbc_id, ii);
      const int sol_idx[3] = { node_index * 4 + 1, node_index * 4 + 2, node_index * 4 + 3 };
      VecSetValues(sol->solution, 3, sol_idx, vals, INSERT_VALUES);
    }
  }

  sol->Assembly_GhostUpdate();
}

void PNonlinear_NS_Solver::update_uniform_inflow_dot_value(
    const double &stime,
    const ALocal_InflowBC * const &infbc,
    PDNSolution * const &dot_sol ) const
{
  const double ramped_freestream_accel = get_ramped_freestream_accel(stime);
  const int num_nbc = infbc->get_num_nbc();

  for(int nbc_id=0; nbc_id<num_nbc; ++nbc_id)
  {
    const Vector_3 outvec = infbc->get_outvec(nbc_id);
    const double vals[3] = {
      -ramped_freestream_accel * outvec.x(),
      -ramped_freestream_accel * outvec.y(),
      -ramped_freestream_accel * outvec.z()
    };

    const int numnode = infbc->get_Num_LD(nbc_id);
    for(int ii=0; ii<numnode; ++ii)
    {
      const int node_index = infbc->get_LDN(nbc_id, ii);
      const int sol_idx[3] = { node_index * 4 + 1, node_index * 4 + 2, node_index * 4 + 3 };
      VecSetValues(dot_sol->solution, 3, sol_idx, vals, INSERT_VALUES);
    }
  }

  dot_sol->Assembly_GhostUpdate();
}

// EOF
