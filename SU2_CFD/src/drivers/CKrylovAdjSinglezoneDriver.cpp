/*!
 * \file driver_adjoint_singlezone.cpp
 * \brief The main subroutines for driving adjoint single-zone problems.
 * \author R. Sanchez
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/drivers/CKrylovAdjSinglezoneDriver.hpp"
#include "../../include/output/tools/CWindowingTools.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CTurboIteration.hpp"
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"

CKrylovAdjSinglezoneDriver::CKrylovAdjSinglezoneDriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                                                 val_nZone,
                                                                                                 MPICommunicator) {


  /*--- Store the number of internal iterations that will be run by the adjoint solver ---*/
  nAdjoint_Iter = config_container[ZONE_0]->GetnInner_Iter();

  /*--- Store the pointers ---*/
  config      = config_container[ZONE_0];
  iteration   = iteration_container[ZONE_0][INST_0];
  solver      = solver_container[ZONE_0][INST_0][MESH_0];
  numerics    = numerics_container[ZONE_0][INST_0][MESH_0];
  geometry    = geometry_container[ZONE_0][INST_0][MESH_0];
  integration = integration_container[ZONE_0][INST_0];

  /*--- Store the recording state ---*/
  TapeOutput = ADJ_OUTPUT::NONE;
  TapeInput = ADJ_INPUT::NONE;

  /*--- Initialize the direct iteration ---*/

  // TODO: Remove error after implementing more solvers
  if (config->GetKind_Solver() != MAIN_SOLVER::DISC_ADJ_EULER)
    SU2_MPI::Error("Solver Kind is not implemented!", CURRENT_FUNCTION);

  switch (config->GetKind_Solver()) {

  case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
  case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;

    if (config->GetBoolTurbomachinery()) {
      direct_iteration = new CTurboIteration(config);
    }
    else { direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::EULER, config); }

    if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
      direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::EULER, config, nDim);
    }
    else { direct_output =  COutputFactory::CreateOutput(MAIN_SOLVER::INC_EULER, config, nDim); }

    MainVariables = RECORDING::SOLUTION_VARIABLES;
    if (config->GetDeform_Mesh()) {
      SecondaryVariables = RECORDING::MESH_DEFORM;
    }
    else { SecondaryVariables = RECORDING::MESH_COORDS; }
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM_EULER : case MAIN_SOLVER::DISC_ADJ_FEM_NS : case MAIN_SOLVER::DISC_ADJ_FEM_RANS :
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;
    direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_EULER, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::FEM_EULER, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: elasticity equation." << endl;
    direct_iteration =  CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_ELASTICITY, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::FEM_ELASTICITY, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFEA_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_HEAT:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: heat equation." << endl;
    direct_iteration = CIterationFactory::CreateIteration(MAIN_SOLVER::HEAT_EQUATION, config);
    direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::HEAT_EQUATION, config, nDim);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJHEAT_SOL;
    break;

  default:
    break;

  }

  /*--- Disable solution updates in all direct solvers since we only need the residuals ---*/

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    for (unsigned short iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      auto solver = solver_container[ZONE_0][INST_0][iMesh][iSol];
      if (solver && !solver->GetAdjoint()) {
        solver->SkipSolutionUpdate = true;
      }
    }
  }

  /*--- Disable multi-grid since we don't update any direct solutions ---*/

  config->SetMGLevels(0);

  // TODO: Temporary fix for solution adjoint extraction

  config->SetRelaxation_Factor_Adjoint(1.0);

  /*--- Initialize RHS vector ---*/

  const auto nVar = GetTotalNumberOfVariables(ZONE_0, true);
  const auto nPoint = geometry->GetnPoint();
  const auto nPointDomain = geometry->GetnPointDomain();

  AdjSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  AdjSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Preprocess history output ---*/

  direct_output->PreprocessHistoryOutput(config, false);

}

CKrylovAdjSinglezoneDriver::~CKrylovAdjSinglezoneDriver() {

  delete direct_iteration;
  delete direct_output;

}

void CKrylovAdjSinglezoneDriver::Preprocess(unsigned long TimeIter) {

  /*--- Set the current time iteration in the config and also in the driver
   * because the python interface doesn't offer an explicit way of doing it. ---*/

  this->TimeIter = TimeIter;
  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- Preprocess the adjoint iteration ---*/

  iteration->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                        solver_container, numerics_container, config_container,
                        surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Clear previous record and check direct solution residuals ---*/
  ClearRecording(true);

  /*--- Prepare the adjoint linear system ---*/
  SetMainAdjointSystem();
}

void CKrylovAdjSinglezoneDriver::Run() {

  SU2_MPI::Error("Run function not implemeted!", CURRENT_FUNCTION);

  CQuasiNewtonInvLeastSquares<passivedouble> fixPtCorrector;
  if (config->GetnQuasiNewtonSamples() > 1) {
    fixPtCorrector.resize(config->GetnQuasiNewtonSamples(),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(),
                          GetTotalNumberOfVariables(ZONE_0,true),
                          geometry_container[ZONE_0][INST_0][MESH_0]->GetnPointDomain());

    if (TimeIter != 0) GetAllSolutions(ZONE_0, true, fixPtCorrector);
  }

  for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {

    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *--- of the previous iteration. The values are passed to the AD tool.
     *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

    config->SetInnerIter(Adjoint_Iter);

    iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdjObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/

    iteration->IterateDiscAdj(geometry_container, solver_container,
                              config_container, ZONE_0, INST_0, false);

    /*--- Monitor the pseudo-time ---*/

    StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                  solver_container, numerics_container, config_container,
                                  surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    /*--- Output files for steady state simulations. ---*/

    if (!config->GetTime_Domain()) {
      iteration->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, INST_0);
    }

    if (StopCalc) break;

    /*--- Correct the solution with the quasi-Newton approach. ---*/

    if (fixPtCorrector.size()) {
      GetAllSolutions(ZONE_0, true, fixPtCorrector.FPresult());
      SetAllSolutions(ZONE_0, true, fixPtCorrector.compute());
    }

  }

}

void CKrylovAdjSinglezoneDriver::Postprocess() {

  switch(config->GetKind_Solver())
  {
    case MAIN_SOLVER::DISC_ADJ_EULER :     case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES :     case MAIN_SOLVER::DISC_ADJ_RANS :
    case MAIN_SOLVER::DISC_ADJ_INC_EULER : case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES : case MAIN_SOLVER::DISC_ADJ_INC_RANS :
    case MAIN_SOLVER::DISC_ADJ_HEAT :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();
      break;

    case MAIN_SOLVER::DISC_ADJ_FEM :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();

      iteration->Postprocess(output_container[ZONE_0], integration_container, geometry_container,
                             solver_container, numerics_container, config_container,
                             surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
      break;

    default:
      break;

  }//switch

}

void CKrylovAdjSinglezoneDriver::SetRecording(ADJ_OUTPUT adj_out, ADJ_INPUT adj_inp, const bool initial){

  AD::Reset();

  /*--- If this is the initial pass, let solver use their existing time integration scheme so
        that Jacobians can be formed. Use EULER_EXPLICIT afterward to avoid reforming Jacobians ---*/

  if (initial) config->SetForceExplicitTimeInt(true);

  /*--- Store the recording state ---*/

  RECORDING kind_recording;
  const bool clear_tape = (adj_out != ADJ_OUTPUT::NONE) && (adj_inp != ADJ_INPUT::NONE);

  if (clear_tape) {
    TapeOutput = ADJ_OUTPUT::NONE;
    TapeInput = ADJ_INPUT::NONE;
    kind_recording = RECORDING::CLEAR_INDICES;
  }
  else {
    TapeOutput = adj_out;
    TapeInput = adj_inp;

    switch (TapeInput)
    {

    case ADJ_INPUT::SOLUTION:
    case ADJ_INPUT::VARIABLES:
      kind_recording = RECORDING::SOLUTION_VARIABLES;
      break;
    case ADJ_INPUT::MESH_COORDS:
      kind_recording = RECORDING::MESH_COORDS;
      break;

    case ADJ_INPUT::MESH_DEFORM:
      kind_recording = RECORDING::MESH_DEFORM;
      break;

    default:
      kind_recording = RECORDING::CLEAR_INDICES;
      break;
    }
  }

  /*--- Prepare for recording by resetting the solution to the initial converged solution. ---*/

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    for (unsigned short iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      auto solver = solver_container[ZONE_0][INST_0][iMesh][iSol];
      if (solver && solver->GetAdjoint()) {
        SU2_OMP_PARALLEL_(if (solver->GetHasHybridParallel()))
        solver->SetRecording(geometry_container[ZONE_0][INST_0][iMesh], config_container[ZONE_0]);
        END_SU2_OMP_PARALLEL
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "\n-------------------------------------------------------------------------\n";
    if (clear_tape) cout << "Clearing the computational graph." << endl;

    switch (TapeOutput) {
      case ADJ_OUTPUT::OUTPUTS: cout << "Recording computational graph of SOLVER OUTPUTS w.r.t. "; break;
      case ADJ_OUTPUT::RESIDUALS: cout << "Recording computational graph of SOLVER RESIDUALS w.r.t. "; break;
      default: break;
    }

    switch (TapeInput) {
      case ADJ_INPUT::SOLUTION: cout << "SOLVER SOLUTION." << endl; break;
      case ADJ_INPUT::VARIABLES: cout << "SOLVER VARIABLES." << endl; break;
      case ADJ_INPUT::MESH_COORDS: cout << "MESH COORDINATES." << endl; break;
      case ADJ_INPUT::MESH_DEFORM: cout << "MESH DEFORMATION." << endl; break;
      default: break;
    }

    if (initial) cout << "Computing residuals to check the convergence of the direct problem." << endl;
  }

  if (!clear_tape) AD::StartRecording();

  /*---Enable recording and register input of the iteration --- */

  if (TapeInput == ADJ_INPUT::SOLUTION){
    iteration->RegisterInputSolution(solver_container, geometry_container, config_container, ZONE_0, INST_0);
  }

  /*--- Set the dependencies of the iteration ---*/

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                             INST_0, kind_recording);

  /*--- Do one iteration of the direct solver ---*/

  DirectRun();

  /*--- Print the direct residual to screen ---*/
  // TODO: Change flag to not always print the residuals
  if (initial) (RECORDING::SOLUTION_VARIABLES);

  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  /*--- Register tape outputs --- */

  switch (TapeOutput)
  {
  case ADJ_OUTPUT::OUTPUTS:
    iteration->RegisterOutputVariables(solver_container, geometry_container, config_container, ZONE_0, INST_0);
    RegisterObjFunction();
    break;

  case ADJ_OUTPUT::RESIDUALS:
    iteration->RegisterOutputResiduals(solver_container, config_container, ZONE_0, INST_0);
    break;

  default:
    break;
  }

  if (!clear_tape && config->GetWrt_AD_Statistics()) {
    AD::PrintStatistics(SU2_MPI::GetComm(), rank == MASTER_NODE);
  }

  AD::StopRecording();

  /*--- Restore default time integration scheme ---*/

  config->SetForceExplicitTimeInt(false);

}

void CKrylovAdjSinglezoneDriver::SetAdjObjFunction(){

  if (TapeOutput != ADJ_OUTPUT::OUTPUTS) return;

  su2double seeding = 1.0;

  if (config->GetTime_Domain()) {
    const auto IterAvg_Obj = config->GetIter_Avg_Objective();
    if (TimeIter < IterAvg_Obj) {
      /*--- Default behavior when no window is chosen is to use Square-Windowing, i.e. the numerator equals 1.0 ---*/
      auto windowEvaluator = CWindowingTools();
      const su2double weight = windowEvaluator.GetWndWeight(config->GetKindWindow(), TimeIter, IterAvg_Obj - 1);
      seeding = weight / IterAvg_Obj;
    }
    else {
      seeding = 0.0;
    }
  }
  if (rank == MASTER_NODE) {
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }
}

void CKrylovAdjSinglezoneDriver::SetObjFunction(){

  ObjFunc = 0.0;

  /*--- Specific scalar objective functions ---*/

  switch (config->GetKind_Solver()) {
  case MAIN_SOLVER::DISC_ADJ_INC_EULER:       case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:      case MAIN_SOLVER::DISC_ADJ_INC_RANS:
  case MAIN_SOLVER::DISC_ADJ_EULER:           case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:          case MAIN_SOLVER::DISC_ADJ_RANS:
  case MAIN_SOLVER::DISC_ADJ_FEM_EULER:       case MAIN_SOLVER::DISC_ADJ_FEM_NS:                 case MAIN_SOLVER::DISC_ADJ_FEM_RANS:

    /*--- Surface based obj. function ---*/

    direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                     config->GetOuterIter(), config->GetInnerIter());
    ObjFunc += solver[FLOW_SOL]->GetTotal_ComboObj();
    break;

  case MAIN_SOLVER::DISC_ADJ_HEAT:
    direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                     config->GetOuterIter(), config->GetInnerIter());
    ObjFunc = solver[HEAT_SOL]->GetTotal_ComboObj();
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM:
    solver[FEA_SOL]->Postprocessing(geometry, config, numerics_container[ZONE_0][INST_0][MESH_0][FEA_SOL], true);

    direct_output->SetHistoryOutput(geometry, solver, config, config->GetTimeIter(),
                                   config->GetOuterIter(), config->GetInnerIter());
    ObjFunc = solver[FEA_SOL]->GetTotal_ComboObj();
    break;

  default:
    break;
  }
}

void CKrylovAdjSinglezoneDriver::DirectRun(){

  /*--- Mesh movement ---*/
  // TODO: Make this work with mesh deformation
  // direct_iteration->SetMesh_Deformation(geometry_container[ZONE_0][INST_0], solver, numerics, config, kind_inp_rec);

  /*--- Zone preprocessing ---*/

  direct_iteration->Preprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Iterate the direct solver ---*/

  direct_iteration->Iterate(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Postprocess the direct solver ---*/

  direct_iteration->Postprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
}

void CKrylovAdjSinglezoneDriver::SetMainAdjointSystem(){

  if ((TapeOutput != ADJ_OUTPUT::NONE) || (TapeInput != ADJ_INPUT::NONE)) ClearRecording();

  /*--- Compute and store RHS of the adjoint problem ---*/

  // Record computational graph of the objective w.r.t. the direct solution
  SetRecording(ADJ_OUTPUT::OUTPUTS, ADJ_INPUT::SOLUTION);
  // Initialize the adjoint of the objective
  SetAdjObjFunction();
  // Extract the derivatives w.r.t. the direct solution
  AD::ComputeAdjoint();

  iteration->ExtractAdjointSolutionVector(AdjSysRes, geometry_container, solver_container,
      config_container, ZONE_0, INST_0);

  // Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();

  /*--- Setup recording for solving the linear system ---*/

  // Clear previous record
  ClearRecording();
  // Record computational graph of the residulas w.r.t. the direct solution
  SetRecording(ADJ_OUTPUT::RESIDUALS, ADJ_INPUT::SOLUTION);

}

void CKrylovAdjSinglezoneDriver::SecondaryRecording(){
  SU2_MPI::Error("Not implemented!", CURRENT_FUNCTION);
  // /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with
  //  *    RECORDING::CLEAR_INDICES as argument ensures that all information from a previous recording is removed. ---*/

  // SetRecording(RECORDING::CLEAR_INDICES);

  // /*--- Store the computational graph of one direct iteration with the secondary variables as input. ---*/

  // SetRecording(SecondaryVariables);

  // /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
  //  *    of the current iteration. The values are passed to the AD tool. ---*/

  // iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  // /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  // SetAdjObjFunction();

  // /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  // AD::ComputeAdjoint();

  // /*--- Extract the computed sensitivity values. ---*/

  // if (SecondaryVariables == RECORDING::MESH_COORDS) {
  //   solver[MainSolver]->SetSensitivity(geometry, config);
  // }
  // else { // MESH_DEFORM
  //   solver[ADJMESH_SOL]->SetSensitivity(geometry, config, solver[MainSolver]);
  // }

  // /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  // AD::ClearAdjoints();

}
