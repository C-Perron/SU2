/*!
 * \file CNewtonIntegration.cpp
 * \brief Newton-Krylov integration.
 * \author P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CNewtonIntegration.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"

using Scalar = CNewtonIntegration::Scalar;

namespace {

class CMatrixFreeProductWrapper final : public CMatrixVectorProduct<Scalar> {
  CNewtonIntegration* integration;
public:
  CMatrixFreeProductWrapper(CNewtonIntegration* i) : integration(i) {}

  /*!
   * \brief Operator for the product operation.
   */
  inline void operator()(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const override {
    integration->MatrixFreeProduct(u, v);
  }
};

class CPreconditionerWrapper final : public CPreconditioner<Scalar> {
  const CNewtonIntegration* integration;
public:
  CPreconditionerWrapper(const CNewtonIntegration* i) : integration(i) {}

  /*!
   * \brief Operator for the preconditioning operation.
   */
  inline void operator()(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const override {
    integration->Preconditioner(u, v);
  }
};
}

CNewtonIntegration::~CNewtonIntegration() { delete preconditioner; }

void CNewtonIntegration::Setup() {
  SU2_ZONE_SCOPED

  auto iparam = config->GetNewtonKrylovIntParam();
  auto dparam = config->GetNewtonKrylovDblParam();

  startupIters = iter = iparam[0];
  startupResidual = dparam[0];
  useDeflation = iparam[3] > 0;
  precondIters = iparam[1];
  precondTol = dparam[1];
  tolRelaxFactor = iparam[2];
  fullTolResidual = dparam[2];
  finDiffStepND = SU2_TYPE::GetValue(dparam[3]);
  nkRelaxation = fmin(SU2_TYPE::GetValue(dparam[4]), 1);
  if (nkRelaxation < 0) {
    autoRelaxation = true;
    nkRelaxation = 1;
  }

  const auto nVar = solvers[FLOW_SOL]->GetnVar();
  const auto nPoint = geometry->GetnPoint();
  const auto nPointDomain = geometry->GetnPointDomain();

  omp_chunk_size = computeStaticChunkSize(nPoint, omp_get_max_threads(), 1024);

  LinSolver.SetxIsZero(true);

  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  if (autoRelaxation || nkRelaxation < 1) {
    LinSysResRelax.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }

  if (!std::is_same<Scalar,su2double>::value) {
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, nullptr);
  }

  /*--- Check if the solver is able to provide a linear preconditioner. ---*/
  if (config->GetKind_TimeIntScheme() != EULER_IMPLICIT) return;

  const auto kindPrec = static_cast<ENUM_LINEAR_SOLVER_PREC>(config->GetKind_Linear_Solver_Prec());

  preconditioner = CPreconditioner<MixedScalar>::Create(kindPrec, solvers[FLOW_SOL]->Jacobian, geometry, config);

  if (!std::is_same<Scalar,MixedScalar>::value) {
    precondIn.Initialize(nPoint, nPointDomain, nVar, nullptr);
    precondOut.Initialize(nPoint, nPointDomain, nVar, nullptr);
  }

  /*--- Only possible with a preconditioner. ---*/
  startupPeriod = (startupIters > 0) || (startupResidual < 0.0);

}

void CNewtonIntegration::PerturbSolution(const CSysVector<Scalar>& dir, Scalar mag) {
  SU2_ZONE_SCOPED

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint) {
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < solvers[FLOW_SOL]->GetnVar(); ++iVar)
      solvers[FLOW_SOL]->GetNodes()->AddSolution(iPoint,iVar, mag*dir(iPoint,iVar));
  }
  END_SU2_OMP_FOR
}

void CNewtonIntegration::ComputeResiduals(ResEvalType type) {
  SU2_ZONE_SCOPED

  /*--- Save the default integration scheme, and force to explicit if required. ---*/
  auto TimeIntScheme = config->GetKind_TimeIntScheme();
  if (type == ResEvalType::EXPLICIT) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetKind_TimeIntScheme(EULER_EXPLICIT);)
  }

  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  if (type == ResEvalType::DEFAULT) {
    solvers[FLOW_SOL]->SetTime_Step(geometry, solvers, config, MESH_0, config->GetTimeIter());
  }

  Space_Integration(geometry, solvers, numerics[FLOW_SOL], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS);

  /*--- Restore default. ---*/
  if (type == ResEvalType::EXPLICIT) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetKind_TimeIntScheme(TimeIntScheme);)
  }

}

void CNewtonIntegration::ComputeFinDiffStep() {
  SU2_ZONE_SCOPED

  static su2double rmsSol;
  su2double rmsSol_loc = 0.0;

  SU2_OMP_MASTER
  rmsSol = 0.0;
  END_SU2_OMP_MASTER

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint)
    for (auto iVar = 0ul; iVar < solvers[FLOW_SOL]->GetnVar(); ++iVar)
      rmsSol_loc += pow(solvers[FLOW_SOL]->GetNodes()->GetSolution(iPoint,iVar), 2);
  END_SU2_OMP_FOR

  atomicAdd(rmsSol_loc, rmsSol);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    su2double t = rmsSol;
    SU2_MPI::Allreduce(&t, &rmsSol, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    finDiffStep = finDiffStepND * max(1.0, sqrt(SU2_TYPE::GetValue(rmsSol) / geometry->GetGlobal_nPointDomain()));
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CNewtonIntegration::MultiGrid_Iteration(CGeometry ****geometry_, CSolver *****solvers_, CNumerics ******numerics_,
                                             CConfig **config_, unsigned short EqSystem, unsigned short iZone,
                                             unsigned short iInst) {
  SU2_ZONE_SCOPED

  config = config_[iZone];
  solvers = solvers_[iZone][iInst][MESH_0];
  geometry = geometry_[iZone][iInst][MESH_0];
  numerics = numerics_[iZone][iInst][MESH_0];

  if (!setup) { Setup(); setup = true; }

  /*--- Remove NK relaxation to compute the current residual. ---*/
  config->SetNewtonKrylovRelaxation(1.0);

  /*--- When using NK relaxation (not fully 2nd order Jacobian products) we need an additional
   * residual evaluation that is used as the reference for finite differences. ---*/
  LinSysRes0 = (!startupPeriod && nkRelaxation < 1) ? &LinSysResRelax : &LinSysRes;

  SU2_OMP_PARALLEL_(if(solvers[FLOW_SOL]->GetHasHybridParallel())) {

  /*--- Save the current solution to be able to perturb it. ---*/

  solvers[FLOW_SOL]->Set_OldSolution();

  /*--- Current residual. ---*/

  ComputeResiduals(ResEvalType::DEFAULT);

  /*--- Compute the approximate Jacobian for preconditioning. ---*/

  solvers[FLOW_SOL]->PrepareImplicitIteration(geometry, solvers, config);

  if (preconditioner) preconditioner->Build();

  auto CopyLinSysRes = [&](int sign, auto& dst) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto i = 0ul; i < dst.GetNElmDomain(); ++i)
      dst[i] = sign * SU2_TYPE::GetValue(solvers[FLOW_SOL]->LinSysRes[i]);
    END_SU2_OMP_FOR
  };
  CopyLinSysRes(1, LinSysRes);

  if (!startupPeriod && nkRelaxation < 1) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetNewtonKrylovRelaxation(nkRelaxation);)
    ComputeResiduals(ResEvalType::EXPLICIT);
    // Here the sign was not flipped by PrepareImplicitIteration.
    CopyLinSysRes(-1, LinSysResRelax);
  }

  su2double residual = 0.0;
  for (auto iVar = 0ul; iVar < LinSysRes.GetNVar(); ++iVar)
    residual += log10(solvers[FLOW_SOL]->GetRes_RMS(iVar)) / LinSysRes.GetNVar();

  /*--- Check if startup period should end after this iteration. ---*/

  bool endStartup = false;

  if (startupPeriod) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      firstResidual = max(firstResidual, residual);
      if (iter) iter -= 1;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
    endStartup = (iter == 0) && (residual - firstResidual < startupResidual);
  }

  /*--- The NK solves are expensive, the tolerance is relaxed while the residuals are high. ---*/

  Scalar toleranceFactor = 1.0;

  if (!startupPeriod && tolRelaxFactor > 1 && fullTolResidual < 0.0) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(firstResidual = max(firstResidual, residual);)
    su2double x = (residual - firstResidual) / fullTolResidual;
    toleranceFactor = 1.0 + (tolRelaxFactor-1)*max(0.0, 1.0-SU2_TYPE::GetValue(x));
  }

  /*--- Solve for the solution update. ---*/

  auto iter = config->GetLinear_Solver_Iter();
  Scalar eps = SU2_TYPE::GetValue(config->GetLinear_Solver_Error());

  auto& linSysSol = GetSolutionVec(solvers[FLOW_SOL]->LinSysSol);

  if (startupPeriod) {
    iter = Preconditioner_impl(LinSysRes, linSysSol, iter, eps);
  } else {
    ComputeFinDiffStep();

    eps *= toleranceFactor;
    if (useDeflation) {
      iter = LinSolver.FGCRODR_LinSolver(LinSysRes, linSysSol, CMatrixFreeProductWrapper(this),
                                         CPreconditionerWrapper(this), eps, iter, eps, false, config);
    } else {
      iter = LinSolver.FGMRES_LinSolver(LinSysRes, linSysSol, CMatrixFreeProductWrapper(this),
                                        CPreconditionerWrapper(this), eps, iter, eps, false, config);
    }
    /*--- Scale back the residual to trick the CFL adaptation. ---*/
    eps /= toleranceFactor;
  }
  SetSolutionResult(solvers[FLOW_SOL]->LinSysSol);

  unsigned short lineSearchIter = 10;

  if (lineSearchIter > 0 && !startupPeriod) {
    auto lineSearchSuccess = ExecuteLineSearch(linSysSol, lineSearchIter);
    // Trigger a CFL adaptation if the line search fails
    if (!lineSearchSuccess) eps = 1.0;
    SetSolutionResult(solvers[FLOW_SOL]->LinSysSol);
  }

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    solvers[FLOW_SOL]->SetIterLinSolver(iter);
    solvers[FLOW_SOL]->SetResLinSolver(eps);

    if (!startupPeriod && autoRelaxation) {
      const su2double adaptTol = config->GetCFL_Adapt() ? config->GetCFL_AdaptParam(4) : 0;
      if (eps > fmax(config->GetLinear_Solver_Error(), adaptTol)) {
        nkRelaxation *= 0.9;
      } else if (eps < 0.9 * fmax(config->GetLinear_Solver_Error(), adaptTol)) {
        nkRelaxation = fmin(fmax(nkRelaxation * 1.05, 0.05), 1);
      }
    }
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /// TODO: Clever back-tracking and CFL adaptation based on residual reduction.

  /// WIP: Line search


  /*--- Update solution. ---*/

  solvers[FLOW_SOL]->CompleteImplicitIteration(geometry, solvers, config);

  /*--- Call the various post processings. ---*/

  solvers[FLOW_SOL]->Preprocessing(geometry, solvers, config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  solvers[FLOW_SOL]->Postprocessing(geometry, solvers, config, MESH_0);

  SU2_OMP_MASTER {
    solvers[FLOW_SOL]->Pressure_Forces(geometry, config);
    solvers[FLOW_SOL]->Momentum_Forces(geometry, config);
    solvers[FLOW_SOL]->Friction_Forces(geometry, config);
  }
  END_SU2_OMP_MASTER

  /*--- At the end of the startup period the CFL is reset to the initial value. ---*/

  if (endStartup) {
    SU2_OMP_MASTER {
      startupPeriod = false;
      firstResidual = residual;
      if (autoRelaxation) nkRelaxation = 0;
    }
    END_SU2_OMP_MASTER
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); ++iPoint)
      solvers[FLOW_SOL]->GetNodes()->SetLocalCFL(iPoint, config->GetCFL(MESH_0));
    END_SU2_OMP_FOR
  }

  }
  END_SU2_OMP_PARALLEL
}

void CNewtonIntegration::MatrixFreeProduct(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) {
  SU2_ZONE_SCOPED

  Scalar factor = finDiffStep / u.norm();

  PerturbSolution(u, factor);

  ComputeResiduals(ResEvalType::EXPLICIT);

  /*--- Finalize product. ---*/
  factor = 1.0 / factor;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
    su2double delta = (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint)) /
                      max(EPS, solvers[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < LinSysRes0->GetNVar(); ++iVar) {
      Scalar perturbRes = SU2_TYPE::GetValue(solvers[FLOW_SOL]->LinSysRes(iPoint,iVar));

      /*--- The global residual had its sign flipped, so we add to get the difference. ---*/
      v(iPoint,iVar) = (perturbRes + (*LinSysRes0)(iPoint,iVar)) * factor;

      /*--- Pseudotime term of the true Jacobian. ---*/
      v(iPoint,iVar) += SU2_TYPE::GetValue(delta) * u(iPoint,iVar);
    }
  }
  END_SU2_OMP_FOR

  CSysMatrixComms::Initiate(v, geometry, config);
  CSysMatrixComms::Complete(v, geometry, config);
}

void CNewtonIntegration::Preconditioner(const CSysVector<Scalar>& u, CSysVector<Scalar>& v) const {
  SU2_ZONE_SCOPED

  if (preconditioner) {
    Scalar eps = SU2_TYPE::GetValue(precondTol);
    Preconditioner_impl(u, v, precondIters, eps);
  }
  else {
    /*--- Approximate diagonal preconditioner. ---*/

    const bool dt1st = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST);
    const bool dt2nd = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
    const su2double dt = config->GetDelta_UnstTimeND() * (dt1st + 1.5 * dt2nd);

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
      su2double delta = (solvers[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint) + dt) /
                        (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint));
      SU2_OMP_SIMD
      for (auto iVar = 0ul; iVar < u.GetNVar(); ++iVar)
        v(iPoint,iVar) = SU2_TYPE::GetValue(delta) * u(iPoint,iVar);
    }
    END_SU2_OMP_FOR

    CSysMatrixComms::Initiate(v, geometry, config);
    CSysMatrixComms::Complete(v, geometry, config);
  }
}

bool CNewtonIntegration::ExecuteLineSearch(CSysVector<Scalar>& linSysSol, unsigned short maxIter) {
  unsigned short iIter;
  bool success;
  Scalar step = 1.0;
  Scalar minStep = 1e-6;

  // Compute the initial residual norm
  lineSearchRes0 = LinSysRes.norm();

  // Restrict the search direction based on the max allowable flow changes
  RestrictLineSearch(linSysSol);

  for (iIter = 0; iIter < maxIter; iIter++) {
    success = IterateLineSearch(linSysSol, step);

    cout << "LS: " << iIter << ", step: " << step << ", res: " << lineSearchRes / lineSearchRes0 << endl;

    if (success || (step <= minStep)) break;

    step = max(minStep, UpdateLineSearch(iIter, step));

  }

  // Report the result of the line search
  if (!success) {
    if (step <= minStep) {
      cout << "Line search failed. Minimum step size reached." << endl;
    }
    else {
      cout << "Line search failed. Maximum iterations reached" << endl;
    }
  }
  else {
    cout << "Line search successful." << endl;
  }

  // Skip update if the initial solution was adequate
  if (iIter > 0) linSysSol *= step;

  return success;

}


void CNewtonIntegration::RestrictLineSearch(CSysVector<Scalar>& linSysSol) {

  // Update relaxation factors for the current solution
  solvers[FLOW_SOL]->ComputeUnderRelaxationFactor(config);

  // Find minimum relaxation factor across all points
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {
    auto urf = solvers[FLOW_SOL]->GetNodes()->GetUnderRelaxation(iPoint);
    for (auto iVar = 0ul; iVar < linSysSol.GetNVar(); ++iVar)
      linSysSol(iPoint,iVar) *= urf;
  }
  END_SU2_OMP_FOR
}

bool CNewtonIntegration::IterateLineSearch(const CSysVector<Scalar>& direction, const Scalar& step) {

  PerturbSolution(direction, step);
  ComputeResiduals(ResEvalType::EXPLICIT);

  // Sufficient decrease condition (Armijo–Goldstein)
  Scalar alpha = 1e-4;

  static su2double sq_res;
  su2double sq_res_loc = 0.0;

  SU2_OMP_MASTER
  sq_res = 0.0;
  END_SU2_OMP_MASTER

  // Compute the residual at the perturbed solution
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < geometry->GetnPointDomain(); ++iPoint) {

    Scalar delta = SU2_TYPE::GetValue(
      (geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint)) /
      max(EPS, solvers[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint))
    );

    SU2_OMP_SIMD
    for (auto iVar = 0ul; iVar < LinSysRes0->GetNVar(); ++iVar) {
      Scalar res_tmp = step * delta * direction(iPoint,iVar);
      res_tmp += SU2_TYPE::GetValue(solvers[FLOW_SOL]->LinSysRes(iPoint,iVar));

      sq_res_loc += res_tmp * res_tmp;
    }
  }
  END_SU2_OMP_FOR

  atomicAdd(sq_res_loc, sq_res);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS
  {
    su2double t = sq_res;
    SU2_MPI::Allreduce(&t, &sq_res, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    lineSearchRes = sqrt(SU2_TYPE::GetValue(sq_res));
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  // Check if the sufficient decrease condition is satisfied
  return (lineSearchRes < ((1 - alpha * step) * lineSearchRes0));

}

Scalar CNewtonIntegration::UpdateLineSearch(const Scalar& iIter, const Scalar& step) {
  Scalar maxFactor = 0.5; // Maximum reduction factor
  Scalar minFactor = 0.1; // Minimum reduction factor

  if (iIter == 0) {
    stepOld = step; // Store the initial step size
    lineSearchResOld = lineSearchRes; // Store the initial residual
    return step * maxFactor; // Reduce step size by maxFactor on first iteration
  }

  // Three-point parabolic interpolation to estimate the next step size

  Scalar ds10 = step;
  Scalar dr10 = lineSearchRes - lineSearchRes0;

  Scalar ds12 = step - stepOld;
  Scalar dr12 = lineSearchRes - lineSearchResOld;  

  Scalar num = ds10 * ds10 * dr12 - ds12 * ds12 * dr10;
  Scalar den = ds10 * dr12 - ds12 * dr10;

  // Check for division by zero
  if (abs(den) < EPS) {
    return step * maxFactor;
  }

  Scalar newStep = min(max(step - 0.5 * num / den, step * minFactor), step * maxFactor);

  cout << step - 0.5 * num / den << endl;

  // Update previous step and residual
  stepOld = step;
  lineSearchResOld = lineSearchRes;

  return newStep;

}
