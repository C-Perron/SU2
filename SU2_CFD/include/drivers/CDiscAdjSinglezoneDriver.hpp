/*!
 * \file CDiscAdjSinglezoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 8.4.0 "Harrier"
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

#pragma once
#include "CSinglezoneDriver.hpp"
#include "../../../Common/include/linear_algebra/CPreconditioner.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"

/*!
 * \class CDiscAdjSinglezoneDriver
 * \ingroup DiscAdj
 * \brief Class for driving single-zone adjoint solvers.
 * \author R. Sanchez
 * \version 8.4.0 "Harrier"
 */
class CDiscAdjSinglezoneDriver : public CSinglezoneDriver {
protected:

protected:
#ifdef CODI_FORWARD_TYPE
  using Scalar = su2double;
#else
  using Scalar = passivedouble;
#endif

  class AdjointProductWrapper : public CMatrixVectorProduct<Scalar> {
  public:
    CDiscAdjSinglezoneDriver* const driver;
    mutable unsigned long iInnerIter = 0;

    AdjointProductWrapper(CDiscAdjSinglezoneDriver* d) : driver(d) {}
    inline void operator()(const CSysVector<Scalar> & u, CSysVector<Scalar> & v) const override {
      driver->SetAllSolutions(ZONE_0, true, u);
      driver->Iterate(iInnerIter, true);
      driver->GetAllSolutions(ZONE_0, true, v);
      v -= u;
      ++iInnerIter;
    }
  };

  class IdentityPreconditioner : public CPreconditioner<Scalar> {
  public:
    inline bool IsIdentity() const override { return true; }
    inline void operator()(const CSysVector<Scalar> & u, CSysVector<Scalar> & v) const override { v = u; }
  };

  unsigned long nAdjoint_Iter;                  /*!< \brief The number of adjoint iterations that are run on the fixed-point solver.*/
  RECORDING RecordingState;                     /*!< \brief The kind of recording the tape currently holds.*/
  RECORDING MainVariables;                      /*!< \brief The kind of recording linked to the main variables of the problem.*/
  RECORDING SecondaryVariables;                 /*!< \brief The kind of recording linked to the secondary variables of the problem.*/
  int MainSolver;                               /*!< \brief Index of the main adjoint solver. */
  su2double ObjFunc;                            /*!< \brief The value of the objective function.*/
  CIteration* direct_iteration;                 /*!< \brief A pointer to the direct iteration.*/

  CConfig *config;                              /*!< \brief Definition of the particular problem. */
  CIteration *iteration;                        /*!< \brief Container vector with all the iteration methods. */
  CIntegration **integration;                   /*!< \brief Container vector with all the integration methods. */
  CGeometry *geometry;                          /*!< \brief Geometrical definition of the problem. */
  CSolver **solver;                             /*!< \brief Container vector with all the solutions. */
  COutput *direct_output;
  CNumerics ***numerics;                        /*!< \brief Container vector with all the numerics. */

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void SetRecording(RECORDING kind_recording);

  /*!
   * \brief Run one iteration of the solver.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void DirectRun(RECORDING kind_recording);

  /*!
   * \brief Set the objective function.
   */
  void SetObjFunction(void);

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdjObjFunction(void);

  /*!
   * \brief Record the main computational path.
   */
  void MainRecording(void);

  /*!
   * \brief Record the secondary computational path.
   */
  void SecondaryRecording(void);

  /*!
   * \brief gets Convergence on physical time scale, (deactivated in adjoint case)
   * \return false
   */
  inline bool GetTimeConvergence() const override { return false; }

  /*!
   * \brief Get RHS vector for the discrete adjoint problem
   */
  bool GetAdjointRHS(CSysVector<Scalar>& rhs);

  /*!
   * \brief Run driver with a Krylov solver
   */
  void RunKrylov(void);

  /*!
   * \brief Run driver with a fixed-point method
   */
  void RunFixedPoint(void);

  /*!
   * \brief Iterate the discrete adjoint problem
   * \param[in] iInnerIter - Current inner iteration
   * \param[in] KrylovMode - Whether this is called from within a Krylov solver
   */
  bool Iterate(unsigned long iInnerIter, bool KrylovMode);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Total number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjSinglezoneDriver(void) override;

  /*!
   * \brief Preprocess the single-zone iteration
   * \param[in] TimeIter - index of the current time-step.
   */
  void Preprocess(unsigned long TimeIter) override;

  /*!
   * \brief Run a single iteration of the discrete adjoint solver with a single zone.
   */
  void Run(void) override;

  /*!
   * \brief Postprocess the adjoint iteration for ZONE_0.
   */
  void Postprocess(void) override;

  ////////////////////////////////////////////////////////////////////////////////
  /* Functions added for more granular control (START) */
  ////////////////////////////////////////////////////////////////////////////////

  /*!
   * \brief TODO
   */
  inline COutput* GetDirectOutput(unsigned short iZone) override {
    if (iZone >= nZone)
      SU2_MPI::Error("Zone " + to_string(iZone) + " is out of bound", CURRENT_FUNCTION);
    return direct_output;
  }

};
