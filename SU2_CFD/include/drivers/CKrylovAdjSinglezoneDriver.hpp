/*!
 * \file CDiscAdjSinglezoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
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

#pragma once
#include "CSinglezoneDriver.hpp"

/*!
 * \class CKrylovAdjSinglezoneDriver
 * \ingroup DiscAdj
 * \brief Class for driving single-zone adjoint solvers.
 * \author R. Sanchez
 * \version 8.1.0 "Harrier"
 */
class CKrylovAdjSinglezoneDriver : public CSinglezoneDriver {
protected:

  enum class ADJ_OUTPUT {
    NONE,
    OUTPUTS,
    RESIDUALS,
  };

  enum class ADJ_INPUT {
    NONE,
    SOLUTION,
    VARIABLES,
    MESH_COORDS,
    MESH_DEFORM,
  };

  unsigned long nAdjoint_Iter;                  /*!< \brief The number of adjoint iterations that are run on the fixed-point solver.*/
  ADJ_OUTPUT TapeOutput;
  ADJ_INPUT TapeInput;
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

  CSysVector<passivedouble> AdjSysRes;
  CSysVector<passivedouble> AdjSysSol;

  /*!
   * \brief Record one iteration of a flow iteration in within multiple zones.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void SetRecording(ADJ_OUTPUT adj_out, ADJ_INPUT adj_inp, const bool initial = false);

  inline void ClearRecording(const bool initial = false) {
    SetRecording(ADJ_OUTPUT::NONE, ADJ_INPUT::NONE, initial);
  }

  /*!
   * \brief Run one iteration of the solver.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void DirectRun(void);

  /*!
   * \brief Set the objective function.
   */
  void SetObjFunction(void);

  inline void RegisterObjFunction(void) {if (rank == MASTER_NODE) AD::RegisterOutput(ObjFunc);}

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdjObjFunction(void);

  /*!
   * \brief Record the main computational path.
   */
  void SetMainAdjointSystem(void);

  /*!
   * \brief Record the secondary computational path.
   */
  void SecondaryRecording(void);

  /*!
   * \brief gets Convergence on physical time scale, (deactivated in adjoint case)
   * \return false
   */
  inline bool GetTimeConvergence() const override { return false; }

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Total number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CKrylovAdjSinglezoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CKrylovAdjSinglezoneDriver(void) override;

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
};
