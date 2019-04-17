#ifndef SU2_SOLVER_REACTIVE
#define SU2_SOLVER_REACTIVE

#include "solver_structure.hpp"
#include "variable_reactive.hpp"

/*! \class CReactiveEulerSolver
 *  \brief Main class for defining a solver for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerSolver: public CSolver {
public:
  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveNSVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:
  static LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  unsigned short nSpecies; /*!< \brief Total number of species. */
  unsigned short nPrimVarLim; /*!< \brief Number of primitive variables to limit. */
  unsigned short nMarker;  /*!< \brief Number of markers for deallocation. */
  std::vector<unsigned long> nVertex; /*!< \brief Store the number of vertices on each marker for deallocation. */

  bool  space_centered,  /*!< \brief True if space centered scheme used. */
        implicit,      /*!< \brief True if euler implicit scheme used. */
        grid_movement, /*!< \brief True if grid movement is used. */
        least_squares,  /*!< \brief True if computing gradients by least squares. */
        second_order, /*!< \brief True if second order recosntruction is applied. */
        limiter; /*!< \brief True if limiting strategy is applied. */

  bool  US_System;  /*!< \brief True if using US units. */

  su2double*** CharacPrimVar;  /*!< \brief Value of the characteristic variables at each boundary. */

  RealVec   Lower_Limit,   /*!< \brief Lower limit conserved variables. */
            Upper_Limit;   /*!< \brief Upper limit conserved variables. */

  su2double Density_Inf,       /*!< \brief Free stream density. */
            Pressure_Inf,		  /*!< \brief Free stream pressure. */
	          Temperature_Inf;   /*!< \brief Trans.-rot. free stream temperature. */

  RealVec   Velocity_Inf,  /*!< \brief Free stream flow velocity. */
            MassFrac_Inf;  /*!< \brief Free stream species mass fraction. */

  RealVec   PrimVar_i,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point i. */
            PrimVar_j,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point j. */
            PrimVar_Vertex;  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at boundary node. */

  RealVec   Prim_i, /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point i. */
            Prim_j, /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point j. */
            Primitive; /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at boundary node. */

  RealVec   Primitive_i, Primitive_j, /*!< \brief Auxiliary nPrimVar vector for storing primitive at point i and j for limiting. */
            Secondary_i, Secondary_j;/*!< \brief Auxiliary nVar vectors for storing pressure derivatives at point i and j for limiting. */

  RealVec   Buffer_Receive_U, /*!< \brief Auxiliary vector to receive information in case of parallel simulation. */
            Buffer_Send_U;    /*!< \brief Auxiliary vector to send information in case of parallel simulation. */

  RealVec Ys_i, Ys_j;  /*!< \brief Auxiliary vectors to store mass fractions at node i and j. */
  RealVec Ys;         /*!< \brief Auxiliary vector to store mass fractions. */

protected:
  unsigned short T_INDEX_PRIM, VX_INDEX_PRIM,
                 P_INDEX_PRIM, RHO_INDEX_PRIM,
                 H_INDEX_PRIM, A_INDEX_PRIM,
                 RHOS_INDEX_PRIM;               /*!< \brief Mapping for position in primitives array. */

  unsigned short RHO_INDEX_SOL, RHOVX_INDEX_SOL,
                 RHOE_INDEX_SOL, RHOS_INDEX_SOL; /*!< \brief Mapping for position in conserved array. */

  unsigned short T_INDEX_GRAD, VX_INDEX_GRAD,
                 P_INDEX_GRAD;                 /*!< \brief Mapping for position in primitives gradient. */

  unsigned short T_INDEX_LIM, VX_INDEX_LIM,
                 P_INDEX_LIM;                 /*!< \brief Mapping for position for limited variables. */

public:
  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerSolver();

	/*!
	 * \overloaded constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerSolver();

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

  /*!
   * \brief Get the pointer to the external library for physical-chemical properties
   */
  inline static LibraryPtr GetLibrary(void) {
    return library;
  }

  /*!
   * \brief Value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the variables at the boundaries.
   */
  inline su2double* GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) override {
    return CharacPrimVar[val_marker][val_vertex];
  }

  /*!
   * \brief Set primitive variables in each point reporting non physical data
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned weighted least squares.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Call MPI to set solution in case of parallel simulation
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Solution(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Call MPI to set limiter of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Call MPI to set gradient of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Gradient(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                    unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep,
                             unsigned short iMesh, unsigned short RunTime_EqSystem) override;

  /*!
   * \brief Compute the spatial integration static const unsigned a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Compute the spatial integration static const unsigned a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                       CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Source term integration.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] second_numerics - Description of the second numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iMesh - Index of the mesh in multigrid computations.
    */
   void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CNumerics* second_numerics,
                        CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Set the free-stream solution all over the domain.
    * \param[in] config - Definition of the particular problem.
    */
   void SetFreeStream_Solution(CConfig* config) override;

   /*!
    * \brief Impose via the residual the Euler wall boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                      CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the far-field boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                     CConfig* config, unsigned short val_marker) override {}

   /*!
    * \brief Impose the inlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                 CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose a supersonic inlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Supersonic_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                            CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose a supersonic outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Supersonic_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                             CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                  CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Set the initial conditions.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] ExtIter - External iteration.
    */
   void SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig *config, unsigned long ExtIter) override;


   /*!
    * \brief Update the solution static const unsigned a Runge-Kutta scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
    */
   void ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep) override;

   /*!
    * \brief Update the solution static const unsigned the explicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
   void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

   /*!
    * \brief Update the solution static const unsigned an implicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
   void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

 protected:
   /*!
  	* \brief Looking for non physical points in the initial solution.
  	* \param[in] config - Definition of the particular problem.
  	*/
   void Check_FreeStream_Solution(CConfig* config);

   /*!
    * \brief Looking for coherence in species order definition.
    * \param[in] config - Definition of the particular problem.
    */
   void Check_FreeStream_Species_Order(CConfig* config);

   /*!
 	  * \brief Reading files in case of restart
 	  * \param[in] geometry - Geometrical definition of the problem.
 	  * \param[in] config - Definition of the particular problem.
 	  */
   virtual void Load_Restart(CGeometry* geometry, CConfig* config);

   /*!
 	  * \brief Reading files in case of restart
 	  * \param[in] geometry - Geometrical definition of the problem.
 	  * \param[in] config - Definition of the particular problem.
 	  */
   void Read_SU2_Restart_Metadata(CGeometry* geometry, CConfig* config);

};

/*! \class CReactiveNSSolver
 *  \brief Main class for defining a solver for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSSolver:public CReactiveEulerSolver {
protected:
  su2double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */

  RealVec Xs_i, Xs_j;  /*!< \brief Auxiliary vectors to store mole fractions at node i and j. */
  RealVec Xs;         /*!< \brief Auxiliary vector to store mole fractions. */

protected:
  unsigned short RHOS_INDEX_GRAD; /*!< \brief Index for position of mole fractions in primitive gradient. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSSolver(): CReactiveEulerSolver(), Viscosity_Inf(), RHOS_INDEX_GRAD() {}

	/*!
	 * \overloaded Constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CReactiveNSSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSSolver() {}

  /*!
   * \brief Set gradient primitive variables static const unsigned Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned weighted least squares.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) override;

   /*!
    * \brief Set the fluid solver nondimensionalization.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] config - Definition of the particular problem.
    */
   void SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
   void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh, unsigned short iRKStep,
                      unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Compute the viscous residuals.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
   void Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Impose an isothermal wall boundary condition (no-slip).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                          CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the regression boundary condition for fuel inflow.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective part.
   * \param[in] visc_numerics - Description of the numerical method for visocus part.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Inflow(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                        CConfig* config, unsigned short val_marker) override;

protected:
  /*!
 	 * \brief Reading files in case of restart
 	 * \param[in] geometry - Geometrical definition of the problem.
 	 * \param[in] config - Definition of the particular problem.
 	 */
  void Load_Restart(CGeometry* geometry, CConfig* config) override;

};

#endif
