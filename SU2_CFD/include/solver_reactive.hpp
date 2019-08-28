#ifndef SU2_SOLVER_REACTIVE
#define SU2_SOLVER_REACTIVE

#include "solver_structure.hpp"
#include "variable_reactive.hpp"

/*! \class CReactiveEulerSolver
 *  \brief Main class for defining a solver for multispecies (chemically reacting or not) inviscid flows.
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

  su2double*** CharacPrimVar;  /*!< \brief Value of the characteristic variables at each boundary.
                                     NOTE: These will be use in case of multiphysics simulations. */

  RealVec   Lower_Limit,   /*!< \brief Lower limit values for conserved variables. */
            Upper_Limit;   /*!< \brief Upper limit values for conserved variables. */

  su2double Density_Inf,       /*!< \brief Free-stream density. */
            Pressure_Inf,		  /*!< \brief Free-stream pressure. */
	          Temperature_Inf;   /*!< \brief Translational free-stream temperature. */

  RealVec   Velocity_Inf,  /*!< \brief Free-stream flow velocity. */
            MassFrac_Inf;  /*!< \brief Free-stream species mass fraction. */

  RealVec   PrimVar_i,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point i. */
            PrimVar_j,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point j. */
            PrimVar_Vertex;  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at boundary node. */

  RealVec   Prim_i,       /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point i. */
            Prim_j,       /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point j. */
            Primitive;    /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at boundary node. */

  RealVec   Primitive_i,    /*!< \brief Auxiliary nPrimVar vector for storing primitive at point i in case of second order. */
            Primitive_j,    /*!< \brief Auxiliary nPrimVar vector for storing primitive at point j in case of second order. */
            Secondary_i,    /*!< \brief Auxiliary nVar vector for storing pressure derivatives at node i for 2nd order in implicit case. */
            Secondary_j;    /*!< \brief Auxiliary nVar vector for storing pressure derivatives at node j for 2nd order in implicit case. */

  RealVec   Buffer_Receive_U, /*!< \brief Auxiliary vector to receive information in case of parallel simulation. */
            Buffer_Send_U;    /*!< \brief Auxiliary vector to send information in case of parallel simulation. */

  RealVec   Ys_i,       /*!< \brief Auxiliary vector to store mass fractions at node i. */
            Ys_j,       /*!< \brief Auxiliary vector to store mass fractions at node j. */
            Ys;         /*!< \brief Auxiliary vector to store mass fractions whenever needed. */

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
	 * \overloaded Constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerSolver();

  /*!
   * \brief Get the pointer to the external library for physical-chemical properties
   */
  inline static LibraryPtr GetLibrary(void) {
    return library;
  }

  /*!
   * \brief Get the value of the characteristic variables at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the variables at the boundaries.
   */
  inline su2double* GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex) override {
    return CharacPrimVar[val_marker][val_vertex];
  }

  /*!
   * \brief Set primitive variables in each point reporting non physical data.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) override;

  /*!
   * \brief Set the gradient of primitive variables using Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the gradient primitive variables using weighted least squares.
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
   * \brief Call MPI to set solution in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Solution(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Call MPI to set old solution in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Solution_Old(CGeometry* geometry, CConfig* config) override;

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
   * \brief Solver preprocessing. This function in particular sets the primitive variables according to current state in all nodes.
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
   * \brief Set the residual in case of dual time simulations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   */
   void SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep,
                             unsigned short iMesh, unsigned short RunTime_EqSystem) override;

  /*!
   * \brief Compute the spatial integration using a centered scheme.
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
   * \brief Compute the spatial integration using an upwind scheme.
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
    * \brief Impose via the residual the Euler wall boundary condition (\bm{u}\cdot\bm{n} = 0).
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
                     CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the subsonic inlet boundary condition.
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
    * \brief Impose the outlet boundary condition (supersonic or subsonic according to current state).
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
    * \brief Set the initial condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] ExtIter - External iteration.
    */
    void SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig *config, unsigned long ExtIter) override;

   /*!
    * \brief Update the solution using a Runge-Kutta scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
    */
    void ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep) override;

   /*!
    * \brief Update the solution using the explicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
    void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

   /*!
    * \brief Update the solution using the implicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
    void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

 protected:
   /*!
  	* \brief Look for eventual non physical points in the initial solution.
  	* \param[in] config - Definition of the particular problem.
  	*/
    void Check_FreeStream_Solution(CConfig* config);

   /*!
    * \brief Check coherence in species order definition between the configuration file and the library.
    * \param[in] config - Definition of the particular problem.
    */
    void Check_FreeStream_Species_Order(CConfig* config);

   /*!
 	  * \brief Read files in case of restart.
 	  * \param[in] geometry - Geometrical definition of the problem.
 	  * \param[in] config - Definition of the particular problem.
 	  */
    virtual void Load_Restart(CGeometry* geometry, CConfig* config);

   /*!
 	  * \brief Read files in case of restart.
 	  * \param[in] geometry - Geometrical definition of the problem.
 	  * \param[in] config - Definition of the particular problem.
 	  */
    void Read_SU2_Restart_Metadata(CGeometry* geometry, CConfig* config);

};

/*! \class CReactiveNSSolver
 *  \brief Main class for defining a solver for multispecies (chemically reacting or not) viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSSolver:public CReactiveEulerSolver {
public:
  using Vec = Eigen::VectorXd;

protected:
  su2double Viscosity_Inf;	/*!< \brief Free-stream viscosity. */

  RealVec Xs_i,       /*!< \brief Auxiliary vectors to store mole fractions at node i. */
          Xs_j,       /*!< \brief Auxiliary vectors to store mole fractions at node j. */
          Xs;         /*!< \brief Auxiliary vector to store mole fractions. */

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
   * \brief Set the gradient of primitive variables using Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the gradient of primitive variables using weighted least squares.
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
   * \brief Solver preprocessing. This function in particular sets the primitive variables according to current state in all nodes.
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
   * \brief Impose a prescribed heat flux at wall (no-slip).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
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
 	 * \brief Read files in case of restart.
 	 * \param[in] geometry - Geometrical definition of the problem.
 	 * \param[in] config - Definition of the particular problem.
 	 */
  void Load_Restart(CGeometry* geometry, CConfig* config) override;

private:
  /*!
   * \brief Compute the diffusive flux along a certain direction solving Stefan-Maxwell equations.
   * \param[in] val_density - Density of the mixture.
   * \param[in] val_alpha - Parameter for artifical diffusion.
   * \param[in] val_Dij - Harmonic average of binary diffusion coefficients.
   * \param[in] val_xs - Molar fractions.
   * \param[in] val_grad_xs - Component along the desired dimension of gradient of molar fractions.
   * \param[in] val_ys - Mass fractions.
   * \return Multispecies diffusion flux along a certain direction
   */
  Vec Solve_SM(const su2double val_density, const su2double val_alpha, const RealMatrix& val_Dij,
               const RealVec& val_xs, const Vec& val_grad_xs, const RealVec& val_ys);
};

#endif
