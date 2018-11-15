#ifndef SU2_SOLVER_REACTIVE
#define SU2_SOLVER_REACTIVE

#include "solver_structure.hpp"
#include "variable_reactive.hpp"

#include <memory>

/*! \class CReactiveEulerSolver
 *  \brief Main class for defining a solver for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerSolver:public CSolver {
public:

  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveEulerVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;
  using SmartArr = CReactiveEulerVariable::SmartArr;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  unsigned short nSpecies; /*!< \brief Total number of species. */

  bool  space_centered,  /*!< \brief True if space centered scheeme used. */
        implicit,      /*!< \brief True if euler implicit scheme used. */
        least_squares,  /*!< \brief True if computing gradients by least squares. */
        second_order, /*!< \brief True if second order recosntruction is applied. */
        limiter; /*!< \brief True if limiting strategy is applied. */

  RealVec   Lower_Limit,   /*!< \brief Lower limit conserved variables. */
            Upper_Limit;   /*!< \brief Upper limit conserved variables. */

  su2double Density_Inf,       /*!< \brief Free stream density. */
            Pressure_Inf,		  /*!< \brief Free stream pressure. */
	          Temperature_Inf;   /*!< \brief Trans.-rot. free stream temperature. */

  RealVec   Velocity_Inf,  /*!< \brief Free stream flow velocity. */
            MassFrac_Inf;  /*!< \brief Free stream species mass fraction. */

  RealVec   Sol_i,  /*!< \brief Auxiliary vector for storing the solution at point i. */
            Sol_j,      /*!< \brief Auxiliary vector for storing the solution at point j. */
            Primitive_i,        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point i. */
            Primitive_j;        /*!< \brief Auxiliary nPrimVar vector for storing the primitive at point j. */

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
	CReactiveEulerSolver(std::shared_ptr<CGeometry> geometry,std::shared_ptr<CConfig> config,unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerSolver() {}

   /*!
 	 * \brief Looking for non physical points in the initial solution
 	 * \param[in] config - Definition of the particular problem.
 	 */
 	void Check_FreeStream_Solution(std::shared_ptr<CConfig> config);

  /*!
	 * \brief Reading files in case of restart
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - Name of the file for the restart
	 */
  //void Read_Restart(std::shared_ptr<CGeometry> geometry,std::shared_ptr<CConfig>,std::string val_filename);

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
   * \brief A virtual member to compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  //virtual void SetPrimitive_Gradient(std::shared_ptr<CConfig> config);

  /*!
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the maximum value of the eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry* geometry, CConfig* config) override;

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
  void Set_MPI_Solution(CGeometry* geometry,CConfig* config) override;

  /*!
   * \brief Call MPI to set limiter of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Limiter(CGeometry* geometry,CConfig* config) override;

  /*!
   * \brief Call MPI to set gradient of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Gradient(CGeometry* geometry,CConfig* config) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned short iRKStep,
                     unsigned short RunTime_EqSystem, bool Output) override;

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
    * \brief Impose via the residual the Euler wall boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config,
                      unsigned short val_marker) override;

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
    * \brief Impose the inlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics, CConfig* config,
                 unsigned short val_marker) override;

   /*!
    * \brief Impose the outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics, CConfig* config,
                  unsigned short val_marker) override;

   /*!
    * \brief Impose the symmetry plane boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                     CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

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
   void ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                             unsigned short iRKStep) override;

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

};

/*! \class CReactiveNSSolver
 *  \brief Main class for defining a solver for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSSolver:public CReactiveEulerSolver {
protected:
  unsigned short nPrimVarAvgGrad; /*!< \brief Number of varaibles for average gradient */

  su2double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
   CReactiveNSSolver(): CReactiveEulerSolver(),nPrimVarAvgGrad(),Viscosity_Inf() {}

	/*!
	 * \overloaded constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CReactiveNSSolver(std::shared_ptr<CGeometry> geometry, std::shared_ptr<CConfig> config, unsigned short iMesh);

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
    * \brief Call MPI to set gradient of primitive variables in case of parallel simulation.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] config - Definition of the particular problem.
    */
   void Set_MPI_Primitive_Gradient(CGeometry* geometry,CConfig* config) override;

  /*!
   * \brief \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   //void SetPrimitive_Gradient(std::shared_ptr<CConfig> config) override;

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
   void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                      unsigned short iMesh, unsigned short iRKStep,
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
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective term.
   * \param[in] visc_numerics - Description of the numerical method for viscous term.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
  */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

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
                          CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) override;

};

#endif
