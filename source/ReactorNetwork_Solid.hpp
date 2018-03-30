

// Interfaces for solvers
#include "solverinterfaces/NetSMOKE_ReactorNetwork_sparseNLSinterface.h"
#include "solverinterfaces/NetSMOKE_ReactorNetwork_sparseODEinterface.h"
#include "solverinterfaces/NetSMOKE_ReactorNetwork_denseNLSinterface.h"
#include "solverinterfaces/NetSMOKE_ReactorNetwork_denseODEinterface.h"
// Math
#include "math/native-nls-solvers/NonLinearSystemSolver"
#include "math/native-ode-solvers/MultiValueSolver"

namespace NetSMOKE
{

	// CONSTRUCTOR
	ReactorNetwork_Solid::ReactorNetwork_Solid(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap, 
												OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMapSolid, OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMapSolid, 
												double rho_solid_input, std::vector<NetSMOKE::Units*> *DeviceMap, std::vector<NetSMOKE::StreamInfo>* StreamMap,
												OpenSMOKE::OpenSMOKEVectorDouble FirstGuessY, std::vector<Eigen::Triplet<int>>* InputTripletsVector, bool input_flag_sequential_solver)
		:
			ReactorNetwork(thermodynamicsMap, kineticsMap, DeviceMap, StreamMap, InputTripletsVector),
			thermodynamicsMapXMLinternalSOLID(*thermodynamicsMapSolid),
			kineticsMapXMLinternalSOLID(*kineticsMapSolid)
	{
		solid_case = true;
		sequential_solver = input_flag_sequential_solver;
		// Print initialize
		NLSiter = 0;
		n_step_video_NLS_ = 50;
		counter_file_video_NLS_ = 0;

		ODEiter = 0;
		n_step_video_ODE_ = 50;
		counter_file_video_ODE_ = 0;

		NS = thermodynamicsMapXMLinternalSOLID.NumberOfSpecies();
		NSgas = thermodynamicsMapXMLinternalSOLID.number_of_gas_species();
		NSsolid = thermodynamicsMapXMLinternalSOLID.number_of_solid_species();
		Nunits = Device.size();
		Nstreams = Stream.size();

		Nreactors = 0;
		Nreactors_HE = 0;
		Nreactors_Solid = 0;
		for (int i = 0; i < Nunits; i++){
			if (Device[i]->GetTag() == "Reactor"){
				++Nreactors;
				if (Device[i]->GetEnergy() == "NonIsothermal") {
					++Nreactors_HE;
				}
				if (Device[i]->GetPhase() == "Solid" || Device[i]->GetPhase() == "Mix") {
					++Nreactors_Solid;
				}
			}
		}

		neq = Nreactors*(NSgas + 2) + Nreactors_Solid*(NSsolid + 1);;

		// Memory allocation
		OpenSMOKE::ChangeDimensions(neq, &y0, true);
		OpenSMOKE::ChangeDimensions(neq, &yf, true);
		OpenSMOKE::ChangeDimensions(NS+3, &residuals, true); // Since it is only one constantly reused, make sure to have enough space for any case

		// Set solid density
		rho_solid = rho_solid_input;

		// Initialization
		y0 = FirstGuessY;

		// Build read-only Flow Matrix
		FlowMatrix.resize(Nunits,Nstreams);
		Eigen::SparseMatrix<int> FlowMatrix_temp(Nunits+1,Nstreams);
		FlowMatrix_temp.setFromTriplets(TripletsVector.begin(),TripletsVector.end());
		FlowMatrix.operator=(FlowMatrix_temp.block(1, 0, Nunits, Nstreams));

		// Find optimal auxiliary devices solving order

		if (sequential_solver == false) {
			ordered_devices = false;
			if (Nreactors < Nunits) {
				ordered_devices = true;
				FindOptimalDeviceOrder();
			}	
		}
		


	};


	// SET VALUES - Assign first guesses

	void ReactorNetwork_Solid::SetValues(OpenSMOKE::OpenSMOKEVectorDouble itervalues)
	{
		// Assign solution vector
		// Equation counter
		int k = 0;
		
		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			// if its a reactor
			if (Device[j]->GetTag() == "Reactor") {
				for (int i = 1; i <= NSgas; ++i) {
					Stream[Device[j]->GetOutlets()[0]].omega_gas[i] = itervalues[++k]; // GAS SPECIES RESIDUAL
				}
				if (Device[j]->GetPhase() != "Gas") {
					for (int i = 1; i <= NSsolid; ++i) {
						Stream[Device[j]->GetOutlets()[0]].omega_solid[i] = itervalues[++k]; // SOLID SPECIES RESIDUAL
					}
				}
				if (Device[j]->GetPhase() != "Gas") {
					Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas = itervalues[++k]; // Gas
					Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_solid = itervalues[++k]; // Solid
					Stream[Device[j]->GetOutlets()[0]].temperature = itervalues[++k]; // T
				}
				else {
					Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas = itervalues[++k]; // Gas
					Stream[Device[j]->GetOutlets()[0]].temperature = itervalues[++k]; // T
				}
				/*
				//TOREMOVE
				// Random debug checks
				std::cout << Unit[j].name + " " + Unit[j].phase << std::endl;
				std::cout << Unit[j].omega_in_gas.SumElements() << std::endl;
				std::cout << Unit[j].omega_in_solid.SumElements() << std::endl;
				*/
			}
		}
		

	};

	

	// NLS EQUATIONS
	int ReactorNetwork_Solid::NLSEquations(OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy)
	{

		// Assign current set of values to the System
		SetValues(y);

		// Solve auxiliary units to provide checking values
		if (ordered_devices == true )
			FlowDevicesOrderedProcessing();
		else
			FlowDevicesProcessing();

			
		// Run the Unit to obtain outlet of each
		ReactorProcessing();

		/// AT THIS POINT, A COMPLETE (BUT NOT CONVERGED) PICTURE OF THE SYSTEM IS ESTABLISHED ///

		// RESIDUAL GENERATION SECTION
		// Equation counter
		int k = 0;

		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			// if its a reactor
			if (Device[j]->GetTag() == "Reactor") {
				Device[j]->GetResiduals(residuals, Stream);

				for (int i = 1; i <= NSgas; ++i) { 
					dy[++k] = residuals[i]; // GAS SPECIES RESIDUAL
				}
				if (Device[j]->GetPhase() != "Gas"){
					for (int i = 1; i <= NSsolid; ++i) { 
						dy[++k] = residuals[i+NSgas]; // SOLID SPECIES RESIDUAL
					}
				}

				// Mass balances and T
				if (Device[j]->GetPhase() != "Gas"){
					dy[++k] = residuals [NS+1]; // Gas
					dy[++k] = residuals[NS+2]; // Solid
					dy[++k]=residuals[NS+3]; // T
				}
				else {
					dy[++k] = residuals [NSgas+1]; // Gas
					dy[++k]=residuals[NSgas+2]; // T
				}								
			}
		}

		// Double check
		if (k != neq)
			OpenSMOKE::FatalErrorMessage("Discrepancy between found between expected n° of equations and equations constructed!");

		// ITERATION COMPLETED

		// TOREMOVE Debug prints
		/*
		NLSiter++;
		std::cout << NLSiter << std::endl;
		*/

		return 0;
	};

	// ODE EQUATIONS
	int ReactorNetwork_Solid::ODEEquations(const double t, OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy)
	{
		// Assign current set of values to the System
		SetValues(y);

		// Solve auxiliary units to provide checking values
		if (ordered_devices == true )
			FlowDevicesOrderedProcessing();
		else
			FlowDevicesProcessing();

		// Run the Unit to obtain outlet of each
		ReactorProcessing();

		/// AT THIS POINT, A COMPLETE (BUT NOT CONVERGED) PICTURE OF THE SYSTEM IS ESTABLISHED ///

		// RESIDUAL GENERATION SECTION
		// Equation counter
		int k = 0;

		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			// if its a reactor
			if (Device[j]->GetTag() == "Reactor") {

				// Equations
				if (tracking == false){
					Device[j]->GetResiduals(residuals, Stream);
				}
				else if (tracking == true){
					Device[j]->RTD(residuals, t, Stream);
				}

				for (int i = 1; i <= NSgas; ++i) { 
					dy[++k] = residuals[i]; // GAS SPECIES RESIDUAL
				}

				if (Device[j]->GetPhase() != "Gas"){
					for (int i = 1; i <= NSsolid; ++i) { 
						dy[++k] = residuals[i+NSgas]; // SOLID SPECIES RESIDUAL
					}
				}

				// Mass balances and T
				if (Device[j]->GetPhase() != "Gas"){
					dy[++k] = residuals [NS+1]; // Gas
					dy[++k] = residuals[NS+2]; // Solid
					dy[++k] = residuals[NS+3]; // T
				}
				else {
					dy[++k] = residuals [NSgas+1]; // Gas
					dy[++k] = residuals[NSgas+2]; // T
				}						
			}
		}

		// Double check
		if (k != neq)
			OpenSMOKE::FatalErrorMessage("Discrepancy between found between expected n° of equations and equations constructed!");

		// ITERATION COMPLETED

		//TOREMOVE Debug prints
		/*
		ODEiter++;
		std::cout << ODEiter << std::endl;
		*/
		
		return 0;
	};

	

	// NLS SOLVE
	void ReactorNetwork_Solid::SolveNLS()  // SERVE A RUNNARE IL SOLVER NLS
	{
		if (use_sparse_solvers == true){
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Sparse NLS solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i) {
				yMin(i) = 0.0;
			}

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd nls_tolerances(neq);

			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					nls_tolerances(k-1)=1.e-9;
				}
				if (Device[i]->GetPhase() != "Gas"){
					for (unsigned int j=0; j<NSsolid; ++j){;
						yMax(++k-1)=1.;
						nls_tolerances(k-1)=1.e-9;
					}
				}
				yMax(++k-1)=100000.; 	// G mass
				nls_tolerances(k-1)=1.e-3;
				if (Device[i]->GetPhase() != "Gas"){
					yMax(++k-1)=100000.;	// S mass
					nls_tolerances(k-1)=1.e-3;
				}
				yMax(++k-1) = 100000.; // T
				nls_tolerances(k-1)=1.;
			}

			// First guesses
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Define the solver
			typedef NlsSMOKE::KernelSparse<NetSMOKE_ReactorNetwork_sparseNLSinterface> kernel;
			NlsSMOKE::NonLinearSolver<kernel> nls_solver;
			nls_solver.SetReactor(this);
			
			// Set first guesses and solver options
			nls_solver.SetFirstGuessSolution(y0_eigen);
			nls_solver.SetMinimumValues(yMin);
			nls_solver.SetMaximumValues(yMax);
			
			// Set tolerances
			nls_solver.SetAbsoluteTolerances(nls_tolerances);

			// Set print
			nls_solver.SetPrint(true);

			// Set sparsity pattern
			nls_solver.SetSparsityPattern(rows_sparsity_pattern, cols_sparsity_pattern, true);

			// Solve the non linear system
			double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			NlsSMOKE::NlsStatus status = nls_solver();
			double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Analyze the solution
			if (status >= 0)
			{
				std::string message("Nls System successfully solved: ");
				if (status == 0)		message += "Start conditions";
				else if (status == 1)	message += "The maximum number of functions calls has been performed";
				else if (status == 2)	message += "The Newton correction has reached the required precision";
				else if (status == 3)	message += "The Quasi Newton correction has reached the required precision";
				else if (status == 4)	message += "The Gradient has reached the required precision";
				else if (status == 5)	message += "The Objective function phiNew has reached the required precision";
				else if (status == 6)	message += "The Objective function phiW has reached the required precision";
				else if (status == 7)	message += "The Objective function phiNew has reached the required precision but solution is dubious";
				else if (status == 8)	message += "Reached the assigned max value for Newton calls";

				std::cout << message << std::endl;

				Eigen::VectorXd f(neq);
				nls_solver.Solution(yf_eigen, f);
				yf.CopyFrom(yf_eigen.data());

				//if (parameters.verbosity_level() > 0)
				nls_solver.NlsSummary(std::cout);


			}
			else
			{
				std::string message("Nls Solver Error: ");
				if (status == -1)		message += "It has been impossible to reach the solution";
				else if (status == -2)	message += "The search has been stopped";
				else if (status == -3)	message += "The object is not initialized";
				else if (status == -4)	message += "It has been impossible to reach the solution in Restart";

				std::cout << message << std::endl;
			}
		}
		if (use_sparse_solvers == false){
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Dense NLS solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i) {
				yMin(i) = 0.0;
			}

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd nls_tolerances(neq);

			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					nls_tolerances(k-1)=1.e-9;
				}
				if (Device[i]->GetPhase() != "Gas"){
					for (unsigned int j=0; j<NSsolid; ++j){;
						yMax(++k-1)=1.;
						nls_tolerances(k-1)=1.e-9;
					}
				}
				yMax(++k-1)=100000.; 	// G mass
				nls_tolerances(k-1)=1.e-3;
				if (Device[i]->GetPhase() != "Gas"){
					yMax(++k-1)=100000.;	// S mass
					nls_tolerances(k-1)=1.e-3;
				}
				yMax(++k-1) = 100000.; // T
				nls_tolerances(k-1)=1.;
			}

			// First guesses
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Define the solver
			typedef NlsSMOKE::KernelDense<NetSMOKE_ReactorNetwork_denseNLSinterface> kernel;
			NlsSMOKE::NonLinearSolver<kernel> nls_solver;
			nls_solver.SetReactor(this);

			
			
			// Set first guesses and solver options
			nls_solver.SetFirstGuessSolution(y0_eigen);
			nls_solver.SetMinimumValues(yMin);
			nls_solver.SetMaximumValues(yMax);
			
			// Set tolerances
			nls_solver.SetAbsoluteTolerances(nls_tolerances);

			// Set print
			nls_solver.SetPrint(true);

			// Solve the non linear system
			double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			NlsSMOKE::NlsStatus status = nls_solver();
			double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Analyze the solution
			if (status >= 0)
			{
				std::string message("Nls System successfully solved: ");
				if (status == 0)		message += "Start conditions";
				else if (status == 1)	message += "The maximum number of functions calls has been performed";
				else if (status == 2)	message += "The Newton correction has reached the required precision";
				else if (status == 3)	message += "The Quasi Newton correction has reached the required precision";
				else if (status == 4)	message += "The Gradient has reached the required precision";
				else if (status == 5)	message += "The Objective function phiNew has reached the required precision";
				else if (status == 6)	message += "The Objective function phiW has reached the required precision";
				else if (status == 7)	message += "The Objective function phiNew has reached the required precision but solution is dubious";
				else if (status == 8)	message += "Reached the assigned max value for Newton calls";

				std::cout << message << std::endl;

				Eigen::VectorXd f(neq);
				nls_solver.Solution(yf_eigen, f);
				yf.CopyFrom(yf_eigen.data());

				//if (parameters.verbosity_level() > 0)
				nls_solver.NlsSummary(std::cout);


			}
			else
			{
				std::string message("Nls Solver Error: ");
				if (status == -1)		message += "It has been impossible to reach the solution";
				else if (status == -2)	message += "The search has been stopped";
				else if (status == -3)	message += "The object is not initialized";
				else if (status == -4)	message += "It has been impossible to reach the solution in Restart";

				std::cout << message << std::endl;
			}
		}


	};



	// ODE SOLVE
	void ReactorNetwork_Solid::SolveODE(const double tf)  // SERVE A RUNNARE IL SOLVER ODE
	{
		if (use_sparse_solvers == true){
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Sparse ODE system solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i) {
				yMin(i) = 0.0;
			}

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd ode_tolerances(neq);
			
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					ode_tolerances(k-1)=1.e-9;
				}
				if (Device[i]->GetPhase() != "Gas"){
					for (unsigned int j=0; j<NSsolid; ++j){;
						yMax(++k-1)=1.;
						ode_tolerances(k-1)=1.e-9;
					}
				}
				yMax(++k-1)=100000.; 	// G mass
				ode_tolerances(k)=1.e-3;
				if (Device[i]->GetPhase() != "Gas"){
					yMax(++k-1)=100000.;	// S mass
					ode_tolerances(k-1)=1.e-3;
				}
				yMax(++k-1) = 100000.; // T
				ode_tolerances(k-1)=1.;
			}

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelSparse<NetSMOKE_ReactorNetwork_sparseODEinterface> sparseOde;
			typedef OdeSMOKE::MethodGear<sparseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);
			ode_solver.SetPrint(true);
			

			// Set sparsity pattern
			ode_solver.SetSparsityPattern(rows_sparsity_pattern, cols_sparsity_pattern);
			
			// Set initial conditions
			ode_solver.SetInitialConditions(0., y0_eigen);

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set tolerances
			ode_solver.SetAbsoluteTolerances(ode_tolerances);

			// Set norm of residuals stop criteria
			
			double dyNorm = 0.00001*(NS*Nreactors);
			ode_solver.SetStopConditionMaximumYPrimeNorm1(dyNorm);
			if (tracking == true) {
				std::cout << "TRACKING ON" << std::endl;
				ode_solver.UnsetStopConditions();
			}

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				std::string message("OpenSMOKE++ Ode System successfully solved: ");

				if (status == 1)		message += "INITIALIZATION_STATE";
				else if (status == 2)	message += "CONTINUATION_STATE";
				else if (status == 3)	message += "INTEGRATION_STOPPED_WHEN_SUM_ABS_Y1_IS_LESS_THAN";

				std::cout << message << std::endl;
				ode_solver.Solution(yf_eigen);
				yf.CopyFrom(yf_eigen.data());
				ode_solver.OdeMethodSummary(std::cout);
				// TO CONTINUE WITH NLS (MIGHT NEED A MORE ELEGANT SOLUTION I GUESS)
				y0.CopyFrom(yf.GetHandle());
			}
			else
			{
				std::string message("OpenSMOKE++ Ode Solver Error: ");
				if (status == -1)	message += "MAX_NUMBER_OF_STEPS_REACHED";
				else if (status == -2)	message += "TOO_STRICT_TOLERANCES";
				else if (status == -3)	message += "ILLEGAL_CONTINUATION_REQUEST";
				else if (status == -4)	message += "MAX_NUMBER_ERRORTEST_FAILURES";
				else if (status == -5)	message += "MAX_NUMBER_CONVERGENCETEST_FAILURES";
				else if (status == -6)	message += "TOO_SMALL_STEP_SIZE";
				else if (status == -7)	message += "ILLEGAL_MAX_INDEPENDENT_VARIABLE";
				else if (status == -9)	message += "ILLEGAL_CONSTRAINTS";
				else if (status == -10)	message += "EXCEPTION_HANDLING_STOP";
				else if (status == -12)	message += "YOU_CANNOT_OVERSHOOT_TCRITIC";

				std::cout << message << std::endl;
			}
		}
		else if (use_sparse_solvers == false){
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Dense ODE system solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i) {
				yMin(i) = 0.0;
			}

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd ode_tolerances(neq);
			
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					ode_tolerances(k-1)=1.e-9;
				}
				if (Device[i]->GetPhase() != "Gas"){
					for (unsigned int j=0; j<NSsolid; ++j){;
						yMax(++k-1)=1.;
						ode_tolerances(k-1)=1.e-9;
					}
				}
				yMax(++k-1)=100000.; 	// G mass
				ode_tolerances(k)=1.e-3;
				if (Device[i]->GetPhase() != "Gas"){
					yMax(++k-1)=100000.;	// S mass
					ode_tolerances(k-1)=1.e-3;
				}
				yMax(++k-1) = 100000.; // T
				ode_tolerances(k-1)=1.;
			}

			// Initial conditions
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Create the solver
			typedef OdeSMOKE::KernelDense<NetSMOKE_ReactorNetwork_denseODEinterface> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);
			ode_solver.SetPrint(true);
			
			// Set initial conditions
			ode_solver.SetInitialConditions(0., y0_eigen);

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set tolerances
			ode_solver.SetAbsoluteTolerances(ode_tolerances);

			// Set norm of residuals stop criteria		
			double dyNorm = 0.00001*(NS*Nreactors);
			ode_solver.SetStopConditionMaximumYPrimeNorm1(dyNorm);
			if (tracking == true) {
				std::cout << "TRACKING ON" << std::endl;
				ode_solver.UnsetStopConditions();
			}

			// Solve the system
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			OdeSMOKE::OdeStatus status = ode_solver.Solve(tf);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			// Check the solution
			if (status > 0)
			{
				std::string message("OpenSMOKE++ Ode System successfully solved: ");

				if (status == 1)		message += "INITIALIZATION_STATE";
				else if (status == 2)	message += "CONTINUATION_STATE";
				else if (status == 3)	message += "INTEGRATION_STOPPED_WHEN_SUM_ABS_Y1_IS_LESS_THAN";

				std::cout << message << std::endl;
				ode_solver.Solution(yf_eigen);
				yf.CopyFrom(yf_eigen.data());
				ode_solver.OdeMethodSummary(std::cout);
				// TO CONTINUE WITH NLS (MIGHT NEED A MORE ELEGANT SOLUTION I GUESS)
				y0.CopyFrom(yf.GetHandle());
			}
			else
			{
				std::string message("OpenSMOKE++ Ode Solver Error: ");
				if (status == -1)	message += "MAX_NUMBER_OF_STEPS_REACHED";
				else if (status == -2)	message += "TOO_STRICT_TOLERANCES";
				else if (status == -3)	message += "ILLEGAL_CONTINUATION_REQUEST";
				else if (status == -4)	message += "MAX_NUMBER_ERRORTEST_FAILURES";
				else if (status == -5)	message += "MAX_NUMBER_CONVERGENCETEST_FAILURES";
				else if (status == -6)	message += "TOO_SMALL_STEP_SIZE";
				else if (status == -7)	message += "ILLEGAL_MAX_INDEPENDENT_VARIABLE";
				else if (status == -9)	message += "ILLEGAL_CONSTRAINTS";
				else if (status == -10)	message += "EXCEPTION_HANDLING_STOP";
				else if (status == -12)	message += "YOU_CANNOT_OVERSHOOT_TCRITIC";

				std::cout << message << std::endl;
			}
		}
	};

	void ReactorNetwork_Solid::SetMWsAndMoleFractions(){

		for (unsigned int i = 0; i < Nstreams; ++i){
			OpenSMOKE::ChangeDimensions(NSgas, &Stream[i].x_gas, true);
			OpenSMOKE::ChangeDimensions(NSsolid, &Stream[i].x_solid, true);	

			if (Stream[i].mass_flow_rate_gas > 1.e-16){		
				double MWmix = 0;
				for (unsigned int j = 1; j <= NSgas; ++j) {
					MWmix += Stream[i].omega_gas[j] / thermodynamicsMapXMLinternalSOLID.MW(j - 1);
				}
				MWmix = 1./MWmix;
				Stream[i].MW_gas = MWmix;
				for (unsigned int j = 1; j <= NSgas; ++j) {
					Stream[i].x_gas[j] = Stream[i].omega_gas[j] * MWmix / thermodynamicsMapXMLinternalSOLID.MW(j - 1) ;
				}	
			}

			if (Stream[i].mass_flow_rate_solid > 1.e-16){
				double MWmix = 0;
				for (unsigned int j = 1; j <= NSsolid; ++j) {
					MWmix += Stream[i].omega_solid[j] / thermodynamicsMapXMLinternalSOLID.MW(NSgas + j - 1);
				}
				MWmix = 1./MWmix;
				Stream[i].MW_solid = MWmix;
				for (unsigned int j = 1; j <= NSsolid; ++j) {
					Stream[i].x_solid[j] = Stream[i].omega_solid[j] * MWmix / thermodynamicsMapXMLinternalSOLID.MW(NSgas + j - 1) ;
				}	
			}
		}

	};

	void ReactorNetwork_Solid::SetStreamPhase(int i)
	{	
		bool contains_gas = false;
		bool contains_solid = false;

		if (Stream[i].mass_flow_rate_gas > 1.e-16)
			contains_gas = true;
		if (Stream[i].mass_flow_rate_solid > 1.e-16)
			contains_solid = true;

		if (contains_gas == true && contains_solid == true){
			Stream[i].phase = "Mix";
		}
		else if (contains_gas == true && contains_solid == false){
			Stream[i].phase = "Gas";
		}
		else if (contains_gas == false && contains_solid == true){
			Stream[i].phase = "Solid";
		}
	};

	void ReactorNetwork_Solid::PrepareStreamSummary(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii) {
		
		// Prepare output file
		widths_of_output_species_.resize(NS);
		for (unsigned int i = 0; i<NS; i++)
			widths_of_output_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXMLinternalSOLID.NamesOfSpecies()[i], NS);


		fOutput.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		fOutput.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Stream # ", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Phase ", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Gas MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Gas massflow[kg/s]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Solid MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Solid massflow[kg/s]", counter);

		for (unsigned int i = 0; i<NS; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMapXMLinternalSOLID.NamesOfSpecies()[i] + "_w", counter);
		for (unsigned int i = 0; i<NS; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMapXMLinternalSOLID.NamesOfSpecies()[i] + "_x", counter);

		fOutput << std::endl;
	};
		
	void ReactorNetwork_Solid::PrintStreamSummary(std::ostream& fOutput)
	{
		
		// Prepare a sorting reference vector
		std::vector<int> StreamNames(Stream.size());
		for (unsigned int i = 0; i < Stream.size(); ++i) {
			StreamNames[i] = Stream[i].name;
		}
		std::sort(StreamNames.begin(), StreamNames.end());
	
		for (unsigned int k = 0; k < Stream.size(); ++k) {
			for (unsigned int i = 0; i < Stream.size(); ++i) {
				if (StreamNames[k] == Stream[i].name) {
					fOutput << std::setw(25) << std::left << Stream[i].name;
					fOutput << std::setw(25) << std::left << Stream[i].phase;
					fOutput << std::setw(25) << std::left << Stream[i].temperature;
					fOutput << std::setw(25) << std::left << Stream[i].pressure;
					fOutput << std::setw(25) << std::left << Stream[i].MW_gas;
					fOutput << std::setw(25) << std::left << Stream[i].mass_flow_rate_gas;
					fOutput << std::setw(25) << std::left << Stream[i].MW_solid;
					fOutput << std::setw(25) << std::left << Stream[i].mass_flow_rate_solid;

					for (unsigned int j = 1; j <= NSgas; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].omega_gas[j];
					for (unsigned int j = 1; j <= NSsolid; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].omega_solid[j];
					for (unsigned int j = 1; j <= NSgas; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].x_gas[j];
					for (unsigned int j = 1; j <= NSsolid; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].x_solid[j];

					fOutput << std::endl;
				}

			}

		}
	};

	void ReactorNetwork_Solid::SetTrack(std::string track_this){
		tracked_species = track_this;
		tracking = true;
		
	};

	void ReactorNetwork_Solid::PrepareTrackHistory(boost::filesystem::path output_folder_track) {
		
		OpenSMOKE::CreateDirectory(output_folder_track);
		Eigen::VectorXi FlowMatrixColumn_temporary;
		for (int i = 0; i<Nstreams; ++i){
			FlowMatrixColumn_temporary.operator=(FlowMatrix.col(i));
			if (FlowMatrixColumn_temporary.maxCoeff()<=0.){
				outlet_streams.push_back(i);
			}
		}

		std::string filename = tracked_species + "_track.out";
		boost::filesystem::path output_track_file = output_folder_track / filename;

		fOutTrack.open(output_track_file.c_str(), std::ios::out);

		fOutTrack.setf(std::ios::scientific);

		unsigned int counter = 1;
		fOutTrack.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutTrack, "t [s] ", counter);
		for (unsigned int i=0; i<outlet_streams.size(); i++ ){
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutTrack, tracked_species +"_w_St" + std::to_string(Stream[outlet_streams[i]].name), counter);
		}
		fOutTrack << std::endl;
		fOutTrack << std::endl;

	};

	void ReactorNetwork_Solid::PrintTrack(double t){

		fOutTrack <<  std::setw(25) << std::left << t;
		for (unsigned int i=0; i<outlet_streams.size(); i++ ){
			if (thermodynamicsMapXMLinternalSOLID.IndexOfSpecies(tracked_species) <= NSgas) {
				fOutTrack << std::setw(25) << std::left << Stream[outlet_streams[i]].omega_gas[thermodynamicsMapXMLinternalSOLID.IndexOfSpecies(tracked_species)];
			}
			else if (thermodynamicsMapXMLinternalSOLID.IndexOfSpecies(tracked_species) > NSgas) {
				fOutTrack << std::setw(25) << std::left << Stream[outlet_streams[i]].omega_solid[thermodynamicsMapXMLinternalSOLID.IndexOfSpecies(tracked_species)];
			}
		}
		fOutTrack << std::endl;

	};
	
} /* END NAMESPACE OPENSMOKE*/
