

// Interfaces for solvers
#include "solverinterfaces/NodusSMOKE_ReactorNetwork_sparseNLSinterface.h"
#include "solverinterfaces/NodusSMOKE_ReactorNetwork_sparseODEinterface.h"
#include "solverinterfaces/NodusSMOKE_ReactorNetwork_denseNLSinterface.h"
#include "solverinterfaces/NodusSMOKE_ReactorNetwork_denseODEinterface.h"

// Math & Algebrae
#include "math/native-nls-solvers/NonLinearSystemSolver"
#include "math/native-ode-solvers/MultiValueSolver"
#include "Eigen/SparseLU"

namespace NodusSMOKE
{
	// CONSTRUCTOR
	ReactorNetwork_Gas::ReactorNetwork_Gas(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap, 
		std::vector<NodusSMOKE::Units*> *DeviceMap, std::vector<NodusSMOKE::StreamInfo>* StreamMap, 
		OpenSMOKE::OpenSMOKEVectorDouble FirstGuessY, std::vector<Eigen::Triplet<int>>* InputTripletsVector, 
		std::vector<Eigen::Triplet<double>>* InputTripletsVector_UnitMatrix, bool input_flag_sequential_solver)
		:
		ReactorNetwork(thermodynamicsMap, kineticsMap, DeviceMap, StreamMap, InputTripletsVector),
		TripletsVector_UnitMatrix(*InputTripletsVector_UnitMatrix)

		//transportMapXMLinternal(*transportMap)

	{

		calls = 0;
		
		// Case switch
		solid_case = false;
		sequential_solver = input_flag_sequential_solver;

		// Print initialize
		NLSiter = 0;
		n_step_video_NLS_ = 50;
		counter_file_video_NLS_ = 0;

		ODEiter = 0;
		n_step_video_ODE_ = 50;
		counter_file_video_ODE_ = 0;


		NS = thermodynamicsMapXMLinternal.NumberOfSpecies();
		NSgas = NS;
		NSsolid = 0;
		Nunits = Device.size();
		Nstreams = Stream.size();

		Nreactors = 0;
		Nreactors_HE = 0;
		for (int i = 0; i < Nunits; i++) {
			if (Device[i]->GetTag() == "Reactor") {
				++Nreactors;
				if (Device[i]->GetEnergy() == "NonIsothermal") {
					++Nreactors_HE;
				}
			}
		}

		neq = NS*Nreactors + Nreactors;

		// Memory allocation
		OpenSMOKE::ChangeDimensions(neq, &y0, true);
		OpenSMOKE::ChangeDimensions(neq, &yf, true);
		OpenSMOKE::ChangeDimensions(NSgas + 2, &residuals, true);



		// Build read-only Flow Matrix
		FlowMatrix.resize(Nunits, Nstreams);
		Eigen::SparseMatrix<int> FlowMatrix_temp(Nunits + 1, Nstreams);
		FlowMatrix_temp.setFromTriplets(TripletsVector.begin(), TripletsVector.end());
		FlowMatrix.operator=(FlowMatrix_temp.block(1, 0, Nunits, Nstreams));

		// Initialization
		y0 = FirstGuessY;

		if (sequential_solver == false) // If sequential solver is not wanted, I can compute with SparseLU all the mass flow and reduce the NLS
		{
		// Find best flowdevices order for no lag if possible
		ordered_devices = false;
			if (Nreactors < Nunits) {
				ordered_devices = true;
				FindOptimalDeviceOrder();
			}
			DistributeMassFlows();
		}

	};


	// DISTRIBUTE MASS FLOWS - Automatically distributes the mass flow within the network which should be a trivial task by solving a linear system

	int ReactorNetwork_Gas::DistributeMassFlows()
	{
		Eigen::VectorXd x(Nunits), b(Nunits);

		// Matrix A will need to be read: is the Unit/Unit matrix ==== we will call it UnitMatrix
		Eigen::SparseMatrix<double> UnitMatrix(Nunits, Nunits);

		if (TripletsVector_UnitMatrix.size() > 1) {

			Eigen::SparseMatrix<double> UnitMatrix_temp(Nunits + 1, Nunits + 1);
			UnitMatrix_temp.setFromTriplets(TripletsVector_UnitMatrix.begin(), TripletsVector_UnitMatrix.end());
			UnitMatrix.operator=(UnitMatrix_temp.block(1, 1, Nunits, Nunits));

		}
		else if (TripletsVector_UnitMatrix.size() <= 1) {

			UnitMatrix.setFromTriplets(TripletsVector_UnitMatrix.begin(), TripletsVector_UnitMatrix.end());

		}

		b.setZero();

		// Identify the inlets and assign them to the correct element of B - made in prediction of the possibility of using a mixer as first unit
		std::vector<std::vector<double>> b_partial(Nunits);
		for (unsigned int i = 0; i<Nstreams;++i){
			if (FlowMatrix.col(i).sum() == 1) {
				for (unsigned int j = 0; j<Nunits; ++j) {
					if (FlowMatrix.coeff(j, i) != 0) {
						b_partial[j].push_back(-Stream[i].mass_flow_rate_gas);
						std::for_each(b_partial[j].begin(), b_partial[j].end(), [&](double n) {
							b[j] += n;
						});
					}
				}
			}
		}

		// Solve the LS by LU factorization == Sparse-dedicated method
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(UnitMatrix);
		x = solver.solve(b);

		// Asssign mass flow rate to the units and from that to streams
		for (int j = 0; j < Nunits; ++j) {
			if (Device[j]->GetTag() == "Splitter") {
				// Find splits of this splitter
				std::vector<double> split_factor;
				for (int to=0; to< Nunits; ++to) { // Check all the units it goes to - hence "to" as iterator xD
					if (UnitMatrix.coeff(to,j) > 0.){
						split_factor.push_back(UnitMatrix.coeff(to,j));
					}
				}
				for (int i = 0; i < split_factor.size(); ++i) {
					Stream[Device[j]->GetOutlets()[i]].mass_flow_rate_gas = x(j)*split_factor[i];
				}
			}
			else {

				Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas = x(j);
			}
		}

		return 0;

	};

	// SET VALUES - Assign first guesses to each stream

	void ReactorNetwork_Gas::SetValues(OpenSMOKE::OpenSMOKEVectorDouble itervalues)
	{
		// Assign solution vector
		// Equation counter
		int k = 0;
		
		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			// if its a reactor
			if (Device[j]->GetTag() == "Reactor") {
				std::cout << Device[j]->GetName() << std::endl;
				for (int i = 1; i <= NSgas; ++i) {
					Stream[Device[j]->GetOutlets()[0]].omega_gas[i] = itervalues[++k]; // GAS SPECIES RESIDUAL
				}
				Stream[Device[j]->GetOutlets()[0]].temperature = itervalues[++k]; // T
				
				/*
				if (tracking == false) {
					for (int i = 1; i <= NSgas; ++i) {
						Stream[Device[j]->GetOutlets()[0]].omega_gas[i] = itervalues[++k]; // GAS SPECIES RESIDUAL
					}
					Stream[Device[j]->GetOutlets()[0]].temperature = itervalues[++k]; // T
				}

				else if (tracking == true) {
					if (Device[j]->GetType() == "PSR"){
						for (int i = 1; i <= NSgas; ++i) {
							Stream[Device[j]->GetOutlets()[0]].omega_gas[i] = itervalues[++k]; // GAS SPECIES RESIDUAL
						}
						Stream[Device[j]->GetOutlets()[0]].temperature = itervalues[++k]; // T
					}
					else if (Device[j]->GetType() == "PFR") {
						// Don't act on these, so skip these values
						for (int i = 1; i<= NSgas + 1; ++i)
							++k;
					}
				}
				*/
			}
		}
		
	};

	// NLS EQUATIONS
	int ReactorNetwork_Gas::NLSEquations(OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy) // COMMENTED DUE TO OBSOLETE - WORKING ON THE ODE SYSTEM
	{

		SetValues(y);

		// Solve auxiliary units to provide checking values
		if (ordered_devices == true )
			FlowDevicesOrderedProcessing();
		else
			FlowDevicesProcessing();

		// Run the Unit to obtain outlet of each
		ReactorProcessing();

		// Residual generator
		int k = 0;
		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			if (Device[j]->GetTag() == "Reactor") { // if its a reactor
				// Equations
				Device[j]->GetResiduals(residuals, Stream);
				for (int i = 1; i <= NSgas; ++i) { 
					dy[++k] = residuals[i]; // GAS SPECIES RESIDUAL
				}
				dy[++k]=residuals[NSgas+2]; // T							
			}
		}		

		return 0;
	};

	// ODE EQUATIONS
	int ReactorNetwork_Gas::ODEEquations(const double t, OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy)
	{

		SetValues(y);

		// Solve auxiliary units to provide checking values
		if (ordered_devices == true )
			FlowDevicesOrderedProcessing();
		else
			FlowDevicesProcessing();

		// Run the Unit to obtain outlet of each
		ReactorProcessing();

		// Residual generator
		int k=0;
		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			if (Device[j]->GetTag() == "Reactor"){ // if its a reactor
				
				// Equations
				if (tracking == false){
					Device[j]->GetResiduals(residuals, Stream);
				}
				else if (tracking == true){
					Device[j]->RTD(residuals, t, Stream);
				}
				for (int i = 1; i <= NSgas; ++i) { 
					dy[++k] = residuals[i]; // GAS SPECIES RESIDUAL
					std::cout << dy[k] << std::endl;
				}
				dy[++k]=residuals[NSgas+2]; // T
									
			}
		}


		return 0;

	};

	// NLS SOLVE
	void ReactorNetwork_Gas::SolveNLS()  // SERVE A RUNNARE IL SOLVER NLS
	{
		if (use_sparse_solvers == true) {
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Sparse NLS solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i)
				yMin(i) = 0.0;

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd nls_tolerances(neq);
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					nls_tolerances(k-1)=1.e-9;
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
			typedef NlsSMOKE::KernelSparse<NodusSMOKE_ReactorNetwork_sparseNLSinterface> kernel;
			NlsSMOKE::NonLinearSolver<kernel> nls_solver;
			nls_solver.SetReactor(this);

			// Set first guesses and solver options
			nls_solver.SetFirstGuessSolution(y0_eigen);
			nls_solver.SetMinimumValues(yMin);
			nls_solver.SetMaximumValues(yMax);
			nls_solver.SetPrint(true);

			// Set tolerances
			nls_solver.SetAbsoluteTolerances(nls_tolerances);
			
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
		else if (use_sparse_solvers == false) {
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Dense NLS solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i)
				yMin(i) = 0.0;

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd nls_tolerances(neq);
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					nls_tolerances(k-1)=1.e-9;
				}
				yMax(++k-1) = 100000.; // T
				nls_tolerances(k-1)=1.e-3;
			}

			// First guesses
			Eigen::VectorXd y0_eigen(y0_.Size());
			y0_.CopyTo(y0_eigen.data());

			// Final solution
			Eigen::VectorXd yf_eigen(y0_eigen.size());

			// Define the solver
			typedef NlsSMOKE::KernelDense<NodusSMOKE_ReactorNetwork_denseNLSinterface> kerneldense;
			NlsSMOKE::NonLinearSolver<kerneldense> nls_solver;
			nls_solver.SetReactor(this);

			// Set first guesses and solver options
			nls_solver.SetFirstGuessSolution(y0_eigen);
			nls_solver.SetMinimumValues(yMin);
			nls_solver.SetMaximumValues(yMax);
			nls_solver.SetPrint(true);

			// Set tolerances
			nls_solver.SetAbsoluteTolerances(nls_tolerances);
		
			
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
	void ReactorNetwork_Gas::SolveODE(const double tf)  // SERVE A RUNNARE IL SOLVER ODE
	{
		if (use_sparse_solvers == true){
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "####################################" << std::endl;
			std::cout << " Dense ODE system solution (OpenSMOKE++)..." << std::endl;
			std::cout << "####################################" << std::endl;

			OpenSMOKE::OpenSMOKEVectorDouble y0_ = y0;

			//// Min and max values
			Eigen::VectorXd yMin(neq);
			for (unsigned int i = 0; i < neq; ++i)
				yMin(i) = 0.0;

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd ode_tolerances(neq);
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					ode_tolerances(k-1)=1.e-9;
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
			typedef OdeSMOKE::KernelSparse<NodusSMOKE_ReactorNetwork_sparseODEinterface> sparseOde;
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

			// Set absolute tolerance
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
				// TO CONTINUE WITH NLS (A BETTER WAY SURELY EXISTS)
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
			for (unsigned int i = 0; i < neq; ++i)
				yMin(i) = 0.0;

			Eigen::VectorXd yMax(neq);
			Eigen::VectorXd ode_tolerances(neq);
			int k=0;
			for (unsigned int i=0; i<Nreactors; ++i){
				for (unsigned int j=0; j<NSgas; ++j){;
					yMax(++k-1)=1.;
					ode_tolerances(k-1)=1.e-9;
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
			typedef OdeSMOKE::KernelDense<NodusSMOKE_ReactorNetwork_denseODEinterface> denseOde;
			typedef OdeSMOKE::MethodGear<denseOde> methodGear;
			OdeSMOKE::MultiValueSolver<methodGear> ode_solver;
			ode_solver.SetReactor(this);
			ode_solver.SetPrint(true);

			// Set initial conditions
			ode_solver.SetInitialConditions(0., y0_eigen);

			// Set minimum and maximum values
			ode_solver.SetMinimumValues(yMin);
			ode_solver.SetMaximumValues(yMax);

			// Set absolute tolerance
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
				// TO CONTINUE WITH NLS (A BETTER WAY SURELY EXISTS)
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


	void ReactorNetwork_Gas::SetMWsAndMoleFractions(){
		for (unsigned int i = 0; i < Nstreams; ++i){
			OpenSMOKE::ChangeDimensions(NS, &Stream[i].x_gas, true);
			thermodynamicsMapXMLinternal.MoleFractions_From_MassFractions(Stream[i].x_gas.GetHandle(), Stream[i].MW_gas, Stream[i].omega_gas.GetHandle());								
		}
	};

	void ReactorNetwork_Gas::PrepareStreamSummary(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii) {
		
		// Prepare output file
		widths_of_output_species_.resize(NS);
		for (unsigned int i = 0; i<NS; i++)
			widths_of_output_species_[i] = OpenSMOKE::CalculateSpeciesFieldWidth(thermodynamicsMapXMLinternal.NamesOfSpecies()[i], NS);


		fOutput.open(output_file_ascii.c_str(), std::ios::out);

		unsigned int counter = 1;
		fOutput.setf(std::ios::scientific);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Stream # ", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "Phase ", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "T[K]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "P[Pa]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "MW[kg/kmol]", counter);
		OpenSMOKE::PrintTagOnASCIILabel(25, fOutput, "massflow[kg/s]", counter);

		for (unsigned int i = 0; i<NS; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMapXMLinternal.NamesOfSpecies()[i] + "_w", counter);
		for (unsigned int i = 0; i<NS; i++)
			OpenSMOKE::PrintTagOnASCIILabel(widths_of_output_species_[i], fOutput, thermodynamicsMapXMLinternal.NamesOfSpecies()[i] + "_x", counter);

		fOutput << std::endl;
	};
		
	void ReactorNetwork_Gas::PrintStreamSummary(std::ostream& fOutput)
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

					for (unsigned int j = 1; j <= NS; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].omega_gas[j];
					for (unsigned int j = 1; j <= NS; j++)
						fOutput << std::setw(widths_of_output_species_[j - 1]) << std::left << Stream[i].x_gas[j];

					fOutput << std::endl;
				}

			}

		}
	};

	void ReactorNetwork_Gas::SetTrack(std::string track_this){
		tracked_species = track_this;
		tracking = true;
		
	};

	void ReactorNetwork_Gas::PrepareTrackHistory(boost::filesystem::path output_folder_track) {
		
		
		OpenSMOKE::CreateDirectory(output_folder_track);
		Eigen::VectorXi FlowMatrixColumn_temporary;
		for (int i = 0; i<Nstreams; ++i) {
			FlowMatrixColumn_temporary.operator=(FlowMatrix.col(i));
			if (FlowMatrixColumn_temporary.maxCoeff() <= 0.) {
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
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutTrack, tracked_species +"_w_St"+ std::to_string(Stream[outlet_streams[i]].name), counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutTrack, tracked_species +"_x_St"+ std::to_string(Stream[outlet_streams[i]].name), counter);
			OpenSMOKE::PrintTagOnASCIILabel(25, fOutTrack, tracked_species +"_massflow_St"+ std::to_string(Stream[outlet_streams[i]].name), counter);
		}
		fOutTrack << std::endl;
		fOutTrack << std::endl;

	};

	void ReactorNetwork_Gas::PrintTrack(double t){

		fOutTrack <<  std::setw(25) << std::left << t;
		for (unsigned int i=0; i<outlet_streams.size(); i++ ){
				fOutTrack << std::setw(25) << std::left << Stream[outlet_streams[i]].omega_gas[thermodynamicsMapXMLinternal.IndexOfSpecies(tracked_species)];
				fOutTrack << std::setw(25) << std::left << Stream[outlet_streams[i]].x_gas[thermodynamicsMapXMLinternal.IndexOfSpecies(tracked_species)];
				fOutTrack << std::setw(25) << std::left << Stream[outlet_streams[i]].omega_gas[thermodynamicsMapXMLinternal.IndexOfSpecies(tracked_species)]*Stream[outlet_streams[i]].mass_flow_rate_gas;
		}
		fOutTrack << std::endl;

	};


} /* End of namespace OpenSMOKE */
