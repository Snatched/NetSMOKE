

// Include OpenMP Header file
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace NetSMOKE
{

	// CONSTRUCTOR
	ReactorNetwork::ReactorNetwork(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap,
		OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap,
		std::vector<NetSMOKE::Units*>* DeviceMap,
		std::vector<NetSMOKE::StreamInfo>* Stream,
		std::vector<Eigen::Triplet<int>>* InputTripletsVector)
		:
		thermodynamicsMapXMLinternal(*thermodynamicsMap),
		kineticsMapXMLinternal(*kineticsMap),
		Device(* DeviceMap),
		Stream(* Stream),
		TripletsVector(*InputTripletsVector)
	{
		tracking = false;

		// Create a backup copy of the associated input
		{
			// Get backup destination
			working_folder = boost::filesystem::current_path();
			backup_folder = working_folder / "backup";
			// Create backup folder
			OpenSMOKE::CreateDirectory(backup_folder);

			// Delete older input backup filename
			boost::filesystem::directory_iterator it(backup_folder);
			boost::filesystem::directory_iterator endit;
			while (it != endit) {
				if (it->path().extension() == ".dic") {
					boost::filesystem::remove(it->path());
				}
				++it;
			}
			// Copy input
			boost::filesystem::path input_file = working_folder / "input.dic";
			boost::filesystem::path backup_input_file = backup_folder / "input_backup.dic";
			boost::filesystem::copy_file(input_file, backup_input_file);
		}

	};

	
	ReactorNetwork::~ReactorNetwork() 
	{

		// CONTEXT: Device vector is a vector of heap-allocated objects. Not clearing this leads to memory corruption.
		for (int i = 0; i < Device.size(); ++i) {
			delete Device[i];
		}
		Device.clear();

	};	

	// REACTOR PROCESSING - Automatically constructs and solve various units within the network through the "units" classes family

	int ReactorNetwork::ReactorProcessing()
	{

		// Process reactors
		for (unsigned int j = 0; j < Nunits; ++j){
			if (Device[j]->GetTag() == "Reactor"){
					Device[j]->Solve(Stream);
			}
		}

		return 0;

	};

	int ReactorNetwork::ReactorSequentialProcessing()
	{

		// Process reactors
		for (unsigned int j = 0; j < Nunits; ++j){
			if (Device[j]->GetTag() == "Reactor"){
				Device[j]->NonIterativeSolve(Stream);
			}
		}

		return 0;

	};


	int ReactorNetwork::FlowDevicesProcessing()
	{

		// Once reactors are set and the outputs have been assigned, solve the trivial units ---> This is done due to the mixers and splitters ideally not requiring initial conditions, they only perform networking operation
		for (int j = 0; j < Nunits; ++j){
			if (Device[j]->GetTag() != "Reactor"){
				Device[j]->Solve(Stream);
			}
		}

		return 0;

	};

	int ReactorNetwork::FlowDevicesOrderedProcessing()
	{

		
		for (int j = 0; j < FlowDevicesOptimalOrder.size(); ++j) {
				Device[FlowDevicesOptimalOrder[j]]->Solve(Stream);	
		}

		return 0;

	};

	int ReactorNetwork::SolveAsSeries()
	{

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "##############################################" << std::endl;
		std::cout << " 			Series of reactor solution ..." << std::endl;
		std::cout << "##############################################" << std::endl;
		
		while ( SolvingOrder.size() < Device.size() ) {
			
			for (int j = 0; j < Device.size(); ++j) {
				
				// Check if the device has already been studied
				int current = j;
		
			
				bool is_it_already_considered = std::any_of(SolvingOrder.begin(), SolvingOrder.end(), [&](int m) { return m == current; });
			
				// Do the stuff only if flow device is not already considered into the optimal order
				if (is_it_already_considered == false) {
					
					std::vector<double> InletValues(Device[j]->GetInlets().size());
					for (int i = 0; i < Device[j]->GetInlets().size(); ++i) {
						InletValues[i]=Stream[Device[j]->GetInlets()[i]].mass_flow_rate_gas+Stream[Device[j]->GetInlets()[i]].mass_flow_rate_solid;
					}

					// Check if one of the inlets is zero - returns true if one of the inlets has mass flow = 0
					bool is_there_a_zero = std::any_of(InletValues.begin(), InletValues.end(), [](double l) { return l == 0.; });
	
					if (is_there_a_zero == false) {
						SolvingOrder.push_back(j);
						if (Device[j]->GetTag() != "Reactor"){
							std::cout << "Solving  " << Device[j]->GetTag() << " " << Device[j]->GetName() << "..." << std::endl;
							Device[j]->Solve(Stream);
						}
						else if (Device[j]->GetTag() == "Reactor"){
							std::cout << "Solving reactor " << Device[j]->GetName() << "..." << std::endl;
							Device[j]->NonIterativeSolve(Stream);
						}
					}
				}
			}
	
		}
		

		return 0;
	};

	int ReactorNetwork::SequentialSubstitution ( int iter )
	{

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "######################################" << std::endl;
		std::cout << " Subsequential substitution method..." << std::endl;
		std::cout << "######################################" << std::endl;
		std::cout << std::endl;

		SetValues(y0);

		int current_iter = 0;
		int iter_print = 10;
		int iter_ = iter;

		std::cout << "Performing " << iter_ << " iteration in direct substitution..." << std::endl;
		std::cout << std::endl;
		while (current_iter != iter_) {

			// Solve auxiliary units to provide checking values
			if (ordered_devices == true )
				FlowDevicesOrderedProcessing();
			else
				FlowDevicesProcessing();

			// Run the Unit to obtain outlet of each
			ReactorSequentialProcessing();	
			++current_iter;	

			if (current_iter%iter_print == 1){
				std::cout << current_iter << " iterations done..." << std::endl;
			}

		}
		std::cout << std::endl;
		std::cout << "Iterations completed, first guess values generated!" << std::endl;
		std::cout << "ODE and NLS subroutine starting..." << std::endl;
		std::cout << std::endl;
		// Import current network status as first guess
		int k = 0;
		for (int j = 0; j < Nunits; ++j) { // FOR EACH UNIT
			// if its a reactor
			if (Device[j]->GetTag() == "Reactor") {

				if (solid_case == true){
					for (int i = 1; i <= NSgas; ++i) { 
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].omega_gas[i]; // GAS SPECIES RESIDUAL
					}
					if (Device[j]->GetPhase() != "Gas"){
						for (int i = 1; i <= NSsolid; ++i) { 
							y0[++k] = Stream[Device[j]->GetOutlets()[0]].omega_solid[i]; // SOLID SPECIES RESIDUAL
						}
					}
					// Mass balances and T
					if (Device[j]->GetPhase() != "Gas"){
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas; // Gas
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_solid; // Solid
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].temperature; // T
					}
					else {
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas; // Gas
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].temperature; // T
					}		
				}
				else if (solid_case == false){
					for (int i = 1; i <= NSgas; ++i) { 
						y0[++k] = Stream[Device[j]->GetOutlets()[0]].omega_gas[i]; // GAS SPECIES RESIDUAL
					}

					//  T
					y0[++k] = Stream[Device[j]->GetOutlets()[0]].temperature; // T
	
				}						
			}
		}
		std::cout << std::endl;
		std::cout << "-- System ready -- " << std::endl;
		std::cout << "ODE and NLS subroutine starting..." << std::endl;
		std::cout << std::endl;


		return 0;

	};


	void ReactorNetwork::FindOptimalDeviceOrder()
	{

		std::cout << "Optimal auxliary device order possible. Finding order..." << std::endl;

		// Give some random values of reactors flows if the simulation is gas phase only (not given by setvalues) - WILL BE OVERWRITTEN BY LINEAR SYSTEM SOLVER FOR MASS FLOW
		if (solid_case == false) {
			for (unsigned int j = 0; j < Nunits; ++j) {
				if (Device[j]->GetTag() == "Reactor") {
					Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_gas = 1.;
					Stream[Device[j]->GetOutlets()[0]].mass_flow_rate_solid = 0.;
				}
			}
		}
		// Find flow devices
		std::vector<int> FlowDevicesVector;
		for (unsigned int j = 0; j < Nunits; ++j) {
			if (Device[j]->GetTag() != "Reactor") {
				FlowDevicesVector.push_back(j);
			}
		}
		

		// Fix streams associated to reactors, using initial guess is good.
		SetValues(y0);

		while ( FlowDevicesOptimalOrder.size() < FlowDevicesVector.size() ) {
		
			for (int j = 0; j < FlowDevicesVector.size(); ++j) {
				
				// Check if the device has already been studied
				int current = FlowDevicesVector[j];
		
			
				bool is_it_already_considered = std::any_of(FlowDevicesOptimalOrder.begin(), FlowDevicesOptimalOrder.end(), [&](int m) { return m == current; });
			
				// Do the stuff only if flow device is not already considered into the optimal order
				if (is_it_already_considered == false) {
					
					std::vector<double> InletValues(Device[FlowDevicesVector[j]]->GetInlets().size());
					for (int i = 0; i < Device[FlowDevicesVector[j]]->GetInlets().size(); ++i) {
						InletValues[i]=Stream[Device[FlowDevicesVector[j]]->GetInlets()[i]].mass_flow_rate_gas + Stream[Device[FlowDevicesVector[j]]->GetInlets()[i]].mass_flow_rate_solid;
					}

					// Check if one of the inlets is zero - returns true if one of the inlets has mass flow = 0
					bool is_there_a_zero = std::any_of(InletValues.begin(), InletValues.end(), [](double l) { return l == 0.; });
	
					if (is_there_a_zero == false) {
						FlowDevicesOptimalOrder.push_back(FlowDevicesVector[j]);
						Device[FlowDevicesVector[j]]->Solve(Stream);
					}
				}
			}

		}

		std::cout << "Optimal solve path found!" << std::endl;


	};
	
	// SPARSITY PATTERN ANALYZER
	void ReactorNetwork::SparsityPattern()
	{	

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "####################################" << std::endl;
		std::cout << " Finding sparsity pattern..." << std::endl;
		std::cout << "####################################" << std::endl;

		// Create pattern matrix
		Eigen::MatrixXi simplified_pattern(Nunits,Nunits);
		simplified_pattern.setZero();

		// Set diagonal to 1
		for (unsigned int i=0; i<Nunits; ++i){
			simplified_pattern(i,i) = 1;
		}

		// Build simplified pattern
		for (unsigned int i=0; i<Nunits; ++i){
			if (Device[i]->GetTag() == "Reactor"){

				bool go_deeper=false;
				int current_reactor = i;
				std::vector<int> connected_units;

				// Check the outlets of the reactor
				for (unsigned int k=0; k<Device[current_reactor]->GetOutlets().size(); ++k){
					for (unsigned int j=0; j<Nunits; ++j){
						if (FlowMatrix.coeff(j,Device[current_reactor]->GetOutlets()[k]) > 0){
							if (Device[j]->GetTag() == "Reactor"){
								simplified_pattern(j,i) = 1;
							}
							else {
								connected_units.push_back(j);
								go_deeper = true;
							}
						}
					}
				}
				
				// If the outlet goes through auxiliary units, trace the stream to the connected reactors
				while (go_deeper == true){	

					for (unsigned int l=0; l<connected_units.size(); ++l){	
						int current_unit = connected_units[l];
						connected_units.erase(connected_units.begin()+l);	
						for (unsigned int k=0; k<Device[current_unit]->GetOutlets().size(); ++k){
							for (unsigned int j=0; j<Nunits; ++j){
								if (FlowMatrix.coeff(j,Device[current_unit]->GetOutlets()[k]) > 0){
									if (Device[j]->GetTag() == "Reactor"){
										simplified_pattern(j,i) = 1;
									}
									else {
										connected_units.push_back(j);									
									}
								}
							}
						}
					}			

					// Stopping criteria
					if (connected_units.size() < 1)
						go_deeper = false;
		
				}	
				
			}
		}
		

		// Cut out auxiliary units columns and rows
		Eigen::MatrixXi simplified_pattern_cut(Nreactors,Nreactors);
		simplified_pattern_cut = simplified_pattern.topLeftCorner(Nreactors, Nreactors);
		std::cout << std::endl;
		std::cout << "Collapsed pattern of the Jacobian!" << std::endl;
		std::cout << simplified_pattern_cut << std::endl;


		// UP UNTIL THIS EVERYTHING IS FINE

		// Initialize expanded sparsity pattern
		Eigen::MatrixXi sparsity_pattern(neq,neq);
		sparsity_pattern.setIdentity();

		// Useful auxiliary counters to build up a proper pattern
		unsigned int current_rows;
		unsigned int current_columns;
		unsigned int current_block_rows;
		unsigned int current_block_columns;

		// Expand pattern cells into NSxNS blocks for species equations and add eventual mass and energy balances
		current_rows = 0;
		for (unsigned int i=0; i<Nreactors; ++i){
			current_block_rows = Device[i]->GetNumberOfEquations();
			current_columns = 0.;
			for (unsigned int j=0; j<Nreactors; ++j){		
				current_block_columns = Device[j]->GetNumberOfEquations();
				if (simplified_pattern_cut.coeff(i,j)>0){ // means i is influenced by j
					sparsity_pattern.block(current_rows, current_columns, current_block_rows, current_block_columns) = Eigen::MatrixXi::Constant(current_block_rows, current_block_columns, 1);
				}		
				current_columns += current_block_columns;
			}
			current_rows += current_block_rows;
		}

		std::cout << std::endl;
		std::cout << "The pattern has been exploded to match the real system!" << std::endl;
		std::cout << "Dimensions " << sparsity_pattern.rows() << "x" << sparsity_pattern.cols() << " as rows*columns." << std::endl;
		std::cout << std::endl;


		// Set non-zero coordinates in the vectors for the solver

		for (unsigned int i=0; i<neq; ++i){
			for (unsigned int j=0; j<neq; ++j){
				if (sparsity_pattern(i,j) > 0){
					rows_sparsity_pattern.push_back(i);
					cols_sparsity_pattern.push_back(j);

				}
			}
		}

		double percentage_of_fillage = 100. *rows_sparsity_pattern.size() / (neq*neq);
		std::cout << "Non-zero elements: " << rows_sparsity_pattern.size() << std::endl;
		std::cout << "Density of the system: " << 100. *rows_sparsity_pattern.size() / (neq*neq) << "%"<< std::endl;
		std::cout << std::endl;
		std::cout << "Rows and columns indexes of non-zero elements have been saved!" << std::endl;
		std::cout << std::endl;

		if (percentage_of_fillage < 15. || neq>8000){
			use_sparse_solvers = true;
			std::cout << "This conditions benefit from using dedicated solvers. All used kernels are set to sparse for the global system." << std::endl;
		}
		else {
			use_sparse_solvers = false;
			std::cout << "It is not worth to use a sparse solver in these conditions! All systems will be solved as dense." << std::endl;
		}

	};

	// GET RESULTS
	void ReactorNetwork::GetResults(OpenSMOKE::OpenSMOKEVectorDouble &solution, std::vector<NetSMOKE::StreamInfo> &StreamFin) // SERVIREBBE A FAR RIFARE I CALCOLI CON IL VETTORE SOLUZIONE PER OTTENERE I RISULTATI FINALI E SALVARLI, PROBABILMENTE SKIPPABILE
	{

		solution = yf;

		for (int j = 0; j < Nstreams; ++j) {
			StreamFin[j] = Stream[j];
		}

	};

	void ReactorNetwork::PrintStreamOnVideo(int j){

		// Prepare name vector
		const std::vector<std::string> NameOfSpecies(thermodynamicsMapXMLinternal.NamesOfSpecies());
		// Stream info display
			std::cout << "==============================="<<std::endl;
			std::cout << "STREAM " << Stream[j].name << std::endl;
			std::cout << "Phase : " << Stream[j].phase << std::endl;
			std::cout << "T : " << Stream[j].temperature << std::endl;
			std::cout << "P : " << Stream[j].pressure << std::endl;
			std::cout << "gas mass flow : " << Stream[j].mass_flow_rate_gas << std::endl;
			std::cout << "solid mass flow : " << Stream[j].mass_flow_rate_solid << std::endl;
			std::cout << std::endl;


	};

	void ReactorNetwork::PrintOutletStreams(){
		Eigen::VectorXi FlowMatrixColumn_temporary;
		for (int i = 0; i<Nstreams; ++i){
			FlowMatrixColumn_temporary.operator=(FlowMatrix.col(i));
			if (FlowMatrixColumn_temporary.maxCoeff()<=0.){
				ReactorNetwork::PrintStreamOnVideo(i);
			}
		}
	};

	void ReactorNetwork::PrintFinalReactorsStatus(boost::filesystem::path output_folder)
	{
		std::cout << "Generating summary files for each reactor. This might take some time..." << std::endl;
		for (int i = 0; i < Nunits; ++i){
			if (Device[i]->GetTag() == "Reactor"){
				Device[i]->PrintStatus(output_folder, Stream);
			}
		}
	};

	int ReactorNetwork::PrintODE(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		ODEiter++;
		{
			// Video output
			if (ODEiter%n_step_video_ODE_ == 1 || n_step_video_ODE_ == 1)
			{
				counter_file_video_ODE_++;
				if (counter_file_video_ODE_%100 == 1)
				{
					
					std::cout << std::endl;
					std::cout << std::setw(7)  << std::left << "#Step";
					std::cout << std::setw(30) << std::left << "Average of residuals";
					std::cout << std::setw(30) << std::left << "Norm of residuals";
					
				}
				
				std::cout << std::endl;
				std::cout << std::fixed << std::setw(7) << std::left << ODEiter;
				std::cout << std::setw(30) << std::left << std::fixed << std::setprecision(7) << dy.SumElements()/dy.Size();
				std::cout << std::setw(30) << std::left << std::fixed << std::setprecision(7) << dy.Norm1();
				std::cout << std::endl;

				// Perform backup
				WriteBackup();
				
			}
			// Track species
			if (tracking == true) {
					SetMWsAndMoleFractions();
					PrintTrack(t);
			}
		}

		return 0;
	}

	int ReactorNetwork::PrintNLS(const OpenSMOKE::OpenSMOKEVectorDouble& dy)
	{
		NLSiter++;
		{
			// Video output
			if (NLSiter%n_step_video_NLS_ == 1 || n_step_video_NLS_ == 1)
			{
				counter_file_video_NLS_++;
				if (counter_file_video_NLS_%100 == 1)
				{
					std::cout << std::endl;
					std::cout << std::setw(7)  << std::left << "#Step";
					std::cout << std::setw(30) << std::left << "Average of residuals";
					std::cout << std::setw(30) << std::left << "Norm of residuals";
				}
				std::cout << std::endl;
				std::cout << std::fixed << std::setw(7) << std::left << NLSiter;
				std::cout << std::setw(30) << std::left << std::fixed << std::setprecision(7) << dy.SumElements()/dy.Size();
				std::cout << std::setw(30) << std::left << std::fixed << std::setprecision(7) << dy.Norm1();
				std::cout << std::endl;

				// Perform backup
				WriteBackup();
				
			}
		}

		return 0;
	}	

	void ReactorNetwork::WriteBackup()
	{

		// Create binary file of system snapshot aka the vector of the current unknowns values
		{
			boost::filesystem::path backup_file = backup_folder / "Backup.bin";
			std::string backup_file_string = backup_file.generic_string();
			std::ofstream fout(backup_file_string, std::ios::binary);
			fout.write(reinterpret_cast<const char*>(&yf[1]), yf.Size() * sizeof(double));
			fout.close();
		}

		// Rename the input backup made by the constructor to match time and date of the binary file
		{

			// Get older input backup filename, which we need to rename
			boost::filesystem::directory_iterator it(backup_folder);
			boost::filesystem::directory_iterator endit;
			boost::filesystem::path old_backup_input;
			while (it != endit) {
				if (it->path().extension() == ".dic") {
					old_backup_input = it->path();
				}
				++it;
			}
			
			// Rename the file with current date and time
			std::string backup_input_name_date = "input_backup_time_" + OpenSMOKE::GetRawCurrentTime() + "_date_" + OpenSMOKE::GetCurrentDate() + ".dic";
			const boost::filesystem::path backup_input_file = backup_folder / backup_input_name_date;
			boost::filesystem::rename(old_backup_input, backup_input_file);
		}

	}

	void ReactorNetwork::InitializeFromBackUp(boost::filesystem::path backup_file) {	
		// Restores the last saved 
		std::string backup_file_string = backup_file.generic_string();
		std::ifstream fin(backup_file_string, std::ios::binary);
		fin.read(reinterpret_cast<char*>(&y0[1]), y0.Size() * sizeof(double));
		fin.close();
	}

	void ReactorNetwork::CloseTrackFile() {
		// Self explanatory
		fOutTrack.close();
	}



} /* END NAMESPACE NETSMOKE*/
