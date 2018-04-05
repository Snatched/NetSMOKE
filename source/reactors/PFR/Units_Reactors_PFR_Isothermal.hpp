

namespace NetSMOKE
{

	Units_Reactors_PFR_Isothermal::Units_Reactors_PFR_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML, NetSMOKE::UnitInfo unit_data)
		:
		thermodynamicsMap(*thermodynamicsMapXML),
		kineticsMap(*kineticsMapXML),
		Units_Reactors_PFR(unit_data)

	{

		// Get T
		temperature_ = unit_data.temperature;

		// Allocate memory
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(),&StreamOut.omega_gas,true);


		// Create static objects
		// Options
		plugflow_options = new OpenSMOKE::PlugFlowReactor_Options();

		// Hide printing - fucking eats cpu

		plugflow_options->SetVerboseVideo(false);
		plugflow_options->SetVerboseASCIIFile(false);
		plugflow_options->SetVerboseOutput(false);
		plugflow_options->SetVerboseXMLFile(false);


		// ODE Parameters
		ode_parameters = new OpenSMOKE::ODE_Parameters();

		// On the fly ROPA
		onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(thermodynamicsMap, kineticsMap);

		// On the fly PostProcessing
		on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(thermodynamicsMap, kineticsMap, plugflow_options->output_path());

		// Polimi soot
		polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(&thermodynamicsMap);

		// RTD booleans
		inert_in = false;
		inert_out = false;
		n_calls = 0;

	};


	int Units_Reactors_PFR_Isothermal::Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {
		
		if (n_calls > 1) {
			Previous_StreamIn = StreamIn;
		}
		else {
			Previous_StreamIn = streams_data_structure[inlets_[0]];
		}
		StreamIn = streams_data_structure[inlets_[0]];
		
		bool constant_pressure = true;
		bool time_independent_variable;
		double velocity;
		if( residence_time_ > 0.){
			time_independent_variable= true;
			velocity = 1.;
		}
		else {
			time_independent_variable = false;	
			double area = 3.14*diameter_*diameter_/4.;
			const double MW = thermodynamicsMap.MolecularWeight_From_MassFractions(StreamIn.omega_gas.GetHandle());
			const double rho = pressure_ * MW / PhysicalConstants::R_J_kmol/temperature_;
			velocity = StreamIn.mass_flow_rate_gas/area/rho;
		}

		// Construct PFR
		OpenSMOKE::PlugFlowReactor_Isothermal plugflow(thermodynamicsMap, kineticsMap,
													 *ode_parameters, *plugflow_options,
													 *onTheFlyROPA, *on_the_fly_post_processing, *polimi_soot,
													 time_independent_variable, constant_pressure,
													 velocity, temperature_, StreamIn.pressure, StreamIn.omega_gas);


		// Solve PSR
		if( residence_time_ > 0.){
			plugflow.Solve(residence_time_);
		}
		else {
			plugflow.Solve(length_);		 
		}

		// Local data storing
		plugflow.GetFinalStatus(StreamOut.temperature, StreamOut.pressure, StreamOut.omega_gas);
		
		StreamOut.mass_flow_rate_gas = StreamIn.mass_flow_rate_gas;
		StreamOut.mass_flow_rate_solid = 0.;
		StreamOut.phase = "Gas";
		// Set outlet phase
		streams_data_structure[outlets_[0]].phase = "Gas";
	
		return 0;
	};

	/* Get Residuals */
	int Units_Reactors_PFR_Isothermal::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		StreamOut_star = streams_data_structure[outlets_[0]];
		// Species composition
		for (int i = 1; i <= thermodynamicsMap.NumberOfSpecies(); ++i)
		{
			residuals[i] = (StreamOut.omega_gas[i] - StreamOut_star.omega_gas[i]);
		}
		// Global mass balance
		residuals[thermodynamicsMap.NumberOfSpecies()+1] = StreamOut.mass_flow_rate_gas - StreamOut_star.mass_flow_rate_gas;
		// Energy balance
		residuals[thermodynamicsMap.NumberOfSpecies()+2] = StreamIn.temperature - StreamOut_star.temperature;
		return 0;
	};

	/* Solve for RTD */
	int Units_Reactors_PFR_Isothermal::RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {
		 
		
		if (inert_in == false) {
			for (unsigned int i = 1; i <= thermodynamicsMap.NumberOfSpecies(); ++i) {
				if (std::abs(StreamIn.omega_gas[i] - Previous_StreamIn.omega_gas[i]) > 1.e-9) {
					inert_in = true;
					t_start = t;
					break;
				}
			}
		}

		if (inert_in == true) {
			if (t >= residence_time_ + t_start) {
				StreamOut.name = streams_data_structure[outlets_[0]].name;
				streams_data_structure[outlets_[0]] = StreamOut;
			}
			else {
				StreamOut = StreamIn;
				StreamOut.name = streams_data_structure[outlets_[0]].name;
				streams_data_structure[outlets_[0]] = StreamOut;
			}
		}


		Units_Reactors_PFR_Isothermal::GetResiduals(residuals, streams_data_structure);
		
		return 0;
		
	};

	/* Solve sequentially */
	int Units_Reactors_PFR_Isothermal::NonIterativeSolve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure)
	{
		
		Units_Reactors_PFR_Isothermal::Solve(streams_data_structure);
		//Get outlet
		streams_data_structure[outlets_[0]].mass_flow_rate_gas = StreamOut.mass_flow_rate_gas;
		streams_data_structure[outlets_[0]].temperature = StreamOut.temperature;
		streams_data_structure[outlets_[0]].omega_gas = StreamOut.omega_gas;
		streams_data_structure[outlets_[0]].pressure = StreamOut.pressure;

		return 0;
	};


	/* Print reactor status */
	void Units_Reactors_PFR_Isothermal::PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		
		
		boost::filesystem::path output = output_folder / name_;
		OpenSMOKE::CreateDirectory(output);


		// Hide printing - fucking eats cpu
		plugflow_options->SetOutputPath(output);
		plugflow_options->SetVerboseVideo(false);
		plugflow_options->SetVerboseASCIIFile(true);
		plugflow_options->SetVerboseOutput(true);
		plugflow_options->SetVerboseXMLFile(false);

		bool constant_pressure = true;
		bool time_independent_variable;
		double velocity;
		if( residence_time_ > 0.){
			time_independent_variable= true;
			velocity = 1.;
		}
		else {
			time_independent_variable = false;	
			double area = 3.14*diameter_*diameter_/4.;
			const double MW = thermodynamicsMap.MolecularWeight_From_MassFractions(StreamIn.omega_gas.GetHandle());
			const double rho = pressure_ * MW / PhysicalConstants::R_J_kmol/temperature_;
			velocity = StreamIn.mass_flow_rate_gas/area/rho;
		}

		StreamIn = streams_data_structure[inlets_[0]];

		// Construct PFR
		OpenSMOKE::PlugFlowReactor_Isothermal plugflow(thermodynamicsMap, kineticsMap,
													 *ode_parameters, *plugflow_options,
													 *onTheFlyROPA, *on_the_fly_post_processing, *polimi_soot,
													 time_independent_variable, constant_pressure,
													 velocity, temperature_, StreamIn.pressure, StreamIn.omega_gas);


		// Solve PSR
		if( residence_time_ > 0.){
			plugflow.Solve(residence_time_);
		}
		else {
			plugflow.Solve(length_);		 
		}


	};

	/* Default destructor */
	Units_Reactors_PFR_Isothermal::~Units_Reactors_PFR_Isothermal()
	{
		delete plugflow_options;
		delete onTheFlyROPA;
		delete on_the_fly_post_processing;
		delete polimi_soot;
		delete ode_parameters;
	};


}
