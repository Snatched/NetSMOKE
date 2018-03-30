

namespace NetSMOKE
{

	Units_Reactors_PSR_Isothermal::Units_Reactors_PSR_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML, NetSMOKE::UnitInfo unit_data) :
		thermodynamicsMap(*thermodynamicsMapXML),
		kineticsMap(*kineticsMapXML),
		Units_Reactors_PSR(unit_data)

	{	
		
		// Import unit info
		temperature_ = unit_data.temperature;

		// Allocate memory
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(),&StreamOut.omega_gas,true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(), &Omega_previous_gas, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(), &Rgas, true);
		n_calls = 0;

		// CONSTRUCT MANDATORY OBJECTS THAT CAN BE RECYCLED
		// Options
		psr_options = new OpenSMOKE::PerfectlyStirredReactor_Options();

		// Hide printing - fucking eats cpu
		psr_options->SetVerboseVideo(false);
		psr_options->SetVerboseASCIIFile(false);
		psr_options->SetVerboseOutput(false);
		psr_options->SetVerboseXMLFile(false);

		// ODE Parameters
		ode_parameters = new OpenSMOKE::ODE_Parameters();

		// On the fly ROPA
		onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(thermodynamicsMap, kineticsMap);

		// On the fly PostProcessing
		on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(thermodynamicsMap, kineticsMap, psr_options->output_path());

		// Polimi soot
		polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(&thermodynamicsMap);

	};


	int Units_Reactors_PSR_Isothermal::Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {

		// Import inlet streams_data_structure from main datastructure
		StreamIn = streams_data_structure[inlets_[0]];
		
		if (n_calls <= 10) {
			Omega_previous_gas = StreamIn.omega_gas;
		}

		if (StreamIn.omega_gas.SumElements() > 1.1) {
			std::cout << "INLETS" << std::endl;
			for (unsigned int i = 1; i <= NCgas; ++i) {
				std::cout << thermodynamicsMap.NamesOfSpecies()[i - 1] << " " << StreamIn.omega_gas[i] << std::endl;
			}
			std::cout << std::endl;
			OpenSMOKE::FatalErrorMessage("Bois we just put unfeasible conditions in reactor " + name_ + ". Sorry for that!");
		}
		
		// Create reactor object
		OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure psr(thermodynamicsMap, kineticsMap,
			*ode_parameters, *psr_options, *onTheFlyROPA, *on_the_fly_post_processing, *polimi_soot,
			temperature_, StreamIn.pressure, Omega_previous_gas,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
			residence_time_, -1., StreamIn.mass_flow_rate_gas);
			
		// Solve PSR
		psr.Solve(1.e8);

		// Data storing
		psr.GetFinalResults(StreamOut.omega_gas, Rgas, StreamOut.rho_gas);
		psr.GetFinalStatus(StreamOut.temperature, StreamOut.pressure, StreamOut.omega_gas);
		StreamOut.mass_flow_rate_gas = StreamIn.mass_flow_rate_gas;

		// Set outlet phase
		streams_data_structure[outlets_[0]].phase = "Gas";
		
		Omega_previous_gas = StreamOut.omega_gas;

		n_calls++;

		return 0;
	};

	/* Get Residuals */
	int Units_Reactors_PSR_Isothermal::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		
		StreamOut_star = streams_data_structure[outlets_[0]];

		// Gas species balances
		for (unsigned int i = 1; i <= thermodynamicsMap.NumberOfSpecies(); ++i){
			residuals[i] = (StreamIn.omega_gas[i] - StreamOut_star.omega_gas[i])/residence_time_ + thermodynamicsMap.MW(i-1)*Rgas[i]/StreamOut.rho_gas;
		}
		// Mass balance
		residuals[thermodynamicsMap.NumberOfSpecies()+1] = StreamIn.mass_flow_rate_gas - StreamOut_star.mass_flow_rate_gas;
		// Energy balance
		residuals[thermodynamicsMap.NumberOfSpecies()+2] = StreamIn.temperature - StreamOut_star.temperature;
		
		return 0;
	};

	/* Solve sequentially */
	int Units_Reactors_PSR_Isothermal::NonIterativeSolve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure)
	{
		n_calls = 0;
		Units_Reactors_PSR_Isothermal::Solve(streams_data_structure);
		streams_data_structure[outlets_[0]].mass_flow_rate_gas = StreamOut.mass_flow_rate_gas;
		streams_data_structure[outlets_[0]].temperature = StreamOut.temperature;
		streams_data_structure[outlets_[0]].omega_gas = StreamOut.omega_gas;
		streams_data_structure[outlets_[0]].pressure = StreamOut.pressure;
		return 0;
	};

	/* Solve for RTD */
	int Units_Reactors_PSR_Isothermal::RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {
		Units_Reactors_PSR_Isothermal::GetResiduals(residuals, streams_data_structure);
		return 0;
	}

	/* Print reactor status */
	void Units_Reactors_PSR_Isothermal::PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		
		boost::filesystem::path output = output_folder / name_;
		OpenSMOKE::CreateDirectory(output);

		// Hide printing - fucking eats cpu
		psr_options->SetOutputPath(output);
		psr_options->SetVerboseVideo(false);
		psr_options->SetVerboseASCIIFile(true);
		psr_options->SetVerboseOutput(true);
		psr_options->SetVerboseXMLFile(false);

		// Import inlet streams_data_structure from main datastructure
		StreamIn = streams_data_structure[inlets_[0]];

		// PSR object
		OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure psr(thermodynamicsMap, kineticsMap,
			*ode_parameters, *psr_options, *onTheFlyROPA, *on_the_fly_post_processing, *polimi_soot,
			temperature_, StreamIn.pressure, StreamIn.omega_gas,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
			residence_time_, -1., StreamIn.mass_flow_rate_gas);
			
		// Solve PSR
		psr.Solve(1.e8);
		
	};

	Units_Reactors_PSR_Isothermal::~Units_Reactors_PSR_Isothermal()
	{
		delete psr_options;
		delete onTheFlyROPA;
		delete on_the_fly_post_processing;
		delete polimi_soot;
		delete ode_parameters;
	};

}
