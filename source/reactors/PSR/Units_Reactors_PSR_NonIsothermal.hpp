
namespace NodusSMOKE{

	Units_Reactors_PSR_NonIsothermal::Units_Reactors_PSR_NonIsothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapXML, NodusSMOKE::UnitInfo unit_data) :
		thermodynamicsMap(*thermodynamicsMapXML),
		kineticsMap(*kineticsMapXML),
		Units_Reactors_PSR(unit_data)

	{

		// Get heat exchange
		UA_ = unit_data.UA;

		// Allocate memory
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(),&StreamOut.omega_gas,true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(), &Omega_previous_gas, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMap.NumberOfSpecies(), &Rgas, true);
		n_calls = 0;

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


	int Units_Reactors_PSR_NonIsothermal::Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {
		
		StreamIn = streams_data_structure[inlets_[0]];

		if (n_calls <= 5) {
			Omega_previous_gas = StreamIn.omega_gas;
		}

		// Construct PSR
		OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure psr(
									thermodynamicsMap,
									kineticsMap,
									*ode_parameters,
									*psr_options,
									*onTheFlyROPA,
									*on_the_fly_post_processing,
									*polimi_soot,
									StreamIn.temperature, StreamIn.pressure, Omega_previous_gas,
									StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
									residence_time_, -1, StreamIn.mass_flow_rate_gas,
									sqrt(UA_),sqrt(UA_), StreamIn.temperature);


		// Solve PSR
		psr.Solve(1.e8);

		// Local data storing
		psr.GetFromationRates_HeatRelease_Density( Rgas, Qr, StreamOut.rho_gas);
		psr.GetFinalStatus(StreamOut.temperature, StreamOut.pressure, StreamOut.omega_gas);
		
		StreamOut.mass_flow_rate_gas = StreamIn.mass_flow_rate_gas;
		Omega_previous_gas = StreamOut.omega_gas;

		// Set outlet phase
		streams_data_structure[outlets_[0]].phase = "Gas";
		
		n_calls++;

		return 0;
	};

	/* Get Residuals */
	int Units_Reactors_PSR_NonIsothermal::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){

		StreamOut_star = streams_data_structure[outlets_[0]];

		for (unsigned int i = 1; i <= thermodynamicsMap.NumberOfSpecies(); ++i){
			residuals[i] = (StreamIn.omega_gas[i] - StreamOut_star.omega_gas[i])/residence_time_ + thermodynamicsMap.MW(i-1)*Rgas[i]/StreamOut.rho_gas;
		}
		
		residuals[thermodynamicsMap.NumberOfSpecies()+1] = StreamIn.mass_flow_rate_gas - StreamOut_star.mass_flow_rate_gas;
		
		
		thermodynamicsMap.SetPressure(pressure_);
		//Hin
		{
			thermodynamicsMap.SetTemperature(StreamIn.temperature);
			OpenSMOKE::OpenSMOKEVectorDouble temp_x_gas(thermodynamicsMap.NumberOfSpecies());
			thermodynamicsMap.MoleFractions_From_MassFractions(temp_x_gas.GetHandle(),StreamIn.MW_gas,StreamIn.omega_gas.GetHandle());
			Hin = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(temp_x_gas.GetHandle());
			StreamIn.x_gas = temp_x_gas;
		}
		//Hout
		{
			thermodynamicsMap.SetTemperature(StreamOut_star.temperature);
			OpenSMOKE::OpenSMOKEVectorDouble temp_x_gas(thermodynamicsMap.NumberOfSpecies());
			thermodynamicsMap.MoleFractions_From_MassFractions(temp_x_gas.GetHandle(),StreamOut_star.MW_gas,StreamOut_star.omega_gas.GetHandle());
			// Outlet Cp
			const double CpMixMolar = thermodynamicsMap.cpMolar_Mixture_From_MoleFractions(temp_x_gas.GetHandle());
			CpMix = CpMixMolar /StreamOut_star.MW_gas;
			// Outlet density and volume
			double cTot_ = pressure_/(PhysicalConstants::R_J_kmol *StreamOut_star.temperature);
			StreamOut_star.rho_gas = cTot_ * StreamOut_star.MW_gas;
			Vgas = StreamOut_star.mass_flow_rate_gas*residence_time_/StreamOut_star.rho_gas;
			// Exchanged heat
			Qe = UA_*(StreamIn.temperature - StreamOut_star.temperature);		// [W]
			// H
			HStarI = OpenSMOKE::Dot(StreamIn.x_gas.Size(), StreamIn.x_gas.GetHandle(), thermodynamicsMap.Species_H_over_RT().data())*(PhysicalConstants::R_J_kmol*StreamOut_star.temperature);	// [J/kmol]
			// HStarI = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(temp_x_gas.GetHandle());
		}
	
		// Temperature
		residuals[thermodynamicsMap.NumberOfSpecies()+2] = ( (Hin-HStarI)/StreamIn.MW_gas/residence_time_ + Qr/StreamOut_star.rho_gas + Qe/(StreamOut_star.rho_gas*Vgas)) / CpMix;

		return 0;
	};

	/* Solve sequentially */
	int Units_Reactors_PSR_NonIsothermal::NonIterativeSolve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure)
	{
		Units_Reactors_PSR_NonIsothermal::Solve(streams_data_structure);
		streams_data_structure[outlets_[0]].mass_flow_rate_gas = StreamOut.mass_flow_rate_gas;
		streams_data_structure[outlets_[0]].temperature = StreamOut.temperature;
		streams_data_structure[outlets_[0]].omega_gas = StreamOut.omega_gas;
		streams_data_structure[outlets_[0]].pressure = StreamOut.pressure;
		return 0;
	};

	/* Solve for RTD */
	int Units_Reactors_PSR_NonIsothermal::RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {
		Units_Reactors_PSR_NonIsothermal::GetResiduals(residuals, streams_data_structure);
		return 0;
	}

	/* Print reactor status */
	void Units_Reactors_PSR_NonIsothermal::PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){
		
		boost::filesystem::path output = output_folder / name_;
		OpenSMOKE::CreateDirectory(output);

		// Hide printing - fucking eats cpu
		psr_options->SetOutputPath(output);
		psr_options->SetVerboseVideo(false);
		psr_options->SetVerboseASCIIFile(true);
		psr_options->SetVerboseOutput(true);
		psr_options->SetVerboseXMLFile(false);


		StreamIn = streams_data_structure[inlets_[0]];

		// Construct PSR
		OpenSMOKE::PerfectlyStirredReactor_NonIsothermal_ConstantPressure psr(
									thermodynamicsMap,
									kineticsMap,
									*ode_parameters,
									*psr_options,
									*onTheFlyROPA,
									*on_the_fly_post_processing,
									*polimi_soot,
									StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
									StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
									residence_time_, -1, StreamIn.mass_flow_rate_gas,
									UA_,
									UA_, StreamIn.temperature);

		// Solve PSR
		psr.Solve(1.e8);
	};

	Units_Reactors_PSR_NonIsothermal::~Units_Reactors_PSR_NonIsothermal()
	{
		delete psr_options;
		delete onTheFlyROPA;
		delete on_the_fly_post_processing;
		delete polimi_soot;
		delete ode_parameters;
	};

}
