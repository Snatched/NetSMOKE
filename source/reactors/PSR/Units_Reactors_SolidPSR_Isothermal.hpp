

namespace NodusSMOKE
{

	Units_Reactors_SolidPSR_Isothermal::Units_Reactors_SolidPSR_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapGasXML, OpenSMOKE::KineticsMap_CHEMKIN *kineticsMapGasXML, OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN *thermodynamicsMapSolidXML, OpenSMOKE::KineticsMap_Solid_CHEMKIN *kineticsMapSolidXML, double rho_solid, NodusSMOKE::UnitInfo unit_data) :
		thermodynamicsMapGas(*thermodynamicsMapGasXML),
		kineticsMapGas(*kineticsMapGasXML),
		thermodynamicsMapSolid(*thermodynamicsMapSolidXML),
		kineticsMapSolid(*kineticsMapSolidXML),
		Units_Reactors_PSR(unit_data)
	{	

		// Get temperature
		temperature_ = unit_data.temperature;
		NCgas = thermodynamicsMapSolid.number_of_gas_species();
		NCsolid = thermodynamicsMapSolid.number_of_solid_species();

		OpenSMOKE::ChangeDimensions(NCgas, &Omega_previous_gas, true);
		OpenSMOKE::ChangeDimensions(NCsolid, &Omega_previous_solid, true);
		n_calls = 0;

		// Allocate memory 
		OpenSMOKE::ChangeDimensions(thermodynamicsMapSolid.number_of_solid_species(), &StreamOut.omega_solid, true);
		OpenSMOKE::ChangeDimensions(thermodynamicsMapGas.NumberOfSpecies(), &StreamOut.omega_gas, true);

		OpenSMOKE::ChangeDimensions(NCgas, &Rgas, true);
		OpenSMOKE::ChangeDimensions(NCsolid, &Rgas_from_solid, true);
		OpenSMOKE::ChangeDimensions(NCsolid, &Rsolid, true);

		rho_solid_ = rho_solid;

		ode_parameters = new OpenSMOKE::ODE_Parameters();
		spsr_options = new OpenSMOKE::SolidPerfectlyStirredReactor_Options();
		spsr_options->SetVerboseVideo(false);
		spsr_options->SetVerboseASCIIFile(false);
		spsr_options->SetVerboseOutput(false);
		spsr_options->SetVerboseXMLFile(false);

	};


	int Units_Reactors_SolidPSR_Isothermal::Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {
		
		StreamIn = streams_data_structure[inlets_[0]];

		if (n_calls < 3) {
			Omega_previous_gas =  StreamIn.omega_gas;
			Omega_previous_solid =  StreamIn.omega_solid;
		}
		
		if (StreamIn.omega_gas.SumElements() > 1.1 || StreamIn.omega_solid.SumElements() > 1.1) {
			std::cout << "INLETS" << std::endl;
			for (unsigned int i = 1; i <= NCgas; ++i) {
				std::cout << thermodynamicsMapGas.NamesOfSpecies()[i - 1] << " " << StreamIn.omega_gas[i] << std::endl;
			}
			for (unsigned int i = 1; i <= NCsolid; ++i) {
				std::cout << thermodynamicsMapSolid.NamesOfSpecies()[NCgas + i - 1] << " " << StreamIn.omega_solid[i] << std::endl;
			}
			std::cout << std::endl;
			OpenSMOKE::FatalErrorMessage("Bois we just put unfeasible conditions in reactor " + name_ + ". Sorry for that!");
		}
		
		// Construct PSR
		OpenSMOKE::SolidPerfectlyStirredReactor_Isothermal_ConstantPressure solidpsr(thermodynamicsMapGas, kineticsMapGas,
			thermodynamicsMapSolid, kineticsMapSolid, *ode_parameters, *spsr_options,
			temperature_, StreamIn.pressure, Omega_previous_gas, Omega_previous_solid,
			StreamIn.temperature, StreamIn.pressure,  StreamIn.omega_gas,  StreamIn.omega_solid,
			StreamIn.mass_flow_rate_gas,  StreamIn.mass_flow_rate_solid, rho_solid_,
			residence_time_, -1.);


		int status = solidpsr.Solve(1.e8);

		
		if (status == -1) {
			std::cout << "INLETS" << std::endl;
			for (unsigned int i = 1; i <= NCgas; ++i) {
				std::cout << thermodynamicsMapGas.NamesOfSpecies()[i - 1] << " " << StreamIn.omega_gas[i] << std::endl;
			}
			for (unsigned int i = 1; i <= NCsolid; ++i) {
				std::cout << thermodynamicsMapSolid.NamesOfSpecies()[NCgas + i - 1] << " " << StreamIn.omega_solid[i] << std::endl;
			}
			OpenSMOKE::FatalErrorMessage("Reactor " + name_ + " couldn't finish integration");
		}

		// Local outlet stream data storage
		solidpsr.GetFinalStatus(StreamOut.temperature, StreamOut.pressure, StreamOut.mass_flow_rate_gas, StreamOut.mass_flow_rate_solid, StreamOut.omega_gas, StreamOut.omega_solid);
		// Recover stuff to calculate residuals
		solidpsr.GetFinalResults( rhoGas, mass_gas, mass_solid, Rgas, Rgas_from_solid, Rsolid);
		
		// Import into "global" structure
		streams_data_structure[outlets_[0]].phase = "Mix";
		StreamOut.name = streams_data_structure[outlets_[0]].name;

		Omega_previous_gas = StreamOut.omega_gas;	
		Omega_previous_solid = StreamOut.omega_solid;		
		n_calls++;
		
		return 0;
	};

	/* Get Residuals */
	int Units_Reactors_SolidPSR_Isothermal::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){
		
		StreamOut_star = streams_data_structure[outlets_[0]];

		// PHENOMENOLOGICAL EQUATIONS	
		Vsolid = mass_solid/rho_solid_;
		
		// Gas species composition
		for (unsigned int i = 1; i <= NCgas; ++i) {
				residuals[i]=(StreamIn.mass_flow_rate_gas*StreamIn.omega_gas[i] - StreamOut_star.mass_flow_rate_gas*StreamOut_star.omega_gas[i])/mass_gas										
							+ thermodynamicsMapGas.MW(i - 1)*(Rgas[i] / rhoGas + Rgas_from_solid[i]*Vsolid/mass_gas);

		}
		

		// Solid species composition
		for (unsigned int i = 1; i <= NCsolid; ++i){	
				residuals[i + NCgas] = (StreamIn.mass_flow_rate_solid*StreamIn.omega_solid[i] - StreamOut_star.mass_flow_rate_solid*StreamOut_star.omega_solid[i]) / mass_solid
									+ thermodynamicsMapSolid.MW(NCgas + i - 1)*Rsolid[i] / rho_solid_;	
		}
				
		// Net gas production in mass		
		for (unsigned int i = 1; i<=NCgas; ++i){
			Rgas_from_solid[i] = Rgas_from_solid[i]*thermodynamicsMapGas.MW(i - 1);
		}
		double net_gas_production = Rgas_from_solid.SumElements() + Rgas.SumElements();
		
		// Net solid production in mass
		for (unsigned int i = 1; i<=NCsolid; ++i){
			Rsolid[i] = Rsolid[i]*thermodynamicsMapSolid.MW(NCgas + i - 1);
		}
		double net_solid_production = Rsolid.SumElements();

		// Gas flow balance
		residuals[NCgas+NCsolid+1] = StreamIn.mass_flow_rate_gas - StreamOut_star.mass_flow_rate_gas + net_gas_production*Vsolid;

		// Solid flow balance
		residuals[NCgas+NCsolid+2] = StreamIn.mass_flow_rate_solid - StreamOut_star.mass_flow_rate_solid + net_solid_production*Vsolid;

		// Energy balance
		residuals[NCgas+NCsolid+3] = StreamIn.temperature - StreamOut_star.temperature;

		return 0;
	};

	/* Solve sequentially */
	int Units_Reactors_SolidPSR_Isothermal::NonIterativeSolve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure)
	{
		
		Units_Reactors_SolidPSR_Isothermal::Solve(streams_data_structure);
		streams_data_structure[outlets_[0]].mass_flow_rate_gas = StreamOut.mass_flow_rate_gas;
		streams_data_structure[outlets_[0]].mass_flow_rate_solid = StreamOut.mass_flow_rate_solid;
		streams_data_structure[outlets_[0]].temperature = StreamOut.temperature;
		streams_data_structure[outlets_[0]].omega_gas = StreamOut.omega_gas;
		streams_data_structure[outlets_[0]].omega_solid = StreamOut.omega_solid;
		streams_data_structure[outlets_[0]].pressure = StreamOut.pressure;

		return 0;
	};

	/* RTD */
	int Units_Reactors_SolidPSR_Isothermal::RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {
		Units_Reactors_SolidPSR_Isothermal::GetResiduals(residuals, streams_data_structure);
		return 0;
	}

	/* Print reactor status */
	
	void Units_Reactors_SolidPSR_Isothermal::PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){
		
		boost::filesystem::path output = output_folder / name_;
		OpenSMOKE::CreateDirectory(output);

		// Hide printing - fucking eats cpu
		spsr_options->SetOutputPath(output);
		spsr_options->SetVerboseVideo(false);
		spsr_options->SetVerboseASCIIFile(true);
		spsr_options->SetVerboseOutput(true);
		spsr_options->SetVerboseXMLFile(false);

		StreamIn = streams_data_structure[inlets_[0]];

		// Construct PSR
		OpenSMOKE::SolidPerfectlyStirredReactor_Isothermal_ConstantPressure solidpsr(thermodynamicsMapGas, kineticsMapGas,
			thermodynamicsMapSolid, kineticsMapSolid,*ode_parameters, *spsr_options,
			temperature_, StreamIn.pressure, StreamIn.omega_gas, StreamIn.omega_solid,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas, StreamIn.omega_solid,
			StreamIn.mass_flow_rate_gas, StreamIn.mass_flow_rate_solid,rho_solid_,
			residence_time_, -1.);
			
		// Solve PSR
		solidpsr.Solve(1.e8);
		
	};

	Units_Reactors_SolidPSR_Isothermal::~Units_Reactors_SolidPSR_Isothermal()
	{
		delete spsr_options;
	};


	

}; // END OF NAMESPACE
