
namespace NodusSMOKE

{

	// Constructor for gas phase case
	Units_Mixer_Isothermal::Units_Mixer_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, NodusSMOKE::UnitInfo unit_data) :
		thermodynamicsMap(thermodynamicsMapXML),
		Units(unit_data)
	{

		// Set up stuff
		n_gas_species = thermodynamicsMap->NumberOfSpecies();
		n_inlets = inlets_.size();
		n_solid_species = 0;
		temperature_ = unit_data.temperature;

		// Memory allocation
		OpenSMOKE::ChangeDimensions(n_gas_species, &GasInletComponentMassFlow, true); // used to compute final omega
		OpenSMOKE::ChangeDimensions(n_gas_species, &omega_final_gas, true);
		contains_gas = true;
		contains_solid = false;


	};

	// Alternative construction for solid species case
	Units_Mixer_Isothermal::Units_Mixer_Isothermal(OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN *thermodynamicsMapXML, NodusSMOKE::UnitInfo unit_data) :
		thermodynamicsMapS(thermodynamicsMapXML),
		Units(unit_data)
	{

		// Set up stuff
		n_gas_species = thermodynamicsMapS->number_of_gas_species();
		n_solid_species = thermodynamicsMapS->number_of_solid_species();
		n_inlets = inlets_.size();
		temperature_ = unit_data.temperature;

		// Memory allocation
		OpenSMOKE::ChangeDimensions(n_inlets, &GasInletComponentMassFlow, true); // used to compute final omega
		OpenSMOKE::ChangeDimensions(n_gas_species, &omega_final_gas, true);
		OpenSMOKE::ChangeDimensions(n_inlets, &SolidInletComponentMassFlow, true); // used to compute final omega
		OpenSMOKE::ChangeDimensions(n_solid_species, &omega_final_solid, true);
		contains_gas = false;
		contains_solid = false;


	};

	// SOLVE THE MIXER
	int Units_Mixer_Isothermal::Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure)
	{

		// Mixture Gas composition
		for (int i = 1; i <= n_gas_species; ++i){
			for (int j = 0; j<n_inlets; ++j) {
				GasInletComponentMassFlow[j+1] = streams_data_structure[inlets_[j]].omega_gas[i]*streams_data_structure[inlets_[j]].mass_flow_rate_gas;
			}
			omega_final_gas[i] = GasInletComponentMassFlow.SumElements();		
		}
		mass_flow_rate_gas_ = omega_final_gas.SumElements();
		if (mass_flow_rate_gas_ > 0.){ // Allows to avoid dividing by zero
			contains_gas = true;
			omega_final_gas = omega_final_gas.operator/=(mass_flow_rate_gas_);
		}
		else {
			contains_gas = false;
		}

		// Mixture Solid composition
		for (int i = 1; i <= n_solid_species; ++i){
			for (int j = 0; j<n_inlets; ++j) {
				SolidInletComponentMassFlow[j+1] = streams_data_structure[inlets_[j]].omega_solid[i]*streams_data_structure[inlets_[j]].mass_flow_rate_solid;
			}
			omega_final_solid[i] = SolidInletComponentMassFlow.SumElements();
		}
		mass_flow_rate_solid_ = omega_final_solid.SumElements();
		if (mass_flow_rate_solid_ > 0.) { // Allows to avoid dividinig by zero
			contains_solid = true;
			omega_final_solid = omega_final_solid.operator/=(mass_flow_rate_solid_);
		}
		else {
			contains_solid = false;
		}
		
		// Set phase
		if (contains_solid == true && contains_gas == false)
			StreamOut.phase = "Solid";
		else if (contains_solid == false && contains_gas == true)
			StreamOut.phase = "Gas";
		else if (contains_solid == true && contains_gas == true)
			StreamOut.phase = "Mix";
		else
			OpenSMOKE::FatalErrorMessage("Mixer " + name_ + " seems to output nothing.");
		

		// Save on outlet streams_data_structure
		StreamOut.mass_flow_rate_gas = mass_flow_rate_gas_;
		StreamOut.omega_gas = omega_final_gas;
		StreamOut.mass_flow_rate_solid = mass_flow_rate_solid_;
		StreamOut.omega_solid = omega_final_solid;
		StreamOut.temperature = temperature_;
		StreamOut.pressure = pressure_;
		StreamOut.name = streams_data_structure[outlets_[0]].name;
		
		// Import over common data structure
		streams_data_structure[outlets_[0]] = StreamOut;

		return 0;
	};

	/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
	int Units_Mixer_Isothermal::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure)
	{
		OpenSMOKE::FatalErrorMessage(" Residuals have been requested for a trivial unit (MIXER). Please check your input file or contact Matteo Mensi at matteo.mensi@mail.polimi.it");
		return 0;
	};

	void Units_Mixer_Isothermal::PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){
		OpenSMOKE::FatalErrorMessage("No status printer for mixers yet");
	}

}
