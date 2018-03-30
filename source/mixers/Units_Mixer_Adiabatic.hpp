
namespace NetSMOKE
{
	// CONSTRUCTOR (REQUIRES INLET STREAM, SPLIT RATIO, )
	Units_Mixer_Adiabatic::Units_Mixer_Adiabatic(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, NetSMOKE::UnitInfo unit_data) :
		thermodynamicsMap(*thermodynamicsMapXML),
		Units(unit_data)
	{


		// Set up stuff

		n_species = thermodynamicsMap.NumberOfSpecies();
		n_inlets = inlets_.size();

		// Memory allocation
		OpenSMOKE::ChangeDimensions(n_inlets, &InletComponentMassFlow, true); // used to compute final omega
		OpenSMOKE::ChangeDimensions(n_species, &omega_final, true);
		OpenSMOKE::ChangeDimensions(n_species, &mole_fraction, true);
		OpenSMOKE::ChangeDimensions(n_species, &xInlet, true);
		OpenSMOKE::ChangeDimensions(n_inlets, &HInlet, true);
		OpenSMOKE::ChangeDimensions(n_inlets, &TFirstGuess_Num, true);
		OpenSMOKE::ChangeDimensions(n_inlets, &TFirstGuess_Denum, true);



	};

	// SOLVE THE MIXER
	int Units_Mixer_Adiabatic::Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure)
	{

		for ( unsigned int i = 1; i<= inlets_.size(); ++i){
			if (streams_data_structure[inlets_[i-1]].phase != "Gas"){
				std::cout << "Stream " << streams_data_structure[inlets_[i-1]].name << " is a non-gas stream entering adiabatic mixer " << name_ << "."<<std::endl;
				OpenSMOKE::FatalErrorMessage("Currently it is not possible to use adiabatic mixers with solid streams!");
			}			
		}

		// Per Species mass balance
		for (int i = 1; i <= n_species; ++i){
			for (int j = 0; j<n_inlets; ++j) {
				InletComponentMassFlow[j+1] = streams_data_structure[inlets_[j]].omega_gas[i]*streams_data_structure[inlets_[j]].mass_flow_rate_gas;
			}
			omega_final[i] = InletComponentMassFlow.SumElements();
		}
		mass_flow_rate_gas_ = omega_final.SumElements();
		omega_final = omega_final.operator/=(mass_flow_rate_gas_);

		// Get Final Mole Fractions & MW mix
		thermodynamicsMap.MoleFractions_From_MassFractions(mole_fraction.GetHandle(), mix_molecular_weight, omega_final.GetHandle());


		// Get inlet enthalpy and first guess temperature
		for (int j = 0; j < inlets_.size(); ++j) {
			// stream enthalpy
			double mw;
			thermodynamicsMap.MoleFractions_From_MassFractions(xInlet.GetHandle(), mw, streams_data_structure[inlets_[j]].omega_gas.GetHandle());
			thermodynamicsMap.SetPressure(streams_data_structure[inlets_[j]].pressure);
			thermodynamicsMap.SetTemperature(streams_data_structure[inlets_[j]].temperature);
			HInlet[j+1] = thermodynamicsMap.hMolar_Mixture_From_MoleFractions(xInlet.GetHandle())*streams_data_structure[inlets_[j]].mass_flow_rate_gas/mw;

			// first guess for T
			Cpstream = thermodynamicsMap.cpMolar_Mixture_From_MoleFractions(xInlet.GetHandle());
			TFirstGuess_Num[j + 1] = streams_data_structure[inlets_[j]].mass_flow_rate_gas*streams_data_structure[inlets_[j]].temperature*Cpstream;
			TFirstGuess_Denum[j+1] = streams_data_structure[inlets_[j]].mass_flow_rate_gas*Cpstream;
		}

		H = HInlet.SumElements()/mass_flow_rate_gas_*mix_molecular_weight; // Inlet enthalpy to the mixer
		T_guess= TFirstGuess_Num.SumElements()/ TFirstGuess_Denum.SumElements(); // First guess for the equilibrium T

		// Energy balance
		T_final = thermodynamicsMap.GetTemperatureFromEnthalpyAndMoleFractions(H, pressure_, mole_fraction.GetHandle(), T_guess);

		// Save on local outlet stream
		StreamOut.mass_flow_rate_gas = mass_flow_rate_gas_;
		StreamOut.omega_gas = omega_final;
		StreamOut.temperature = T_final;
		StreamOut.pressure = pressure_;
		StreamOut.phase = "Gas";
		StreamOut.name = streams_data_structure[outlets_[0]].name;

		// Import
		streams_data_structure[outlets_[0]] = StreamOut;

		return 0;
	};

		/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
	int Units_Mixer_Adiabatic::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure)
	{
		OpenSMOKE::FatalErrorMessage(" Residuals have been requested for a trivial unit (MIXER). Please check your input file or contact Matteo Mensi at matteo.mensi@mail.polimi.it");
		return 0;
	};

	void Units_Mixer_Adiabatic::PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		OpenSMOKE::FatalErrorMessage(" No status printer for mixer yet");
	}

}
