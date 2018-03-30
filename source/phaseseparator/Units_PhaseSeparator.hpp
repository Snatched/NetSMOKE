namespace NetSMOKE

{
	// CONSTRUCTOR
	Units_PhaseSeparator::Units_PhaseSeparator(OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermo,NetSMOKE::UnitInfo unit_data) :
		thermodynamicsMapSOLID(*thermo),
		Units(unit_data)
	{

		// identify the streams_data_structure

		n_outlets_ = 2;
		splitted_phase_.resize(n_outlets_);
		splitted_phase_=unit_data.splitted_phase;
		
		// Allocate memory
		StreamOut.resize(n_outlets_);
		for (int i =0; i<n_outlets_; ++i){
			OpenSMOKE::ChangeDimensions(thermo->number_of_gas_species(),&StreamOut[i].omega_gas, true);
			OpenSMOKE::ChangeDimensions(thermo->number_of_solid_species(),&StreamOut[i].omega_solid, true);
		}

	};

	// PHASE SPLITTING OPERATION
	int Units_PhaseSeparator::Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure)
	{
		// set local working variables
		StreamIn = streams_data_structure[inlets_[0]];
		StreamOut[0] = streams_data_structure[outlets_[0]];
		StreamOut[1] = streams_data_structure[outlets_[1]];
		StreamOut[0].name = streams_data_structure[outlets_[0]].name;
		StreamOut[1].name = streams_data_structure[outlets_[1]].name;
		StreamOut[0].phase = splitted_phase_[0];
		StreamOut[1].phase = splitted_phase_[1];
		
		// Perform splitting operation
		for (int i = 0; i < n_outlets_; ++i) {
			StreamOut[i].temperature = StreamIn.temperature;
			//if solid do
			if (StreamOut[i].phase == "Solid"){
				StreamOut[i].omega_solid = StreamIn.omega_solid;
				StreamOut[i].mass_flow_rate_solid = StreamIn.mass_flow_rate_solid;

				StreamOut[i].mass_flow_rate_gas = 0.;

			}
			// if gas do
			if (StreamOut[i].phase == "Gas"){
				StreamOut[i].omega_gas = StreamIn.omega_gas;
				StreamOut[i].mass_flow_rate_gas = StreamIn.mass_flow_rate_gas;

				StreamOut[i].mass_flow_rate_solid = 0.;
			}
		}

		// Return the info to the main
		streams_data_structure[outlets_[0]]=StreamOut[0];
		streams_data_structure[outlets_[1]]=StreamOut[1];

		return 0;
	};

	/* Residual obtainer - Arguments are not used but needed for homogeneity */
	int Units_PhaseSeparator::GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
			
		OpenSMOKE::FatalErrorMessage(" Residuals have been requested for a trivial unit (PHASESPLITTER). Please check your input file or contact Matteo Mensi at matteo.mensi@mail.polimi.it");

		return 0;
	};

	/* Printer - RETURNS ERROR */
	void Units_PhaseSeparator::PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		OpenSMOKE::FatalErrorMessage("Can't print phase splitters yet brah");
	}

}
