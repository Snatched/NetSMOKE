
namespace NetSMOKE
{

	// CONSTRUCTOR
	Units_Splitter::Units_Splitter(UnitInfo unit_data) :
	Units(unit_data)
	
	{

		n_outlets_ = unit_data.outlets.size();
		split_factor_.resize(n_outlets_);
		split_factor_ = unit_data.split_factor;
		n_equations_ = 0;

	};

	// SPLITTING OPERATION
	int Units_Splitter::Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {

		for (int i = 1; i <= n_outlets_; ++i) {
			//TONOTICE Can this be done in a more robust way?
			streams_data_structure[outlets_[i-1]].phase = streams_data_structure[inlets_[0]].phase;

			streams_data_structure[outlets_[i-1]].mass_flow_rate_gas = streams_data_structure[inlets_[0]].mass_flow_rate_gas*split_factor_[i-1];
			streams_data_structure[outlets_[i-1]].mass_flow_rate_solid = streams_data_structure[inlets_[0]].mass_flow_rate_solid*split_factor_[i-1];

			streams_data_structure[outlets_[i-1]].omega_gas = streams_data_structure[inlets_[0]].omega_gas;
			streams_data_structure[outlets_[i-1]].omega_solid = streams_data_structure[inlets_[0]].omega_solid;

			streams_data_structure[outlets_[i-1]].temperature = streams_data_structure[inlets_[0]].temperature;
			
		}
		return 0;
	};

	/* Residual obtainer - RETURNS ERROR IF YOU CALL RESIDUALS FOR TRIVIAL UNITS */
	int Units_Splitter::GetResiduals( OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		OpenSMOKE::FatalErrorMessage(" Residuals have been requested for a trivial unit (SPLITTER). Please check your input file or contact Matteo Mensi at matteo.mensi@mail.polimi.it");
		return 0;
	};

	/* Status printer - RETURNS ERROR */
	void Units_Splitter::PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){
		OpenSMOKE::FatalErrorMessage("No print function implemented for trivial units yet");
	};

}
