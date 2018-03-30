/*-----------------------------------------------------------------------*\
|																		  |
|			 _   _      _    _____ __  __  ____  _  ________         	  |
|			| \ | |    | |  / ____|  \/  |/ __ \| |/ /  ____|        	  |
|			|  \| | ___| |_| (___ | \  / | |  | | ' /| |__   			  |
|			| . ` |/ _ \ __|\___ \| |\/| | |  | |  < |  __|  		  	  |
|			| |\  |  __/ |_ ____) | |  | | |__| | . \| |____ 		 	  |
|			|_| \_|\___|\__|_____/|_|  |_|\____/|_|\_\______|		 	  |
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef NETSMOKE_RNM_UNITS_H
#define NETSMOKE_RNM_UNITS_H

namespace NetSMOKE
{
	/*
	The purpose of this class is to provide a common framework to deal with
	different kind of reactors and devices for a reactor network simulation.
	It is called inside the ReactorNetwork class.
	*/

	class Units
	{
	public:
		/* Default contsructor */
		Units(NetSMOKE::UnitInfo unit_data)
		{
			name_ = unit_data.name;
			type_ = unit_data.type;
			tag_ = unit_data.tag;
			energy_ = unit_data.energy;
			inlets_.resize(unit_data.inlets.size());
			inlets_ = unit_data.inlets;
			outlets_.resize(unit_data.outlets.size());
			outlets_ = unit_data.outlets;
			pressure_ = unit_data.pressure;
			phase_ = unit_data.phase;
			n_equations_ = unit_data.Nequations;
		}

		// Information accessers
		inline std::string GetName() { return name_; }							// Get device name
		inline std::string GetTag() { return tag_; }							// Get device tag - which kind of device it is
		inline std::string GetType() { return type_; }							// Get device type
		inline std::vector<int> GetOutlets() { return outlets_; }		// Get vector of corresponding outlet streams
		inline std::vector<int> GetInlets() { return inlets_; }		// Get vector of corresponding inlet streams
		inline unsigned int GetNumberOfEquations() { return n_equations_; }		// Get number of equations associated with the device for the global system
		inline std::string GetPhase() { return phase_; }						// Get phase of the reactor, no result if no reactor
		inline std::string GetEnergy() { return energy_; }						// Get energy type, isothermal or non isothermal

		/* Virtual unit solver */
		virtual int Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Virtual residual obtainer - ONLY FOR REACTORS OR NON TRIVIAL UNITS, THESE WILL GO INTO THE CONSTRUCTED NLS/ODE */
		virtual int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Virtual status printer- ONLY FOR REACTORS FOR NOW */
		virtual void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {};

		/* Virtual sequential solver */
		virtual int NonIterativeSolve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure){ return 0; };

		/* RTD aware solver */
		virtual int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals,  const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {return 0;};

	protected:

		// Identification variables
		std::string name_;
		std::string type_;
		std::string tag_;
		std::string energy_;
		std::string phase_;
		std::vector<int> inlets_;
		std::vector<int> outlets_;
		unsigned int n_equations_;

		// Valued variables
		double pressure_;
		double mass_flow_rate_gas_;
		double mass_flow_rate_solid_;

	};

}

#endif /* NETSMOKE_RNM_UNITS_H */
