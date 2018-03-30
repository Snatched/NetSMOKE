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

#ifndef NETSMOKE_UNITS_REACTORS_PFR_H
#define NETSMOKE_UNITS_REACTORS_PFR_H

// Parent class
#include "Units.h"

// Relative reactor
#include "idealreactors/plugflow/PlugFlowReactor"

namespace NetSMOKE
{

	class Units_Reactors_PFR : public Units {
	public:
		/* Default constructor */
		Units_Reactors_PFR(NetSMOKE::UnitInfo unit_data) : 
			Units(unit_data)
				{
					residence_time_ = unit_data.residence_time;
				};

		/* Solve */
		virtual int Solve() { return 0; };

		/* Get Residuals */
		virtual int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Print reactor status */
		virtual void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {};

		/* RTD aware solver */
		virtual int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {return 0;};


	protected:
		// Generic
		double residence_time_;
		double volume;

		// Working structs
		NetSMOKE::StreamInfo StreamIn;
		NetSMOKE::StreamInfo StreamOut;
		NetSMOKE::StreamInfo StreamOut_star; // dictated by solver

		// Indexes
		int Nstreams;
		int unit_index;  // Used to identify the index of the unit

		// Static things
		// Options
		OpenSMOKE::PlugFlowReactor_Options* plugflow_options;

		// ODE Parameters
		OpenSMOKE::ODE_Parameters* ode_parameters;

		// On the fly ROPA
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA;

		// On the fly PostProcessing
		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing;

		// Polimi soot
		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot;
	};

}

#endif /* NETSMOKE_UNITS_REACTORS_PFR_H */
