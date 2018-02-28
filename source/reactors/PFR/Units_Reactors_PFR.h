/*-----------------------------------------------------------------------*\
|																		  |
|	 _   _           _            _____ __  __  ____  _  ________  		  |
|	| \ | |         | |          / ____|  \/  |/ __ \| |/ /  ____| 		  |
|	|  \| | ___   __| |_   _ ___| (___ | \  / | |  | | ' /| |__    		  |
|	| . ` |/ _ \ / _` | | | / __|\___ \| |\/| | |  | |  < |  __|   		  |
|	| |\  | (_) | (_| | |_| \__ \____) | |  | | |__| | . \| |____  		  |
|	|_| \_|\___/ \__,_|\__,_|___/_____/|_|  |_|\____/|_|\_\______|		  |                                                              |
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef NODUSSMOKE_UNITS_REACTORS_PFR_H
#define NODUSSMOKE_UNITS_REACTORS_PFR_H

// Parent class
#include "Units.h"

// Relative reactor
#include "idealreactors/plugflow/PlugFlowReactor"

namespace NodusSMOKE
{

	class Units_Reactors_PFR : public Units {
	public:
		/* Default constructor */
		Units_Reactors_PFR(NodusSMOKE::UnitInfo unit_data) : 
			Units(unit_data)
				{
					residence_time_ = unit_data.residence_time;
				};

		/* Solve */
		virtual int Solve() { return 0; };

		/* Get Residuals */
		virtual int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Print reactor status */
		virtual void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {};

		/* RTD aware solver */
		virtual int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {return 0;};


	protected:
		// Generic
		double residence_time_;
		double volume;

		// Working structs
		NodusSMOKE::StreamInfo StreamIn;
		NodusSMOKE::StreamInfo StreamOut;
		NodusSMOKE::StreamInfo StreamOut_star; // dictated by solver

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

#endif /* NODUSSMOKE_UNITS_REACTORS_PFR_H */
