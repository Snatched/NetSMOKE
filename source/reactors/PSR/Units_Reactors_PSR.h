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

#ifndef NODUSSMOKE_UNITS_REACTORS_PSR_H
#define NODUSSMOKE_UNITS_REACTORS_PSR_H

// Parent class
#include "Units.h"

// Relative reactor
#include "idealreactors/psr/PerfectlyStirredReactor"

namespace NodusSMOKE
{

	class Units_Reactors_PSR : public Units {
	public:
		/* Default constructor */
		Units_Reactors_PSR(NodusSMOKE::UnitInfo unit_data) : 
			Units(unit_data)
				{
					residence_time_ = unit_data.residence_time;
				};

		/* Solve */
		virtual int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Get Residuals */
		virtual int GetResiduals( OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Print reactor status */
		virtual void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure){};

		/* RTD aware solver */
		virtual int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure) {return 0;};
		
	protected:
	
		// for oscillation integrator
		bool is_oscillating;

		// Working structures
		NodusSMOKE::StreamInfo StreamOut;
		NodusSMOKE::StreamInfo StreamIn;
		NodusSMOKE::StreamInfo StreamOut_star;	// given by solver

		// Indexes
		int Nstreams;
		int unit_index;  // Used to identify the index of the unit
		int NCgas;
		int NCsolid;

		// defined variables
		double residence_time_;
		double temperature_;

		// densities
		double rhoGas;

		// storage for production rates
		OpenSMOKE::OpenSMOKEVectorDouble Rgas;

		// temporary omega storage
		OpenSMOKE::OpenSMOKEVectorDouble temp_omega_gas_out;

		// storage for closer initial conditions
		OpenSMOKE::OpenSMOKEVectorDouble Omega_previous_gas;
		
		int n_calls;
	};
	
}

#endif /* NODUSSMOKE_UNITS_REACTORS_PSR_H */
