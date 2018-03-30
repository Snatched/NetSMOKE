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

#ifndef NETSMOKE_UNITS_REACTORS_PSR_H
#define NETSMOKE_UNITS_REACTORS_PSR_H

// Parent class
#include "Units.h"

// Relative reactor
#include "idealreactors/psr/PerfectlyStirredReactor"

namespace NetSMOKE
{

	class Units_Reactors_PSR : public Units {
	public:
		/* Default constructor */
		Units_Reactors_PSR(NetSMOKE::UnitInfo unit_data) : 
			Units(unit_data)
				{
					residence_time_ = unit_data.residence_time;
				};

		/* Solve */
		virtual int Solve(std::vector<NetSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Get Residuals */
		virtual int GetResiduals( OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) { return 0; };

		/* Print reactor status */
		virtual void PrintStatus(boost::filesystem::path output_folder, std::vector<NetSMOKE::StreamInfo> &streams_data_structure){};

		/* RTD aware solver */
		virtual int RTD(OpenSMOKE::OpenSMOKEVectorDouble &residuals, const double t, std::vector<NetSMOKE::StreamInfo> &streams_data_structure) {return 0;};
		
	protected:
	
		// for oscillation integrator
		bool is_oscillating;

		// Working structures
		NetSMOKE::StreamInfo StreamOut;
		NetSMOKE::StreamInfo StreamIn;
		NetSMOKE::StreamInfo StreamOut_star;	// given by solver

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

#endif /* NETSMOKE_UNITS_REACTORS_PSR_H */
