/*-----------------------------------------------------------------------*\
|																		  |
|	_   _           _            _____ __  __  ____  _  ________  		  |
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



#include "vector"
#include "string"
#include "OpenSMOKEpp"

#ifndef NODUSSMOKE_DATASTRUCTURES_H
#define NODUSSMOKE_DATASTRUCTURES_H

/* The purpose of this header is to provide the common structure for data storage and transmission */

namespace NodusSMOKE
{
	// Structure for units
	struct UnitInfo {
		std::string name;
		std::string tag;
		double pressure;
		int flow_matrix_index;

		// According to me COHERENT WITH MATRIXES
		std::vector<int> outlets; // Needed to assign correctly unit flowrates to streams, need to use it till better way found
		std::vector<int> inlets; // Inlet streams (presented as corresponding indexes in the Stream vector of StreamInfo variables
		
		// According to silvio NOT COHERENT WITH ALL OTHER DATA STRUCTURES
		std::vector<int> inlets_names;
		std::vector<int> outlets_names;
		std::vector<std::string> splitted_phase;

		// For reactors
		std::string type;
		std::string energy; // also used in mixers now PogChamp
		std::string phase;
		double residence_time;
		int Nequations;

		// Used only in massflows distribution linear system for gas phase only case
		double mass_flow_rate_solid;
		double mass_flow_rate_gas;

			// Isothermal
			double temperature;

			// NonIsothermal
			double OutletT;
			double UA;
			double T_environment;

			// For iterating on T
			double InletT;

		// For mixers - No particular input information besides energy and temperature

		// For splitters
		std::vector<double> split_factor;


	};

	// Structure for streams

	struct StreamInfo {

		int name;
		std::string phase;

		double temperature;
		double rho_gas;
		double MW_gas;
		double rho_solid;
		double MW_solid;
		double pressure;

		double mass_flow_rate_solid;
		double mass_flow_rate_gas;
		OpenSMOKE::OpenSMOKEVectorDouble omega_gas;
		OpenSMOKE::OpenSMOKEVectorDouble x_gas;
		OpenSMOKE::OpenSMOKEVectorDouble omega_solid;
		OpenSMOKE::OpenSMOKEVectorDouble x_solid;

	};
}
#endif /* NODUSSMOKE_DATASTRUCTURES_H */
