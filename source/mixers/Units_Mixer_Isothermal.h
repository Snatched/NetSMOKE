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

#ifndef NODUSSMOKE_UNITS_MIXER_ISOTHERMAL_H
#define NODUSSMOKE_UNITS_MIXER_ISOTHERMAL_H

// A class to evaluate isothermal mixers in a RNM

#include "Units.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace NodusSMOKE
{

class Units_Mixer_Isothermal : public Units
{
public: // FUNCTIONS

	// CONSTRUCTOR
	Units_Mixer_Isothermal(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, NodusSMOKE::UnitInfo unit_data);

	// ALTERNATIVE CONSTRUCTOR
	Units_Mixer_Isothermal(OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN *thermodynamicsMapXML, NodusSMOKE::UnitInfo unit_data);

	// MIXING OPERATION
	int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

	int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

	void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

private: // INTERNAL VARIABLES

	// For species balance
	OpenSMOKE::OpenSMOKEVectorDouble GasInletComponentMassFlow;
	OpenSMOKE::OpenSMOKEVectorDouble SolidInletComponentMassFlow;
	OpenSMOKE::OpenSMOKEVectorDouble omega_final_gas;
	OpenSMOKE::OpenSMOKEVectorDouble omega_final_solid;

	// Indexes
	int n_gas_species;
	int n_solid_species;
	int n_inlets;

	// Inlet structures
	OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMap = nullptr;
	OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN *thermodynamicsMapS = nullptr;
	NodusSMOKE::StreamInfo StreamOut;

	// Booleans
	bool contains_gas;
	bool contains_solid;

	// Temperature
	double temperature_;

};

}

#include "Units_Mixer_Isothermal.hpp"

#endif /* NODUSSMOKE_UNITS_MIXER_ISOTHERMAL_H */
