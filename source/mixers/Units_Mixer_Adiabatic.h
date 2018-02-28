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

#ifndef NODUSSMOKE_UNITS_MIXER_ADIABATIC_H
#define NODUSSMOKE_UNITS_MIXER_ADIABATIC_H

// A class to evaluate adiabatic mixers in a RNM

#include "Units.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

namespace NodusSMOKE
{

	class Units_Mixer_Adiabatic : public Units
	{
	public: // FUNCTIONS

		// CONSTRUCTOR
		Units_Mixer_Adiabatic(OpenSMOKE::ThermodynamicsMap_CHEMKIN *thermodynamicsMapXML, NodusSMOKE::UnitInfo unit_data);

		// MIXING OPERATION
		int Solve(std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		int GetResiduals(OpenSMOKE::OpenSMOKEVectorDouble &residuals, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

		void PrintStatus(boost::filesystem::path output_folder, std::vector<NodusSMOKE::StreamInfo> &streams_data_structure);

	private: // INTERNAL VARIABLES

		// For species balance
		OpenSMOKE::OpenSMOKEVectorDouble InletComponentMassFlow;
		OpenSMOKE::OpenSMOKEVectorDouble omega_final;
		OpenSMOKE::OpenSMOKEVectorDouble mole_fraction;
		double mix_molecular_weight;
		// Energy balances
		double T_guess;
		double T_final;
		double H;
		double Cpstream;
		OpenSMOKE::OpenSMOKEVectorDouble HInlet;
		OpenSMOKE::OpenSMOKEVectorDouble xInlet;
		OpenSMOKE::OpenSMOKEVectorDouble TFirstGuess_Num;
		OpenSMOKE::OpenSMOKEVectorDouble TFirstGuess_Denum;

		// Indexes
		int n_species;
		int n_inlets;

		// Inlet structures
		OpenSMOKE::ThermodynamicsMap_CHEMKIN &thermodynamicsMap;
		NodusSMOKE::StreamInfo StreamOut;

	};

}

#include "Units_Mixer_Adiabatic.hpp"

#endif /* NODUSSMOKE_UNITS_MIXER_ADIABATIC_H */
