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

/*
Collection of auxiliary functions used in the NetSMOKE solver
These functions are not needed for correct functioning of NetSMOKE as a framework
but are useful in the author's use case.
*/

#ifndef NETSMOKE_FUNCTIONS_H
#define NETSMOKE_FUNCTIONS_H

// For first guess generators
#include "idealreactors/psr/PerfectlyStirredReactor"
#include "solidreactors/solidpsr/SolidPerfectlyStirredReactor"

namespace NetSMOKE
{

    void NetSMOKE_Logo(){

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
        
		std::cout << R"(
		 _   _      _    _____ __  __  ____  _  ________ 
		| \ | |    | |  / ____|  \/  |/ __ \| |/ /  ____|
		|  \| | ___| |_| (___ | \  / | |  | | ' /| |__   
		| . ` |/ _ \ __|\___ \| |\/| | |  | |  < |  __|  
		| |\  |  __/ |_ ____) | |  | | |__| | . \| |____ 
		|_| \_|\___|\__|_____/|_|  |_|\____/|_|\_\______|   
		)"
		<< '\n';
                                                                                                                                          
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "           	 Matteo Mensi <matteo.mensi@mail.polimi.it>        		   " << std::endl;
		std::cout << "           Department of Chemistry, Materials and Chemical Engineering       " << std::endl;
		std::cout << "                            Politecnico di Milano                            " << std::endl;
		std::cout << "                        http://www.opensmoke.polimi.it/                      " << std::endl;
		std::cout << "                      http://creckmodeling.chem.polimi.it/                   " << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                                  WARNING                                    " << std::endl;
		std::cout << "    This version of NetSMOKE can be used for educational purposes only.    " << std::endl;
		std::cout << "       The software is and remains the sole property of Matteo Mensi.        " << std::endl;
		std::cout << "           Whenever NetSMOKE is used to produce any publication,     	   " << std::endl;
		std::cout << "           a detailed reference to the OpenSMOKE++ and NetSMOKE            " << std::endl;
		std::cout << "                  code should be reported (see User's Guide).                " << std::endl;
		std::cout << "    Use for commercial purposes is not permitted. For any commercial issue   " << std::endl;
		std::cout << "         please contact Alberto Cuoci (email: alberto.cuoci@polimi.it)       " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                            LIMITED WARRANTY                                 " << std::endl;
		std::cout << "     This software is provided \"as is\" and without warranties as to        " << std::endl;
		std::cout << "  performance of merchantability or any other warranties whether expressed   " << std::endl;
		std::cout << "    or implied. Because of the various hardware and software environments    " << std::endl;
		std::cout << "   into which this library may be installed, no warranty of fitness for a    " << std::endl;
		std::cout << "   particular purpose is offered. The user must assume the entire risk of    " << std::endl;
		std::cout << "                          using  the library.                                " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
    }

	int GasPhaseCase_FirstGuessGenerator(std::vector<NetSMOKE::UnitInfo> Unit, std::vector<NetSMOKE::StreamInfo> Stream, std::vector<int> From,
	OpenSMOKE::OpenSMOKEVectorDouble &first_guess,
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
	OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML)	
	{

		int NS = thermodynamicsMapXML->NumberOfSpecies();
		int Nstreams = Stream.size();
		int Nunits = Unit.size();
		int Nreactors = 0;
		for (int i = 0; i < Unit.size(); ++i) {
			if (Unit[i].tag == "Reactor") {
				Nreactors++;
			}
		}

		OpenSMOKE::OpenSMOKEVectorDouble omega_first(NS);
		double first_guess_T;

		std::vector<double> residence_time_vect;
		for (unsigned int i = 0; i < Nunits; ++i) {
			if (Unit[i].tag == "Reactor")
				residence_time_vect.push_back(Unit[i].residence_time);
		}

		// I use the minimum residence time reduced
		int target = std::distance(residence_time_vect.begin(), std::min_element(residence_time_vect.begin(), residence_time_vect.end()));
		double tau_1st = 0.5*residence_time_vect[target];

		std::vector<NetSMOKE::StreamInfo> stream_inlet_vect;
		NetSMOKE::StreamInfo StreamIn;
		StreamIn.pressure = Stream[0].pressure;

		OpenSMOKE::ChangeDimensions(NS, &StreamIn.omega_gas, true);
		for (unsigned int i = 0; i < Nstreams; i++) {
			if (From[i] == 0) {
				stream_inlet_vect.push_back(Stream[i]);
			}
		}

		StreamIn.temperature = stream_inlet_vect[0].temperature;

		if (stream_inlet_vect.size() == 1) {
			StreamIn = stream_inlet_vect[0];
		}
		else {
			OpenSMOKE::OpenSMOKEVectorDouble StreamInMassFlow(stream_inlet_vect.size());
			for (unsigned int i = 1; i <= stream_inlet_vect.size(); ++i) {
				StreamInMassFlow[i] = stream_inlet_vect[i - 1].mass_flow_rate_gas;
			}
			StreamIn.mass_flow_rate_gas = StreamInMassFlow.SumElements();

			OpenSMOKE::OpenSMOKEVectorDouble InletComponentMassFlow(NS);
			for (int i = 1; i <= NS; ++i) {
				for (int j = 0; j < stream_inlet_vect.size(); ++j) {
					InletComponentMassFlow[j + 1] = stream_inlet_vect[j].omega_gas[i] * stream_inlet_vect[j].mass_flow_rate_gas;
				}
				StreamIn.omega_gas[i] = InletComponentMassFlow.SumElements() / StreamIn.mass_flow_rate_gas;
			}
			StreamIn.temperature = stream_inlet_vect[0].temperature;
		}

		OpenSMOKE::PerfectlyStirredReactor_Options* psr_options = new OpenSMOKE::PerfectlyStirredReactor_Options();
		psr_options->SetVerboseVideo(false);
		psr_options->SetVerboseASCIIFile(false);
		psr_options->SetVerboseOutput(false);
		psr_options->SetVerboseXMLFile(false);
		OpenSMOKE::ODE_Parameters* ode_parameters = new OpenSMOKE::ODE_Parameters();
		OpenSMOKE::OnTheFlyROPA* onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML, *kineticsMapXML);

		OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, psr_options->output_path());

		OpenSMOKE::PolimiSoot_Analyzer* polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);

		OpenSMOKE::PerfectlyStirredReactor_Isothermal_ConstantPressure psr(*thermodynamicsMapXML, *kineticsMapXML,
			*ode_parameters, *psr_options, *onTheFlyROPA, *on_the_fly_post_processing, *polimi_soot,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas,
			tau_1st, -1., StreamIn.mass_flow_rate_gas);
		// Solve PSR
		psr.Solve(1.e8);
		double P;
		psr.GetFinalStatus(first_guess_T, P, omega_first);

		std::vector<double> first_guess_dynamic_vector;
		for (unsigned int i = 0; i < Nreactors; ++i)
		{
			for (unsigned int j = 1; j <= NS; ++j)
				first_guess_dynamic_vector.push_back(omega_first[j]);
			if (Unit[i].energy == "NonIsothermal")
				first_guess_dynamic_vector.push_back(first_guess_T);
			else if (Unit[i].energy == "Isothermal")
				first_guess_dynamic_vector.push_back(Unit[i].temperature);
		}
		first_guess.CopyFrom(first_guess_dynamic_vector.data());

		return 0;
	};

	int SolidPhaseCase_FirstGuessGenerator(double solid_density, std::vector<NetSMOKE::UnitInfo> Unit, std::vector<NetSMOKE::StreamInfo> Stream, std::vector<int> From,
	OpenSMOKE::OpenSMOKEVectorDouble &first_guess,
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
	OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML,
	OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMapSolidXML,
	OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMapSolidXML)
	{
		int NSsolid = thermodynamicsMapSolidXML->number_of_solid_species();
		int NSgas = thermodynamicsMapSolidXML->number_of_gas_species();
		int NS = NSsolid + NSgas;
		int Nstreams = Stream.size();
		int Nunits = Unit.size();
		int Nreactors = 0;
		for (int i = 0; i < Unit.size(); ++i) {
			if (Unit[i].tag == "Reactor") {
				Nreactors++;
			}
		}


		OpenSMOKE::OpenSMOKEVectorDouble omega_first_solid(NSsolid);
		OpenSMOKE::OpenSMOKEVectorDouble omega_first_gas(NSgas);
		double first_guess_flow;
		double first_guess_flow_gas;
		double first_guess_flow_solid;
		double first_guess_T;



		std::vector<double> residence_time_vect;
		for (unsigned int i = 0; i < Nunits; ++i) {
			if (Unit[i].tag == "Reactor")
				residence_time_vect.push_back(Unit[i].residence_time);
		}

		// I use the minimum residence time reduced
		int target = std::distance(residence_time_vect.begin(), std::min_element(residence_time_vect.begin(), residence_time_vect.end()));
		double tau_1st = 0.5*residence_time_vect[target];

		std::vector<StreamInfo> stream_inlet_vect;
		StreamInfo StreamIn;
		StreamIn.pressure = Stream[0].pressure;
		OpenSMOKE::ChangeDimensions(NSgas, &StreamIn.omega_gas, true);
		OpenSMOKE::ChangeDimensions(NSsolid, &StreamIn.omega_solid, true);
		for (unsigned int i = 0; i < Nstreams; i++) {
			if (From[i] == 0) {
				stream_inlet_vect.push_back(Stream[i]);
			}
		}

		StreamIn.temperature = stream_inlet_vect[0].temperature;

		if (stream_inlet_vect.size() == 1) {
			StreamIn = stream_inlet_vect[0];
		}
		else {
			// Create a set of mass flows
			OpenSMOKE::OpenSMOKEVectorDouble GasStreamInMassFlow(stream_inlet_vect.size());
			OpenSMOKE::OpenSMOKEVectorDouble SolidStreamInMassFlow(stream_inlet_vect.size());
			for (unsigned int i = 1; i <= stream_inlet_vect.size(); ++i) {
				GasStreamInMassFlow[i] = stream_inlet_vect[i - 1].mass_flow_rate_gas;
				SolidStreamInMassFlow[i] = stream_inlet_vect[i - 1].mass_flow_rate_solid;
			}
			StreamIn.mass_flow_rate_gas = GasStreamInMassFlow.SumElements();
			StreamIn.mass_flow_rate_solid = SolidStreamInMassFlow.SumElements();
			// Create a set of compositions
			OpenSMOKE::OpenSMOKEVectorDouble GasInletComponentMassFlow(NSgas);
			OpenSMOKE::OpenSMOKEVectorDouble SolidInletComponentMassFlow(NSsolid);
			for (int i = 1; i <= NSgas; ++i) {
				for (int j = 0; j < stream_inlet_vect.size(); ++j) {
					GasInletComponentMassFlow[j + 1] = stream_inlet_vect[j].omega_gas[i] * stream_inlet_vect[j].mass_flow_rate_gas;
				}
				StreamIn.omega_gas[i] = GasInletComponentMassFlow.SumElements() / StreamIn.mass_flow_rate_gas;
			}
			for (int i = 1; i <= NSsolid; ++i) {
				for (int j = 0; j < stream_inlet_vect.size(); ++j) {
					SolidInletComponentMassFlow[j + 1] = stream_inlet_vect[j].omega_solid[i] * stream_inlet_vect[j].mass_flow_rate_solid;
				}
				StreamIn.omega_solid[i] = SolidInletComponentMassFlow.SumElements() / StreamIn.mass_flow_rate_solid;
			}
			StreamIn.temperature = stream_inlet_vect[0].temperature;
		}

		OpenSMOKE::SolidPerfectlyStirredReactor_Options *spsr_options = new  OpenSMOKE::SolidPerfectlyStirredReactor_Options();
		OpenSMOKE::ODE_Parameters *ode_parameters = new OpenSMOKE::ODE_Parameters();
		spsr_options->SetVerboseOutput(false);
		spsr_options->SetVerboseVideo(false);
		spsr_options->SetVerboseASCIIFile(false);
		spsr_options->SetVerboseXMLFile(false);

		OpenSMOKE::SolidPerfectlyStirredReactor_Isothermal_ConstantPressure psr(*thermodynamicsMapXML, *kineticsMapXML,
			*thermodynamicsMapSolidXML, *kineticsMapSolidXML, *ode_parameters, *spsr_options,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas, StreamIn.omega_solid,
			StreamIn.temperature, StreamIn.pressure, StreamIn.omega_gas, StreamIn.omega_solid,
			StreamIn.mass_flow_rate_gas, StreamIn.mass_flow_rate_solid, solid_density,
			tau_1st, -1.);

		first_guess_flow_gas = StreamIn.mass_flow_rate_gas;
		first_guess_flow_solid = StreamIn.mass_flow_rate_solid;
		// Solve PSR
		psr.Solve(1.e8);
		double P, dummyMgas, dummyMsolid;
		psr.GetFinalStatus(first_guess_T, P, dummyMgas, dummyMsolid, omega_first_gas, omega_first_solid);

		std::vector<std::string> ReactorPhases;
		std::vector<double> first_guess_dynamic_vector;

		for (unsigned int i = 0; i < Unit.size(); ++i) {
			if (Unit[i].tag == "Reactor") {
				ReactorPhases.push_back(Unit[i].phase);
			}
		}

		for (unsigned int i = 0; i < Nreactors; ++i)
		{
			// Assign composition
			if (ReactorPhases[i] == "Mix") {
				for (unsigned int j = 1; j <= NS; ++j) {
					if (j <= NSgas)
						first_guess_dynamic_vector.push_back(omega_first_gas[j]);
					else if (j>NSgas)
						first_guess_dynamic_vector.push_back(omega_first_solid[j - NSgas]);
				}
				// Assign mass flow
				first_guess_dynamic_vector.push_back(first_guess_flow_gas);
				first_guess_dynamic_vector.push_back(first_guess_flow_solid);

			}
			else if (ReactorPhases[i] == "Gas") {
				for (unsigned int j = 1; j <= NSgas; ++j) {
					if (j <= NSgas)
						first_guess_dynamic_vector.push_back(omega_first_gas[j]);
				}
				// Assign mass flow
				first_guess_dynamic_vector.push_back(first_guess_flow_gas);
			}
			else if (ReactorPhases[i] == "Solid") {
				for (unsigned int j = 1; j <= NS; ++j) {
					if (j <= NSgas)
						first_guess_dynamic_vector.push_back(0.);
					else if (j>NSgas)
						first_guess_dynamic_vector.push_back(omega_first_solid[j - NSgas]);
				}
				// Assign mass flow
				first_guess_dynamic_vector.push_back(0.);
				first_guess_dynamic_vector.push_back(first_guess_flow_solid);
			}

			// Assign T
			if (Unit[i].energy == "NonIsothermal")
				first_guess_dynamic_vector.push_back(first_guess_T);
			else if (Unit[i].energy == "Isothermal")
				first_guess_dynamic_vector.push_back(Unit[i].temperature);

		}

		first_guess.CopyFrom(first_guess_dynamic_vector.data());
		return 0;
	}

	int CreateUnitObjectsFromRawData(std::vector<NetSMOKE::UnitInfo> Unit,
									double rho_solid,
									std::vector<NetSMOKE::Units*> &Device,
									OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
									OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML,
									OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMapSolidXML,
									OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMapSolidXML)
	{	
		// Allocate Memory
		Device.resize(Unit.size());
		int Nunits = Unit.size();
		
		// Construct reactor and units objects
		for (int j = 0; j < Nunits; ++j) {
			if (Unit[j].tag == "Reactor") {
				if (Unit[j].type == "PSR") {
					if (Unit[j].energy == "Isothermal") {
						if (Unit[j].phase == "Gas"){
							Device[j] = new NetSMOKE::Units_Reactors_PSR_Isothermal(thermodynamicsMapXML, kineticsMapXML, Unit[j] );
						}
						else if (Unit[j].phase == "Solid" || Unit[j].phase == "Mix"){
							Device[j] = new NetSMOKE::Units_Reactors_SolidPSR_Isothermal(thermodynamicsMapXML, kineticsMapXML, thermodynamicsMapSolidXML, kineticsMapSolidXML, rho_solid, Unit[j]);
						}
					}
					else if (Unit[j].energy == "NonIsothermal") {
						Device[j] = new NetSMOKE::Units_Reactors_PSR_NonIsothermal(thermodynamicsMapXML, kineticsMapXML, Unit[j]);

					}
				}
				if (Unit[j].type == "PFR") {
					if (Unit[j].energy == "Isothermal") {
						Device[j] = new NetSMOKE::Units_Reactors_PFR_Isothermal(thermodynamicsMapXML, kineticsMapXML, Unit[j]);

					}
					else if (Unit[j].energy == "NonIsothermal") {
						Device[j] = new NetSMOKE::Units_Reactors_PFR_NonIsothermal(thermodynamicsMapXML, kineticsMapXML, Unit[j]);

					}
				}
			}
			else if (Unit[j].tag == "Mixer"){
				if (Unit[j].energy == "Adiabatic"){
					Device[j] = new NetSMOKE::Units_Mixer_Adiabatic(thermodynamicsMapXML, Unit[j]);

				}
				else if (Unit[j].energy == "Isothermal"){
					if (rho_solid > 0.) {
						Device[j] = new NetSMOKE::Units_Mixer_Isothermal(thermodynamicsMapSolidXML, Unit[j]);
					}
					else {
						Device[j] = new NetSMOKE::Units_Mixer_Isothermal(thermodynamicsMapXML, Unit[j]);
					}

				}
			}
			else if (Unit[j].tag == "Splitter"){
				Device[j] = new NetSMOKE::Units_Splitter(Unit[j]);
			}
			else if (Unit[j].tag == "PhaseSplitter"){
				Device[j] = new NetSMOKE::Units_PhaseSeparator(thermodynamicsMapSolidXML, Unit[j]);
			}
		}
		return 0;
	}
}

#endif /* NETSMOKE_FUNCTIONS_H  */