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

#ifndef NETSMOKE_REACTORNETWORK_SOLID_H
#define NETSMOKE_REACTORNETWORK_SOLID_H

// Include the Reactor Network Main Class
#include "ReactorNetwork.h"
#include "NetSMOKE_Units"
#include "NetSMOKE_Utilities"


	//!  A class for simulating multi phase (gas/solid) reactor networks
	/*!
	The purpose of this class is to simulate reactor network which have multiphase streams.
	*/

namespace NetSMOKE
{

	class ReactorNetwork_Solid : public ReactorNetwork 

	{


		public:

			/**
			*@brief Default constructor
			*@param thermodynamicsMap Thermodynamic interpreted XML file
			*@param kineticsMap Kinetics interpreted XML file
			*@param thermodynamicsMapSOLID Solid thermodynamics for phase recognition
			*@param kineticsMapSOLID Solid kinetics for solid phase reactions
			*@param rho_solid solid density
			*@param DeviceMap vector of network objects
			*@param Stream structure containing the streams informations stored in a vector
			*@param FirstGuessY vector containing the first guess solution
			*@param InputTripletsVector Matrix stream/units
			*@input_flag_sequential_solver Boolean to skip some unnecessary step when the sequential serial solver is used
			*/
			ReactorNetwork_Solid
			(
				OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMap,
				OpenSMOKE::KineticsMap_CHEMKIN* kineticsMap,
				OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMapSolid,
				OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMapSolid,
				double rho_solid,
				std::vector<NetSMOKE::Units*> *DeviceMap,
				std::vector<NetSMOKE::StreamInfo>* StreamMap,
				OpenSMOKE::OpenSMOKEVectorDouble FirstGuessY,
				std::vector<Eigen::Triplet<int>>* InputTripletsVector,
				bool input_flag_sequential_solver
			);

				//  Set values - used to change the values of the reactor structure by the solver
				/**
				*@brief  Assign the current solution vector to the network
				*@param itervalues current solution vector provided by the solver
				*/
				void SetValues(OpenSMOKE::OpenSMOKEVectorDouble itervalues);


				// Network equations in NLS - equations for all the species in all the Unit composing the NLS
				/**
				*@brief  Equation accessible to an NLS solver
				*@param y current solution vector
				*@param dy current residuals vector
				*/
				int NLSEquations(OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy);

				// Network equations in ODE- equations for all the species in all the Unit composing the ODE system
				/**
				*@brief  Equations accessible to an ODE solver
				*@param t current time
				*@param y current solution vector
				*@param dy current residuals vector
				*/
				int ODEEquations(const double t, OpenSMOKE::OpenSMOKEVectorDouble &y, OpenSMOKE::OpenSMOKEVectorDouble &dy);

				// Print - not sure but needed to fit for the solver error request
				/**
				*@brief  Print function to display steps performed by the solver
				*@brief y current solution vector
				*/
				int Print(const OpenSMOKE::OpenSMOKEVectorDouble& y);

				// Network Solving NLS - Run the solver of the NLS
				/**
				*@brief  Constructs and runs the NLS solver
				*/
				void SolveNLS();

				// Network Solving ODE - Run the solver of the ODE system
				/**
				*@brief  Construct and runs the ODE solver
				*@param tf time of integration
				*/
				void SolveODE(const double tf);
				
				// Set stream phase
				/**
				*@brief Checks the composition of the stream to set the phase between solid, gas and mix accordingly
				*/
				void SetStreamPhase(int i);

				// Set MWs and rho - assign to each stream a proper molecular weight, density of mix and mole fraction
				/**
				*@brief Assign to each stream a proper molecular weight, density of mix and mole fraction
				*/
				void SetMWsAndMoleFractions();

				// Print output file - generates a summary output containing all the streams and their composition
				/**
				*@brief Calls the fileprinting of the stream status
				*/
				void PrepareStreamSummary(std::ofstream& fOutput, const boost::filesystem::path output_file_ascii);
				void PrintStreamSummary(std::ostream& fOutput);

				// Test stuff rtd
				void SetTrack(std::string track_this);
				void PrintTrack(double t);
				void PrepareTrackHistory(boost::filesystem::path output_folder_track);


		private:

			// Solid/gas heterogeneous reactions maps
			OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN thermodynamicsMapXMLinternalSOLID;
			OpenSMOKE::KineticsMap_Solid_CHEMKIN kineticsMapXMLinternalSOLID;

			// Solid particle density
			double rho_solid;


	};


}


#include "ReactorNetwork_Solid.hpp"

#endif /* NETSMOKE_REACTORNETWORK_SOLID_H */
