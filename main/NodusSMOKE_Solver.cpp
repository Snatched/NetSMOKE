/*-----------------------------------------------------------------------*\
|																		  |
|	 _   _           _            _____ __  __  ____  _  ________  		  |
|	| \ | |         | |          / ____|  \/  |/ __ \| |/ /  ____| 		  |
|	|  \| | ___   __| |_   _ ___| (___ | \  / | |  | | ' /| |__    		  |
|	| . ` |/ _ \ / _` | | | / __|\___ \| |\/| | |  | |  < |  __|   		  |
|	| |\  | (_) | (_| | |_| \__ \____) | |  | | |__| | . \| |____  		  |
|	|_| \_|\___/ \__,_|\__,_|___/_____/|_|  |_|\____/|_|\_\______|		  |                                                              
|                                                                         |
|   Author: Matteo Mensi <matteo.mensi@mail.polimi.it>                    |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "OpenSMOKEpp" 				// OpenSMOKE definitions
#include "vector"  					// std vector
#include "iostream" 				// input output stream
#include "NodusSMOKE"  				// NodusSMOKE comprehensive header
#include "NodusSMOKE_Functions.h" 	// Extra functions to help with the main.cpp
#include "ReactorNetwork.h"
#include "ReactorNetwork_Gas.h"

int main(int argc, char** argv) {

	NodusSMOKE::NodusSMOKE_Logo();

	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	// Important paths
	boost::filesystem::path working_folder = boost::filesystem::current_path();
	boost::filesystem::path kinetic_folder = working_folder / "kinetics";
	boost::filesystem::path output_folder = working_folder / "output";
	boost::filesystem::path backup_folder = working_folder / "backup";
	boost::filesystem::path track_folder = working_folder / "RTD";
	OpenSMOKE::CreateDirectory(output_folder);

	// Prepare structures
	std::vector<NodusSMOKE::UnitInfo> Unit;
	std::vector<NodusSMOKE::StreamInfo> Stream;
	std::vector<int> From;
	std::vector<int> To;
	std::vector<double> relative_split;

	// Case flag
	bool solid_case;
	bool analytical_case;
	bool VideoPrint;
	bool LegacyPrint;
	double solid_density = 0.;
	bool rtd;
	bool restart;
	std::string track;

	// Prepare maps storage for kinetics
	// Gas phase mechanism
	OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML = nullptr;
	OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML = nullptr;

	// Solid phase mechanism
	OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN* thermodynamicsMapSolidXML = nullptr;
	OpenSMOKE::KineticsMap_Solid_CHEMKIN* kineticsMapSolidXML = nullptr;

	// ------------------------------------------------------------------------------------------- //
	//                                     Kinetic Mechanism		                 		       //
	// ------------------------------------------------------------------------------------------- //

	rapidxml::xml_document<> doc_gas; // Initializing the XML file
	std::vector<char> xml_string_gas;
	OpenSMOKE::OpenInputFileXML(doc_gas, xml_string_gas, kinetic_folder / "kinetics.xml");

	// Local gas map
	thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc_gas, true);
	kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc_gas, true);

	// If a solid stream or working phase has been given, interpret also the kinetic mechanism
	rapidxml::xml_document<> doc_solid;
	std::vector<char> xml_string_solid;


	if (boost::filesystem::exists( kinetic_folder / "kinetics.solid.xml")) {

		OpenSMOKE::OpenInputFileXML(doc_solid, xml_string_solid, kinetic_folder / "kinetics.solid.xml");

		// Local solid mix map
		thermodynamicsMapSolidXML = new OpenSMOKE::ThermodynamicsMap_Solid_CHEMKIN(doc_solid);
		kineticsMapSolidXML = new OpenSMOKE::KineticsMap_Solid_CHEMKIN(*thermodynamicsMapSolidXML, doc_solid, true);
	}

	// ------------------------------------------------------------------------------------------- //
	//                                     Read Input				                 		       //
	// ------------------------------------------------------------------------------------------- //

	{
		// Create the input reading system object
		NodusSMOKE::InputReader iosystem(working_folder, *thermodynamicsMapXML, *kineticsMapXML, *thermodynamicsMapSolidXML, *kineticsMapSolidXML);

		// Run the reader and generate the PDF image of the network
		iosystem.ReadAndDraw();

		// Set the case according to presence of solid species in the inlets
		iosystem.CheckSimulationType(solid_case, solid_density, analytical_case, VideoPrint, LegacyPrint, restart, rtd, track);

		// Create local maps inside the reader object to allow for data analysis
		iosystem.CreateInfoFromMaps();

		// Import corrected data structures to the main
		iosystem.ImportData(Stream, Unit, From, To, relative_split);

	}

	if (solid_case == true) {
		std::cout << std::endl;
		std::cout << "SOLID SPECIES DETECTED: solid/gas simulation will be performed." << std::endl;
		std::cout << std::endl;
	}
	else {
		std::cout << std::endl;
		std::cout << "GAS SPECIES ONLY: gas phase simulation will be performed." << std::endl;
		std::cout << std::endl;
	}

	if (analytical_case == true) {
		std::cout << std::endl;
		std::cout << "DECLARED RECYCLE ABSENCE: if the unit in the input are ordered as required by a by-hand solution, the simulation will be fast even for complex systems." << std::endl;
		std::cout << std::endl;
	}

	// ------------------------------------------------------------------------------------------- //
    //                                    Prepare device object				                 	   //
	// ------------------------------------------------------------------------------------------- //

	std::vector<NodusSMOKE::Units*> Device;

	NodusSMOKE::CreateUnitObjectsFromRawData(Unit, solid_density, Device, thermodynamicsMapXML, kineticsMapXML, thermodynamicsMapSolidXML, kineticsMapSolidXML);

	/* At this point the vector of units as objects is
	ready to use. These objects are the ones
	fed to the ReactorNetwork classes, in the solver (or main) for the sake of consistency the UnitInfo data structure
	is used. Instead of changing that as well, it would be just better to have an input reader that creates and initialize
	objects on the go with no middle-man. */


	// ------------------------------------------------------------------------------------------- //
	//                                     Indexing					                 		       //
	// ------------------------------------------------------------------------------------------- //

	const unsigned int Nunits = Unit.size();
	const unsigned int Nstreams = Stream.size();
	unsigned int Nreactors;
	unsigned int Nreactors_HE;
	unsigned int Nreactors_Solid;
	unsigned int neq;
	unsigned int NS = 0;
	unsigned int NSgas = 0;
	unsigned int NSsolid = 0;

	Nreactors = 0;
	Nreactors_HE = 0;
	Nreactors_Solid = 0;
	for (int i = 0; i < Nunits; i++) {
		if (Unit[i].tag == "Reactor") {
			++Nreactors;
			if (Unit[i].energy == "NonIsothermal") {
				++Nreactors_HE;
			}
			if (Unit[i].phase != "Gas") {
				++Nreactors_Solid;
			}
		}
	}

	// Number of species
	if (solid_case == false) {
		NS = thermodynamicsMapXML->NumberOfSpecies();
		NSgas = NS;
		NSsolid = 0;
	}
	else {
		NS = thermodynamicsMapSolidXML->NumberOfSpecies();
		NSgas = thermodynamicsMapSolidXML->number_of_gas_species();
		NSsolid = thermodynamicsMapSolidXML->number_of_solid_species();
	}

	// Number of equations
	if (solid_case == true) {
		neq = Nreactors*(NSgas + 2) + Nreactors_Solid*(NSsolid + 1);
	}
	else {
		neq = Nreactors*(NS + 1);
	}


	//CREATE TRIPLETS STORAGE STRUCTURE - MOST MEMORY OCCUPATION EFFICIENT METHOD TO BUILD OUR MATRIXES

	typedef Eigen::Triplet<int> triplet_struct_int;

	std::vector<Eigen::Triplet<int>> triplets_vector_FlowMatrix;

	for (unsigned int i = 0; i < Nstreams; ++i) { // rows
		triplets_vector_FlowMatrix.push_back(triplet_struct_int(To[i], i, +1));
		triplets_vector_FlowMatrix.push_back(triplet_struct_int(From[i], i, -1));
	}

	// ------------------------------------------------------------------------------------------- //
	//                              First guess generation                                         //
	// ------------------------------------------------------------------------------------------- //

	// First guess declaration
	OpenSMOKE::OpenSMOKEVectorDouble first_guess(neq);

	if (analytical_case == false || restart == true) {
	std::cout << "Building first guess solution..." << std::endl;
		// FIRST GUESS GENERATION
		/*----------------------------------------*/
		/*			if gas only					  */
		/*----------------------------------------*/
		if (solid_case == false) {
			NodusSMOKE::GasPhaseCase_FirstGuessGenerator(Unit, Stream, From, first_guess, thermodynamicsMapXML, kineticsMapXML);
		}

		/*----------------------------------------*/
		/*			if solid					  */
		/*----------------------------------------*/
		if (solid_case == true) {
			NodusSMOKE::SolidPhaseCase_FirstGuessGenerator(solid_density, Unit, Stream, From, first_guess, thermodynamicsMapXML, kineticsMapXML, thermodynamicsMapSolidXML, kineticsMapSolidXML);
		}

		std::cout << "DONE" << std::endl;

	}

	// File print preparation
	std::ofstream ASCIIStreamSummary;



	// Initialize Network

	if (solid_case == false) {

		// ------------------------------------------------------------------------------------------- //
		//                                    Gas Phase Network		                                   //
		// ------------------------------------------------------------------------------------------- //

		// Construct Unit/Unit triplet sparse matrix, useful for one time solving of a LS with direct method
		// This allows to know all mass flow rates of the system

		typedef Eigen::Triplet<double> triplet_struct_double;

		std::vector<Eigen::Triplet<double>> triplets_vector_UnitMatrix;

		if (Nunits <= 1) {
			triplets_vector_UnitMatrix.push_back(triplet_struct_double(0, 0, -1.));
		}
		else if (Nunits > 1) {
			for (unsigned int i = 0; i < relative_split.size(); ++i) { // rows
				triplets_vector_UnitMatrix.push_back(triplet_struct_double(To[i], From[i], relative_split[i]));
			}
			for (unsigned int i = 0; i < Nunits + 1; ++i) { // diagonal
				triplets_vector_UnitMatrix.push_back(triplet_struct_double(i, i, -1.));
			}
		}

	
		NodusSMOKE::ReactorNetwork_Gas myRNM(thermodynamicsMapXML, kineticsMapXML, &Device, &Stream, first_guess, &triplets_vector_FlowMatrix, &triplets_vector_UnitMatrix, analytical_case);
		std::cout << "Everything is set-up. The solver is being run..." << std::endl;
		// Solve Network

		if (restart == true) {
			boost::filesystem::path backup_file = backup_folder / "Backup.bin";
			myRNM.InitializeFromBackUp(backup_file);
		}

		double tStartSolver = OpenSMOKE::OpenSMOKEGetCpuTime();
		if (analytical_case == true) {
			myRNM.SolveAsSeries();
		}
		else if (rtd == true) {
			myRNM.SetTrack(track);
			myRNM.PrepareTrackHistory(track_folder);
			myRNM.SparsityPattern();
			myRNM.SolveODE(1.e8);
			myRNM.CloseTrackFile();
		}
		else {
			myRNM.SparsityPattern();
			myRNM.SequentialSubstitution(100);
			myRNM.SolveODE(1.e8);
			myRNM.SolveNLS();
		}
		double tEndSolver = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time required: " << tEndSolver - tStartSolver << std::endl;

		myRNM.WriteBackup();

		// Printing
		myRNM.SetMWsAndMoleFractions();

		if (VideoPrint == true) {
			myRNM.PrintOutletStreams();
		}

		myRNM.PrepareStreamSummary(ASCIIStreamSummary, output_folder / "StreamSummary.out");
		ASCIIStreamSummary.setf(std::ios::scientific);
		myRNM.PrintStreamSummary(ASCIIStreamSummary);
		ASCIIStreamSummary.close();

		if (LegacyPrint == true) {
			myRNM.PrintFinalReactorsStatus(output_folder);
		}
	}


	else if (solid_case == true) {

		// ------------------------------------------------------------------------------------------- //
		//                              Multiphase Gas & Solid Network                                 //
		// ------------------------------------------------------------------------------------------- //

		NodusSMOKE::ReactorNetwork_Solid myRNM(thermodynamicsMapXML, kineticsMapXML, thermodynamicsMapSolidXML, kineticsMapSolidXML, solid_density, &Device, &Stream, first_guess, &triplets_vector_FlowMatrix, analytical_case);
		std::cout << "Everything is set-up. The solver is being run..." << std::endl;
		// Solve Network
		if (restart == true) {
			myRNM.InitializeFromBackUp(backup_folder);
		}

		double tStartSolver = OpenSMOKE::OpenSMOKEGetCpuTime();
		if (analytical_case == true) {
			myRNM.SolveAsSeries();
		}
		else if (rtd == true) {
			myRNM.SetTrack(track);
			myRNM.PrepareTrackHistory(track_folder);
			myRNM.SparsityPattern();
			myRNM.SolveODE(1.e8);
			myRNM.CloseTrackFile();
		}
		else {
			myRNM.SparsityPattern();
			myRNM.SequentialSubstitution(100);
			myRNM.SolveODE(1.e8);
			myRNM.SolveNLS();
		}
		double tEndSolver = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time required: " << tEndSolver - tStartSolver << std::endl;

		myRNM.WriteBackup();

		// Force correct stream phases to display
		for (int i = 0; i<Nstreams; ++i) {
			myRNM.SetStreamPhase(i);
		}

		// Printing
		myRNM.SetMWsAndMoleFractions();

		if (VideoPrint == true) {
			myRNM.PrintOutletStreams();
		}

		myRNM.PrepareStreamSummary(ASCIIStreamSummary, output_folder / "StreamSummary.out");
		ASCIIStreamSummary.setf(std::ios::scientific);
		myRNM.PrintStreamSummary(ASCIIStreamSummary);
		ASCIIStreamSummary.close();

		if (LegacyPrint == true) {
			myRNM.PrintFinalReactorsStatus(output_folder);
		}


	}


	return OPENSMOKE_SUCCESSFULL_EXIT;

};
