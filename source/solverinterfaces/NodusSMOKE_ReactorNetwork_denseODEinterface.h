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

#ifndef NODUSSMOKE_REACTORNETWORK_DENSEODEINTERFACE_H
#define NODUSSMOKE_REACTORNETWORK_DENSEODEINTERFACE_H

#include "ReactorNetwork.h"

// ODE Solver
#include "math/native-ode-solvers/MultiValueSolver"
#include "NodusSMOKE_ReactorNetwork_denseODEinterface.h"

class NodusSMOKE_ReactorNetwork_denseODEinterface
{
public:

	NodusSMOKE_ReactorNetwork_denseODEinterface() {};

	void SetReactor(NodusSMOKE::ReactorNetwork* reactor)
	{
		reactor_ = reactor;
	};

protected:

	unsigned int ne_;

	void MemoryAllocation()
	{
		OpenSMOKE::ChangeDimensions(ne_, &y_, true);
		OpenSMOKE::ChangeDimensions(ne_, &dy_, false);
	};

	virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
	{
		y_.CopyFrom(Y.data());
		reactor_->ODEEquations(t, y_, dy_); // Insert variables in the reactor
		dy_.CopyTo(DY.data());
	};

	virtual void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) {};

	void Print(const double t, const Eigen::VectorXd &Y)
	{	
		reactor_->PrintODE(t, dy_);
	};

private:

	NodusSMOKE::ReactorNetwork* reactor_;
	OpenSMOKE::OpenSMOKEVectorDouble  y_;
	OpenSMOKE::OpenSMOKEVectorDouble dy_;
};

#endif /* NODUSSMOKE_REACTORNETWORK_DENSEODEINTERFACE_H */