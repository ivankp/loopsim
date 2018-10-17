//-------------------------------------------------------BEGINHEADER--
// This file is part of LoopSim (arXiv:1006.2144).
//
// Copyright 2009-2011 Mathieu Rubin, Gavin P. Salam and Sebastian Sapeta
//
// This software is provided as a snapshot of the ongoing LoopSim
// research project. As such, you may not redistribute it without the
// authors' permission.
//
// If you use LoopSim as part of your scientific work, you should
// discuss and agree with the LoopSim authors how best to acknowledge
// LoopSim in your work (e.g. whether through a reference, or through
// joint authorship with the LoopSim authors). 
//
// To help guide LoopSim's proper use in its current early stage of
// development, a condition of use of LoopSim is that any results that
// you obtain with it must be shown to the LoopSim authors before
// they are made public.
//
// $Id: FlavourPlugin.cc 802 2011-10-07 18:56:01Z salam $
//-------------------------------------------------------ENDHEADER----

#include "FlavourPlugin.hh"
#include "fastjet/NNH.hh"
#include "Flavour.hh"
#include <vector>
#include <limits>

using namespace std;

class FlavourBriefJet {
public:
  void init(const fastjet::PseudoJet & jet, FlavourPlugin * fp) {
    _jet = jet;
    _flav.set_flavour_from_user_index(jet.user_index());
    _R2 = pow(fp->R(),2);
  }

  double distance(const FlavourBriefJet * jet) const {
    //if i and j are flavour incompatible, then return infinity
    if (!(this->_flav).can_be_recombined(jet->_flav)) {
      return numeric_limits<double>::max();
    }
    double dij = _jet.squared_distance(jet->_jet)/_R2;
    //If one of them is a boson, and the other one a parton
    if ((this->_flav).is_a_vector_boson()!=(jet->_flav).is_a_vector_boson()) {
      double diff; // = E_ti^2*DeltaR_ij^2/R^2-pt_j^2
      if ((this->_flav).is_a_vector_boson()) {
	diff = (this->_jet).mperp2()*dij - (jet->_jet).perp2();
      }
      else {
	diff = (jet->_jet).mperp2()*dij - (this->_jet).perp2();
      }
      if (diff>0) dij = numeric_limits<double>::max();
    }
    return dij;
  }

  double beam_distance() const {
    return 1.0;
  }

private:
  fastjet::PseudoJet _jet;
  Flavour _flav;
  double _R2;
};


string FlavourPlugin::description() const {
  ostringstream desc;
  desc << "Flavour algorithm for LoopSim";
  return desc.str();
}


void FlavourPlugin::run_clustering(fastjet::ClusterSequence & cs) 
  const {
  int njets = cs.jets().size();
  //define a non-const FlavourPlugin object in order to be used with
  //NNH constructor
  FlavourPlugin copy_this = *this;
  fastjet::NNH<FlavourBriefJet,FlavourPlugin> nnh(cs.jets(),&copy_this);

  while (njets>0) {
    int i,j,k;
    double dij = nnh.dij_min(i,j);
    if (j>=0) {
      cs.plugin_record_ij_recombination(i,j,dij,k);
      nnh.merge_jets(i,j,cs.jets()[k],k);
    }else {
      cs.plugin_record_iB_recombination(i,dij);
      nnh.remove_jet(i);
    }
    njets--;
  }

}
