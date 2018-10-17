#ifndef __RESHUFFLE_PLUGIN_HH__
#define __RESHUFFLE_PLUGIN_HH__

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
// $Id: ReshufflePlugin.hh 802 2011-10-07 18:56:01Z salam $
//-------------------------------------------------------ENDHEADER----

#include "fastjet/ClusterSequence.hh"
#include <string>
#include <cassert>

/// plugin that allows one to generate a cluster sequence that has the
/// same clusterings as some "original" sequence, but with the input
/// particles remapped according to the orig_to_new array, where
/// orig_to_new[i] indicates the new index of the particle that
/// originally had index i
class ReshufflePlugin : public fastjet::JetDefinition::Plugin {
public:
  /// constructor
  ReshufflePlugin(const fastjet::ClusterSequence * original_cs, 
		  const std::vector<int> & orig_to_new) :
    _orig_cs(original_cs), _orig_to_new(orig_to_new), _n(original_cs->n_particles()) {}

  virtual std::string description() const {
    return _orig_cs->jet_def().description() + " (reshuffled)";
  }

  virtual void  run_clustering(fastjet::ClusterSequence &) const ;

  virtual double R() const {return _orig_cs->jet_def().R();} 

private:
  /// deal with conversions from a cluster_hist_index in the original
  /// parent history to a jet index in the new history
  int jetp_new(int orig_hist_index) const {
    
    if        (orig_hist_index < 0) {
      // negative entries remain unchanged
      return orig_hist_index;
    } else if (orig_hist_index < int(_n)) {
      // entries < n are original particles and need reshuffling
      // (but the history and jet indices coincide)
      return _orig_to_new[orig_hist_index];
    } else {
      // entries >= n are recombined particles and 
      // get the jet index corresponding to this history index
      // (numbering will be the same in the old and new CS).
      return _orig_cs->history()[orig_hist_index].jetp_index;
    }
  }

  const fastjet::ClusterSequence * _orig_cs;

  std::vector<int> _orig_to_new;
  unsigned _n;
};


void ReshufflePlugin::run_clustering(fastjet::ClusterSequence & new_cs) const {

  assert(new_cs.jets().size() == _n);
  assert(new_cs.jets().size() == _n);

  const std::vector<fastjet::ClusterSequence::history_element> & orig_history = _orig_cs->history();
  for (unsigned i = _n; i < orig_history.size(); i++) {
    
    int jetp1 = jetp_new(orig_history[i].parent1);
    int jetp2 = jetp_new(orig_history[i].parent2);

//    std::cout  << "TR: "
//	       << orig_history[i].parent1 << " -> " << jetp1 << "\t\t" 
//      	       << orig_history[i].parent2 << " -> " << jetp2 << std::endl;
    int newjet_k;

    if (jetp2 >= 0) {
      new_cs.plugin_record_ij_recombination(jetp1, jetp2, orig_history[i].dij, 
					    //					    _orig_cs->jets()[orig_history[i].jetp_index],
					    newjet_k);
    } else {
      new_cs.plugin_record_iB_recombination(jetp1, orig_history[i].dij);
    }
  }
}



#endif // __RESHUFFLE_PLUGIN_HH__
