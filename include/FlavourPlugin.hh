#ifndef __FLAVOUR_PLUGIN_HH__
#define __FLAVOUR_PLUGIN_HH__

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
// $Id: FlavourPlugin.hh 802 2011-10-07 18:56:01Z salam $
//-------------------------------------------------------ENDHEADER----

#include <string>
#include <fastjet/ClusterSequence.hh>

namespace loopsim {

class FlavourPlugin : public fastjet::JetDefinition::Plugin {
public:
  ///R = radius of the algorithm
  FlavourPlugin(double R){_R = R;}
  virtual std::string description() const;
  virtual void run_clustering(fastjet::ClusterSequence &) const;
  virtual double R() const {return _R;}
private:
  double _R;
};

}

#endif // __FLAVOUR_PLUGIN_HH__
