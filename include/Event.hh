#ifndef __EVENT_HH__
#define __EVENT_HH__

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
// $Id: Event.hh 2251 2014-10-28 15:00:15Z sapeta $
//-------------------------------------------------------ENDHEADER----

#include <vector>
#include "Flavour.hh"
#include "FlavourPlugin.hh"
#include <iostream>
#include <limits>

#define fj fastjet

//----------------------------------------------------------------------
/// class that helps dealing with flavour and decay products
class LSParticle : public fj::PseudoJet {
public:
  LSParticle(){};
  LSParticle(double px, double py, double pz, double E, int flavour,
	     std::vector<fj::PseudoJet> decay_products = std::vector<fj::PseudoJet>());
  LSParticle(const fj::PseudoJet & p, int flavour,
	     std::vector<fj::PseudoJet> decay_products = std::vector<fj::PseudoJet>());
  
  const Flavour & flavour() const {return _f;}
  const std::vector<fj::PseudoJet> & decay_products() const {return _decay_products;}
  void set_flavour(Flavour f) {_f=f;}
  void set_decay_products(std::vector<fj::PseudoJet> dp) {_decay_products=dp;}
  
  // only changes the (px,py,pz,E) components without modifying the user_index
  // and the cluster_hist_index (we also check that the flavour is not changed).
  void reset_4vector(const fj::PseudoJet & p);

  /// recoil the decay products of this particle either given the vector k
  /// for the boost (in which case k_is_known=true), or given the old momentum
  /// for this particle (i.e. before reshuffling,
  /// in which case k_is_known=false). [Intended for internal use]
  void recoil_decay_products();

protected:
  Flavour _f;
  std::vector<fj::PseudoJet> _decay_products;
};

/// class for holding one event, i.e a set of particles and a
/// corresponding weight.
class Event {
public:
  std::vector<LSParticle> particles;
  double weight;
  std::vector<double> weights;

  std::vector<LSParticle> Zbosons() const;
  std::vector<LSParticle> Wbosons() const;
  std::vector<LSParticle> partons() const;
  std::vector<LSParticle> Hbosons() const;
  //std::vector<Flavour> flavours;
  //the momenta wrt which we make the boost once some beam particles 
  //become virtual
  //std::vector<fastjet::PseudoJet> k; 

  /// return the sum of all the particles' momenta
  fj::PseudoJet sum() const;
  /// return the scalar sum of the transverse momenta
  double HT() const;
};

/// print events in the form:
//  px,  py,  pz, E (if cart = true)
//  pt, rap, phi, m (if cart = false)
void print_event(const Event & ev, const bool cart = true, 
                 const bool with_decay_products = true,
		 std::ostream* ostr = (&std::cout));

void print_jet(const fj::PseudoJet & p);

void print_jet_ppvars(const fj::PseudoJet & p);


#endif //__EVENT_HH__
