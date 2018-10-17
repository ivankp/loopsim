#ifndef __TREELEVEL_HH__
#define __TREELEVEL_HH__

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
// $Id: TreeLevel.hh 802 2011-10-07 18:56:01Z salam $
//-------------------------------------------------------ENDHEADER----

#include <vector>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SharedPtr.hh"
#include "Flavour.hh"
#include "FlavourPlugin.hh"
#include <cassert>
#include <iostream>
#include <limits>
#include "Event.hh"

#define fj fastjet


enum STATUS {
  undef=-4, //< undefined (but will be soon defined!)
  born=-3,  //< born particle
  nc=-2,    //< non-clustering (particle that will not become virtual but which is not a born particle)
  beam=-1   /// particle that clusters only with the beam (so becomes virtual)
};

class TreeLevel {

public:
  
  /// nborn = number of born particles that we want in the event
  /// R = radius of the Flavour Algorithm
  TreeLevel(const Event & ev, int nborn, double R = 1, bool spread_virtuals = false);

  /// Returns the maximum number of loops that has been generated (GPS interpretation? MR: that's true)
  int number_of_virtuals() {return _nvirtuals;}
  const std::vector<std::vector<Event> > & extract_all_events();
  /// extract all events with iloops number of loops
  const std::vector<Event> & extract_all_events(int iloops);
  /// return the event that was given in the constructor, ie _all_events[0][0]
  const Event & extract_original_event();
  /// return the event _all_events[i][j];
  const Event & extract_event(int i, int j);

  /// Returns number of events at i loops (total number if iloops is omitted).
  ///
  /// Returns -1 if iloops >= 0 and the requested number of loops does
  /// not exist (GPS: why not 0? MR: I don't know, I have strange ideas sometimes!)
  int number_of_events(int iloops = -1);


  void print_all_events();
  void print_status();

private:

  /// function that sets the _status for all particles
  /// R is the radius of the Flavour Algorithm
  void init(double R);
  // 
  /// function that updates _status information for particles indicated
  /// by the user indices of p1 and p2 (latter may be NULL), based on what
  /// happened in the p1,p2 clustering.
  ///
  /// on input ip is the number of particles whose status was already defined
  /// and on exit may have been updated
  ///
  /// jet_index is expected to be equal to (p1+p2).user_index() on input
  ///
  /// p1 and p2 are (pointers to) the pseudojets that come from one
  /// part of the clustering history.
  ///
  /// [NB: is called by init()]
  // 
  void fill_status(int & ip, int jet_index, const fj::PseudoJet * p1,
		   const fj::PseudoJet * p2 = NULL);


  /// Function that finds all particles that can become virtual
  /// and put their position in the initial event into the vector
  /// _pseudo_virtuals while giving them a user_index equal to their status.
  /// To make analysis easier, change ordering of particles, placing
  /// pseudo virtual particles at the end of the vector _no_loop_event.particles
  ///
  /// [NB: is called by init()]
  void finish_init();



  /// function that works through the list of particles
  /// that can be made virtual (_pseudo_virtual), and for each
  /// particle two options are followed (recursing in each case):
  /// 
  ///  - it is left as is
  ///  - it is made virtual
  ///
  /// The recursion's end point is reached when ip = _nvirtual and
  /// at that stage the recombinations are actually made (carefully, cf
  /// function make_virtual(...))
  /// 
  /// \param iloop is the number of loops decided on so far
  ///        (GPS+SS: presumably redundant)
  /// \param ip is the index of the particle being considered
  /// \param virtuals is the list of indices of the particles
  ///        to be made virtual
  ///
  void generate_all_relevant_loops(int iloops, int ip, 
				   std::vector<int> virtuals);


  /// Take all particles that are virtuals in the event, and try to make
  /// them virtual in a "physical" way:
  /// 1- If a virtual particle is looped over a final state parton p, then
  /// recombine it with p and put the resulting 4-vector mass to 0
  /// 2- Take all virtual particles which are looped over a parton from the
  /// beam, and make a boost of all real (recombined) particles in order
  /// to get a null total pt vector.
  ///
  /// [NB: called by generate_all_relevant_loops]
  Event make_virtual(const std::vector<int> & virtuals);

  /// alternative version that spreads "looped" momenta over
  /// all real children of their partner.
  Event make_virtual_spread(const std::vector<int> & virtuals);

  
  /// given an event evt, sort out its recoils, including transverse
  /// momentum recoil and setting particles back to the correct mass.
  ///
  /// @param evt   :  the event prior to sorting out recoils
  /// @param vbeam :  the amount of momentum that has been recombined with the beam
  void do_recoils(Event & evt, const fj::PseudoJet & vbeam);

  /// check that there's not missing pt, and if there is, output 
  /// a warning
  void verify_no_missing_pt(const Event & ev, const char * label = 0, bool print_orig = false) const;

  Event _no_loop_event;
  //_all_events[i] contains all events at i loops.
  std::vector<std::vector<Event> > _all_events;
  int _nparticles, _nborn;
  bool _spread_virtuals;

  /// contains index of particles whose status is either "beam" or >=0
  /// i.e. particles that can become virtual
  /// This index corresponds to the place in the _no_loop_event.particles
  /// vector.
  std::vector<int> _pseudo_virtuals;
  int _nvirtuals; /// = _pseudo_virtuals.size()
  /// number of particles that will stay real in all versions
  /// of this event
  int _nreals; /// = _nparticles - _nvirtuals

  /// Status for all particles in the no_loop_event.
  ///
  /// status[i] is one of {undef,born,nc,beam} (all -ve) or is >= 0.
  /// If status is >= 0, this indicates that it clusters with the
  /// particle whose index=status.
  //
  // Not sure it is really useful to keep it as a member of the class.
  // We'll see later...
  std::vector<int> _status;

  /// contains the "original" cluster sequence
  fastjet::SharedPtr<fastjet::ClusterSequence> _cs;
  /// contains the "Reshuffling" jet definition
  fastjet::SharedPtr<fastjet::JetDefinition::Plugin>   _reshuffle_plugin;
  fastjet::JetDefinition           _new_jet_def;
  /// contains the "Reshuffled" cluster sequence
  fastjet::SharedPtr<fastjet::ClusterSequence> _new_cs;

  static constexpr double _zero = 1e-8;
  static constexpr double _missing_pt_tolerance = 1e-4;
  // a small value is useful for diagnosing issues in funny events
  //const static double _missing_pt_tolerance = 1e-13;
};


/// class to help us establish auxiliary info in the clustering
/// Note that it relies on the specific form of encoding the CS history in FJ.
class AuxHist {
public:
  AuxHist(const fj::ClusterSequence & cs, int hindex) : hist_index(hindex) {
    const fj::ClusterSequence::history_element & hist = cs.history()[hindex];
    if (hist.parent1 < 0) {
      //-- the history entry for an original particle
      //aux_distance = 0;
      // slightly modified version of Mathieu's suggestion to
      // force minimal aux_distances when parents are ill-defined.
      aux_distance = -std::numeric_limits<double>::max();
    } else if (hist.parent2 < 0) {
      //-- the history entry for a particle-beam merging
      // GPSToCheck: to generalise beyond C/A, this should probably be
      // aux_distance = cs.jets()[hist.jetp_index].mperp2();
      aux_distance = cs.jets()[hist.parent1].mperp2();
    } else {
      //-- the history entry for a particle-particle merging
      // GPSToCheck: to generalise beyond C/A, this should probably be
      // int jet1_index = cs.history[hist.parent1].jetp_index;
      // int jet2_index = cs.history[hist.parent1].jetp_index;
      // aux_distance = cs.jets()[jet1_index].kt_distance(cs.jets()[jet2_index]);
      aux_distance = cs.jets()[hist.parent1].kt_distance(cs.jets()[hist.parent2]);
    }
  }

  int     hist_index;
  double  aux_distance;
};
// comparison operator
bool operator<(const AuxHist & aux1, const AuxHist & aux2);


typedef fastjet::JetDefinition::DefaultRecombiner DefRecomb;
// forward def

//----------------------------------------------------------------------
/// class that does same things as the JetDefinition::DefaultRecombiner
/// plus propagates the most energetic (in the Et sense) particle's index
/// while keeping track of flavour information
/// Beware that it assumes that the user_index's 2 most significant bytes
/// encode flavour. (See Flavour::translate(int) for more information).
/// MEAF: MostEnergeticAndFlavour
class MEAFRecombiner : public  DefRecomb {
public:
  MEAFRecombiner(fastjet::RecombinationScheme recomb_scheme = 
                    fastjet::E_scheme) : 
    DefRecomb(recomb_scheme) {};

  virtual std::string description() const {return DefRecomb::description()
      +" (with propagation of the most energetic particle's index)";}

  /// recombine pa and pb and put result into pab
  virtual void recombine(const fastjet::PseudoJet & pa, 
                         const fastjet::PseudoJet & pb, 
                         fastjet::PseudoJet & pab) const {
    DefRecomb::recombine(pa,pb,pab);
    pab.set_user_index(recombined_user_index(pa,pb));

  }

  /// returns the user_index associated with the recombination of pa and pb.
  ///
  /// The user index has its flavour part set according to the
  /// Flavour::recombine function. Its jet index part is equal to that
  /// of the parent that is "closer" to the child, as determined by a
  /// combination of flavour considerations (flavour::is_closer) and
  /// kinematics.
  int recombined_user_index(const fastjet::PseudoJet & pa,
			    const fastjet::PseudoJet & pb) const;

};



fj::PseudoJet PxPyPzM(const double px, const double py, 
		      const double pz, const double M = 0.0);


#endif //__TREELEVEL_HH__
