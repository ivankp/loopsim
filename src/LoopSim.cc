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
// $Id: LoopSim.cc 2005 2013-08-23 15:12:43Z sapeta $
//-------------------------------------------------------ENDHEADER----

/** \mainpage The LoopSim library documentation

    <b> TODO before beta3 release: </b>

      \li rename create-tarball-6.6-0.0.5.sh according to the oficial MCFM 
          version
      \li archive create-tarball-6.6.sh and related files/entries in doc
      \li remove MCFM 6.6 from this wiki
    
    <hr> 

    This is the documentation for the LoopSim library, implementing
    the method described in 
    <a href="http://arxiv.org/abs/arXiv:1006.2144">arXiv:1006.2144</a>
    to generate unitarity-based approximations to loop contributions.

    \section Installation

    - tar xzvf loopsim-1.0beta3.tar.gz 
    - cd loopsim-1.0beta3
 
    - In the makefile, adjust INSTALLDIR (and FASTJETCONFIG, in case you do not
      have it in your PATH)
      - optional: MCFM 6.6-0.0.5 interface
        - adjust MCFMHOME and LHADIR in the main LoopSim makefile (compiling
          with LHAPDF allows for more flexibility)
        - in your MCFM installation directory:
           - ./Install
           - set LHAPDFLIB in the MCFM makefile and create the link 
             to PDFsets directory in Bin
           - make mcfmlib
      - optional: MCFM 6.6 interface [obsolete]
        - adjust MCFMHOME and LHADIR in the main LoopSim makefile (compiling
          with LHAPDF allows for more flexibility)
        - install mcfm 6.6 mods using mcfm-6.6-mods/install-mods.py script:
          $ ./mcfm-6.6-mods/install-mods.py path_to_your_MCFM_installation
        - in your MCFM installation directory:
           - ./Install
           - set LHAPDFLIB in the MCFM makefile and create the link 
             to PDFsets directory in Bin
           - make mcfmlib

    - make
    - make install   


    \section Example

    The example directory includes a simple analysis of \f$W^-+{\rm jet}\f$
    inclusive production. Two of them are stand-alone and use pre-generated 
    event files. Others make use of MCFM, which generates the events on 
    the fly. The main program is wjets.cc, the main analysis class is 
    LSAnalysisLHC. The send-jobs scripts allow to produce results at LO, NLO,
    LO and \f$\rm \bar nNLO\f$. The input directory contains card files for MCFM
    and LoopSim. The latter has the following parameters:

       \param [outfilename] name of the output file with all histograms 
       \param [nborn]       number of Born particles
       \param [RLS]         jet radius of the C/A alg. used to attribute 
                            emission  sequence to an event
       \param [sv]          spread momentum of looped particles in a collinear
                            safe way (default should be 1) 
       \param [makels]      make LoopSim on the given set of events from 
                            LHE file
       \param [verbose]     print some extra information

    <em>
    The distributions obtained from these examples serve for illustrative
    purposes only!  
    
    For efficiency reasons we set the ptjet_min parameter of MCFM to
    10 GeV, which may affect the lowest bins of the distributions as discussed
    in <a href="http://arxiv.org/abs/arXiv:1307.2252">arXiv:1307.2252</a>.
    Similarly, we set the ncall1, ncall2, itmx1, itmx2 parameters to moderate
    values so the distributions may show non-negligible fluctuations, especially
    at nNLO. The fastest converging results, in terms of statistics, are those
    for \f$H_{\rm T,tot}\f$ and  \f$H_{\rm T,jet}\f$. 
    To get paper-quality distributions you need to increase the ncall and itmx
    parameters and perform a set of MCFM runs with different seeds. It is also
    recommended to produce optimized MCFM grids first.

    (The approximate run times correspond to 2 core 2.8GHz CPU with 8GB of RAM.)
    </em>

    \subsection standalone Stand-alone runs: 

    The events are stored in Les Houches Event files, which can be downloaded
    from <a href="http://www.hepforge.org/downloads/loopsim">www.hepforge.org/downloads/loopsim</a>.

    For \f$W^-+{\rm jet}\f$ inclusive at nLO execute:

       - $ ./wjets input/input_LS_Wm1j_sf1_lordseed0.dat < Wm_1jet_lord_sf1_example_seed0.lhe

       - $ ./wjets input/input_LS_Wm2j_sf1_lordseed0.dat < Wm_2jet_lord_sf1_example_seed0.lhe

    The nLO result is the sum of the histograms in the two above output files.

    \subsection mcfminterfaced MCFM-interfaced runs: 
    
    The events are generated on the fly by MCFM, then piped to LoopSim, and
    finally passed to the analysis class.

    For \f$W^-+{\rm jet}\f$ inclusive at nLO execute:

       - send-jobs-Wm1jlordsf1seed0.sh (~1min)
       - send-jobs-Wm2jlordsf1seed0.sh (~2min)

    The nLO result is the sum of the histograms in the two above output files.

    For \f$W^-+{\rm jet}\f$ inclusive at nNLO execute:

       - send-jobs-Wm1jlordsf1seed0.sh (~1min)  

       - send-jobs-Wm1jrealsf1seed0.sh (~5min) 
       - send-jobs-Wm1jvirtsf1seed0.sh (~1min)

       - send-jobs-Wm2jrealglusf1seed0.sh (~1.5h) 
       - send-jobs-Wm2jrealqrksf1seed0.sh (~1.5h)
       - send-jobs-Wm2jvirtglusf1seed0.sh (~20min)
       - send-jobs-Wm2jvirtqrksf1seed0.sh (~20min)

    The nNLO result is the sum of the histograms in the seven above output
    files.


 **/

#include "LoopSim.hh"
#include <iostream>
#include <cstdlib>

using namespace std;


LoopSim::LoopSim(int order, int iloops, const Event & ev, 
		 double R, int nborn, bool spread_virtuals) {
  if (iloops>order) {
    cout << "ERROR in constructor LoopSim::LoopSim(*): ";
    cout << "the number of loops in the current event cannot be greater "; 
    cout << "than the order in perturbation theory you are ";
    cout << "currently evaluating!" << endl;
    exit(-1);
  }

  _nborn = nborn;
  _spread_virtuals = spread_virtuals;
  _order = order;
  _nloops_gen = _order-iloops;
  _R = R;

  _tree_levels.resize(_nloops_gen+1);
  _tree_levels[0].push_back(TreeLevel(ev,nborn,R,_spread_virtuals));
  for (int iloop_start = 1; iloop_start<=_nloops_gen; iloop_start++) {
    // fill _tree_levels[iloop_start]
    generate_all_tree_levels(iloop_start);
  }
  init_indexes();
}


void LoopSim::generate_all_tree_levels(int ilstart) {
  for (int j=0; j<ilstart; j++) {
    int iloops = ilstart-j;
    for (unsigned int k=0; k<_tree_levels[j].size(); k++) {
      // we have to check that there are enough virtual particles in order
      // to extract events with iloops number of loops
      if (_tree_levels[j][k].number_of_virtuals()>=iloops) {
	const vector<Event> & evt = 
	  _tree_levels[j][k].extract_all_events(iloops);
	generate_some_tree_levels(evt,ilstart);
      }
    }
  }
}

void LoopSim::generate_some_tree_levels(const vector<Event> & evt,int ilstart) {
  for (unsigned int i=0; i<evt.size(); i++) {
    Event this_evt;
    // copy evt[i] in this_evt, changing weight sign
    for (unsigned j=0; j<evt[i].particles.size(); j++) {
      this_evt.particles.push_back(evt[i].particles[j]);
      //this_evt.flavours.push_back(evt[i].flavours[j]);
    }
    this_evt.weight = -evt[i].weight;
    //this_evt.k = evt[i].k;
    _tree_levels[ilstart].push_back(TreeLevel(this_evt,_nborn,_R,_spread_virtuals));
  }
}

void LoopSim::init_indexes() {
  _tl_i = -1;
  _tl_j = 0;
  _ev_i = 0;
  _ev_j = 0;
  _tl_i_max = _nloops_gen;
  _tl_j_max = 0;
  _ev_i_max = 0;
  _ev_j_max = 0;
  _done = false;
}

const Event & LoopSim::extract_next_event() {
  if (!_done) {
    cout << "ERROR: Function LoopSim::there_is_a_next_event() ";
    cout << "must be called!" << endl;
    exit(-1);
  }
  _done = false;
  return _tree_levels[_tl_i][_tl_j].extract_event(_ev_i,_ev_j);
}

bool LoopSim::there_is_a_next_event() {
  _done = true;
  _ev_j++;
  if (_ev_j>_ev_j_max) {
    _ev_j=0;
    _ev_i++;
    if (_ev_i>_ev_i_max) {
      _ev_i=0;
      _tl_j++;
      if (_tl_j>_tl_j_max) {
	_tl_j=0;
	_tl_i++;
	if (_tl_i>_tl_i_max) {
	  return false;
	}
	_tl_j_max = _tree_levels[_tl_i].size()-1;
	if (_tl_j_max<0) return false;
      }
      _ev_i_max = _tree_levels[_tl_i][_tl_j].number_of_virtuals();
    }
    _ev_j_max = _tree_levels[_tl_i][_tl_j].number_of_events(_ev_i)-1;
  }
  return true;
}


void LoopSim::print_all_events(bool cart) {

  //save indexes
  int tl_i = _tl_i, tl_j = _tl_j, ev_i = _ev_i, ev_j = _ev_j;
  int tl_i_max = _tl_i_max, tl_j_max = _tl_j_max;
  int ev_i_max = _ev_i_max, ev_j_max = _ev_j_max;

  //in case some events were already extracted before
  init_indexes();

  int old_tl_i = _tl_i, old_tl_j = _tl_j;
  while (this->there_is_a_next_event()) {
    const Event & ev = this->extract_next_event();
    if (old_tl_i!=_tl_i || old_tl_j!=_tl_j) {
      cout << "all events extracted from _tree_levels[" << _tl_i << "]";
      cout << "[" << _tl_j << "]:" << endl << endl;
      old_tl_i = _tl_i;
      old_tl_j = _tl_j;
      _tree_levels[_tl_i][_tl_j].print_status();
    }
    print_event(ev,cart);
    cout << endl;
  }

  //put indexes to original values
  _tl_i = tl_i;
  _tl_j = tl_j;
  _ev_i = ev_i;
  _ev_j = ev_j;
  _tl_i_max = tl_i_max;
  _tl_j_max = tl_j_max;
  _ev_i_max = ev_i_max;
  _ev_j_max = ev_j_max; 

}

