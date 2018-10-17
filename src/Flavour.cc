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
// $Id: Flavour.cc 2253 2014-10-30 17:34:18Z sapeta $
//-------------------------------------------------------ENDHEADER----

#include <iostream>
#include <cassert>
#include <cstdlib>
#include "Flavour.hh"

using namespace std;

//reference: "QCD and Collider Physics", PDG 2013 for Higgs
const double Z_BOSON_MASS = 91.19;
const double W_BOSON_MASS = 80.33;
const double TOP_QUARK_MASS = 178.0;
const double HIGGS_BOSON_MASS = 125.9;

// constructors
Flavour::Flavour(int pdg_id) { _set_flavour(pdg_id); }
Flavour::Flavour(int pdg_id, double m) { _set_flavour(pdg_id,m); }

// the old ones for backward compatibility
Flavour::Flavour(char flav) {
  if (flav=='q') {
    _set_flavour(81);
  } else if (flav=='b') {
    _set_flavour(5);
  } else if (flav=='t') {
    _set_flavour(6);
  } else if (flav=='Z') {
    _set_flavour(23);
  } else if (flav=='W') {
    _set_flavour(24);
  }
}

Flavour::Flavour(char flav, double m) {
  if (flav=='q') {
    _set_flavour(81,m);
  } else if (flav=='b') {
    _set_flavour(5,m);
  } else if (flav=='t') {
    _set_flavour(6,m);
  } else if (flav=='Z') {
    _set_flavour(23,m);
  } else if (flav=='W') {
    _set_flavour(24,m);
  }
}

// do not follow Gavin's temptation to remove this routine (since we
// have another one below with masses specified), in view of future
// work e.g. that might include W+b = top (then we would need to
// specify a top mass; though we'd probably want to discard the result
// of the recombination anyway -- diagrams with missing vector bosons
// don't count...).
void Flavour::_set_flavour(int flav) {
  _pdg_id = flav;
  _set_all_flags_to_false();
  if (_pdg_id==81 || _pdg_id==5) {
    _is_a_parton = true;
    _mass = 0.0; //all quarks are supposed massless here.
  } else if (_pdg_id==6) {
    _is_a_parton = true;
    _mass = TOP_QUARK_MASS;
  } else if (_pdg_id==23) {
    _is_a_Z = true;
    _mass = Z_BOSON_MASS;
  } else if (_pdg_id==24) {
    _is_a_Wp = true;
    _mass = W_BOSON_MASS;
  } else if (_pdg_id==-24) {
    _is_a_Wm = true;
    _mass = W_BOSON_MASS;
  } else if (_pdg_id==22) {
    _is_a_gamma = true;
    _mass = 0;
  } else if (_pdg_id==25) {
    _is_a_H = true;
    _mass = HIGGS_BOSON_MASS;
  } else {
    cout << "ERROR in Flavour class (default mass): Unrecognized flavour " << flav << "!" << endl;
    exit(-1);
  }
}

void Flavour::_set_flavour(int flav, double m) {
  _pdg_id = flav;
  _mass = m;
  _set_all_flags_to_false();
  if (_pdg_id==81 || _pdg_id==5 || _pdg_id==6) {
    _is_a_parton = true;
  } else if (_pdg_id==23 ) {
    _is_a_Z = true;
  } else if (_pdg_id == 22) {
    _is_a_gamma = true;
  } else if (_pdg_id==24) {
    _is_a_Wp = true;
  } else if (_pdg_id==-24) {
    _is_a_Wm = true;
  } else if (_pdg_id==25 ) {
    _is_a_H = true;
  } else {
    cout << "ERROR in Flavour class (specified mass): Unrecognized flavour " << flav << "!" << endl;
    exit(-1);
  }  
}


//       |gamma|  q  |  W+ |  W- |  Z  |  H
// -------------------------------------------
// gamma |  -  |  +  |  +  |  +  |  -  |  -
// -------------------------------------------
//    q  |/////|  -  |  +  |  +  |  +  |  -
// -------------------------------------------
//    W+ |/////|/////|  -  |  +  |  +  |  +
// -------------------------------------------
//    W- |/////|/////|/////|  -  |  +  |  +
// -------------------------------------------
//    Z  |/////|/////|/////|/////|  -  |  +
// -------------------------------------------
//    H  |/////|/////|/////|/////| ////|  -
bool Flavour::can_be_recombined(const Flavour & other_flav) const {
  if ((this->flavour()==23 && other_flav.flavour()==23) ||
      (this->flavour()==22 && other_flav.flavour()==22) ||
      (this->flavour()==22 && other_flav.flavour()==23) ||
      (this->flavour()==23 && other_flav.flavour()==22) ||
      (this->flavour()==25 && other_flav.flavour()==81) ||
      (this->flavour()==81 && other_flav.flavour()==25) ||
      (this->flavour()==25 && other_flav.flavour()==25) ||
      (this->flavour()==25 && other_flav.flavour()==22) ||
      (this->flavour()==22 && other_flav.flavour()==25) ||
      (this->flavour()==24 && other_flav.flavour()==24) ||
      (this->flavour()==-24 && other_flav.flavour()==-24)) return false;
  else return true;
}

/// To deal with masses, I follow the following rule to recombine
/// flavour A with flavour B, giving flavour C = A.recombine(B):
/// 1) If C=A then I give flavour C the mass of particle A (and not the default)
/// 2) Else if C=B then I five flavour C the mass of particle B
/// 3) Else (if C!=A and C!=B), then I give C its default mass
///
Flavour Flavour::recombine(const Flavour & other_flav) const {

  if (!this->can_be_recombined(other_flav)) {
    cout << "ERROR! " << _pdg_id << " and " << other_flav.flavour();
    cout << " cannot be recombined!" << endl;
    cout << "A test should be made before trying to recombine 2 flavours. ";
    cout << "Calling function can_be_recombined() is highly recommended!";
    cout << endl;
    exit(-1);
  }
  
  double mA = this->m();
  double mB = other_flav.m();
  if (this->flavour()==81) {
      switch(other_flav.flavour()) {
      case 81: return Flavour(81,mA);
      case 23: return Flavour(81,mA);
      case 22: return Flavour(81,mA);
      case 24: return Flavour(81,mA);
      case -24: return Flavour(81,mA);
      }
  } else if (this->flavour()==23) {
      switch(other_flav.flavour()) {
        case 81: return Flavour(81,mB);
        case 24: return Flavour(24,mB);
        case -24: return Flavour(-24,mB);
        case 25: return Flavour(23,mB);
      }
  } else if (this->flavour()==22) {
      switch(other_flav.flavour()) {
        case 81: return Flavour(81,mB);
        case 24: return Flavour(24,mB);
        case -24: return Flavour(-24,mB);
      }
  } else if (this->flavour()==24) {
      switch(other_flav.flavour()) {
        case 81: return Flavour(81,mB);
        case 23: return Flavour(24,mA);
        case 22: return Flavour(24,mA);
        case -24: return Flavour(23);
        case 25: return Flavour(24,mA);
      }
  } else if (this->flavour()==-24) {
      switch(other_flav.flavour()) {
        case 81: return Flavour(81,mB);
        case 23: return Flavour(-24,mA);
        case 22: return Flavour(-24,mA);
        case 24: return Flavour(23);
        case 25: return Flavour(-24,mA);
      }
  } else if (this->flavour()==25) {
      switch(other_flav.flavour()) {
        case 23: return Flavour(23,mB);
        case 24: return Flavour(24,mB);
        case -24: return Flavour(-24,mB);
      }
  } 

  //We should not reach this point
  cout << "ERROR: A flavour is missing in function recombine() " << endl;
  cout << "while recombining " << this->flavour() << " and " << other_flav.flavour() << endl;
  cout << "in Flavour.cc. " << endl;
  exit(-1);

}

/// translate a fastjet user_index into this flavour
/// user_index = flavour_index+jet_index
/// jet_index is in [0x00000000, 0x0000FFFF]
/// flavour_index is in [0xFFFF0000, 0xFFFFFFFF] and it is just pdg id of
/// a particle <<-shifted by 16 bits
void Flavour::set_flavour_from_user_index(int user_index) {
  *(this)  = Flavour(user_index >> 16);
}


/// returns an integer in which the flavour information has been
/// added (top two bytes) to the non-flavour part of the supplied index.
///
/// This is also what should be taken as the user index for particle
/// jet_index;
/// Warning: removes previous flavours information.
int Flavour::flavour_plus_index(int jet_index) const {
  return (jet_index & flavour_annihilator) + (_pdg_id << 16);
}



int Flavour::is_closer(const Flavour & fA, const Flavour & fB) {
  if (fA.is_of_the_same_kind(fB)) {
    if (fA.flavour()==fB.flavour()) return 0;
    else if (this->flavour()==fA.flavour()) return 1;
    else if (this->flavour()==fB.flavour()) return -1;
    else return 0;
  } else { //fA and fB are very different
    if (this->is_of_the_same_kind(fA)) return 1;
    else if (this->is_of_the_same_kind(fB)) return -1;
    else return 0; //should not happen with only 2 different kinds up to now
  }
}



void Flavour::_set_all_flags_to_false() {
  _is_a_parton = false;
  _is_a_Z = false;
  _is_a_Wp = false;
  _is_a_Wm = false;
  _is_a_gamma= false;
  _is_a_H = false;
}


void delete_flavour(int & user_index) {
  user_index = user_index & flavour_annihilator;
}
