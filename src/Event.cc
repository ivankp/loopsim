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
// $Id: Event.cc 2252 2014-10-28 16:58:54Z sapeta $
//-------------------------------------------------------ENDHEADER----

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include "loopsim/Event.hh"

namespace loopsim {

using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
LSParticle::LSParticle(double px, double py, double pz, double E, int flavour,
		       vector<fj::PseudoJet> decay_products) {
  this->reset(px,py,pz,E);
  _f=Flavour(flavour,this->m());
  _decay_products = decay_products;
}

//----------------------------------------------------------------------
LSParticle::LSParticle(const fj::PseudoJet & p, int flavour,
	   vector<fj::PseudoJet> decay_products) {
  this->reset(p);
  _f=Flavour(flavour,this->m());
  _decay_products = decay_products;
}

//----------------------------------------------------------------------
void LSParticle::reset_4vector(const fj::PseudoJet & p) {
  int user_index = this->user_index();
  int cluster_hist_index = this->cluster_hist_index();
  int f = this->flavour().flavour();
  this->reset(p);
  this->set_user_index(user_index);
  this->set_cluster_hist_index(cluster_hist_index);
  assert(f==this->flavour().flavour());
}


//----------------------------------------------------------------------
/// recoil the decay products, if there are some, of this particle according
/// to its old momentum oldp or the boost vector k
void LSParticle::recoil_decay_products() {

  vector<fj::PseudoJet> dp = this->decay_products();
  int ndp = dp.size();
  if (ndp==0) return; //nothing to recoil
  fj::PseudoJet tot(0,0,0,0); //old particle's momentum = sum_i dp[i]
  for (int idp=0; idp<ndp; idp++) tot += dp[idp];
  fj::PseudoJet pref_init = PtYPhiM(0,tot.rap(),0,1.0);
  fj::PseudoJet pref_end = PtYPhiM(0,this->rap(),0,1.0);
  fj::PseudoJet p0 = PtYPhiM(tot.perp(),0,tot.phi(),tot.m());
  fj::PseudoJet p1 = PtYPhiM(this->perp(),0,this->phi(),this->flavour().m());
  double Ek = p1.E()+p0.E();
  double C = pow(p1.px()-p0.px(),2)+pow(p1.py()-p0.py(),2);
  C = C/pow(Ek,2);
  double G = 2.0/(1.0+C);
  double kx = G*(p1.px()-p0.px());
  double ky = G*(p1.py()-p0.py());
  fj::PseudoJet k = fj::PseudoJet(kx,ky,0,Ek);

  //transversally boost the decay products according to k
  for (int i=0; i<ndp; i++) {
    dp[i].unboost(pref_init);
    dp[i].boost(k);
    dp[i].boost(pref_end);
  }
  //save the new decay products in this particle
  this->set_decay_products(dp);
}


//----------------------------------------------------------------------
vector<LSParticle> Event::Zbosons() const {
  vector<LSParticle> Zbosons;
  for (unsigned j=0; j<particles.size(); j++) {
    if (particles[j].flavour().is_a_Z()) {
      Zbosons.push_back(particles[j]);
    }
  }
  return Zbosons;
}

//----------------------------------------------------------------------
vector<LSParticle> Event::Wbosons() const {
  vector<LSParticle> Wbosons;
  for (unsigned j=0; j<particles.size(); j++) {
    if (particles[j].flavour().is_a_W()) {
      Wbosons.push_back(particles[j]);
    }
  }
  return Wbosons;
}

//----------------------------------------------------------------------
vector<LSParticle> Event::Hbosons() const {
  vector<LSParticle> Hbosons;
  for (unsigned j=0; j<particles.size(); j++) {
    if (particles[j].flavour().is_a_H()) {
      Hbosons.push_back(particles[j]);
    }
  }
  return Hbosons;
}


//----------------------------------------------------------------------
vector<LSParticle> Event::partons() const {
  vector<LSParticle> partons;
  for (unsigned j=0; j<particles.size(); j++) {
    if (particles[j].flavour().is_a_parton()) {
      partons.push_back(particles[j]);
    }
  }
  return partons;
}


//----------------------------------------------------------------------
PseudoJet Event::sum() const {
  PseudoJet total(0,0,0,0);
  for (unsigned i = 0; i < particles.size(); i++) {
    total += particles[i];
  }
  return total;
}

//----------------------------------------------------------------------
double Event::HT() const {
  double ht = 0;
  for (unsigned i = 0; i < particles.size(); i++) {
    ht += particles[i].perp();
  }
  return ht;
}


void print_jet(const PseudoJet & p) {
  cout << p.px() << "  " << p.py() << "  " << p.pz() << "  " << p.E() << endl;
}

void print_jet_ppvars(const PseudoJet & p) {
  cout << p.perp() << "  " << p.rap() << "  " << p.phi() << "  " << p.m() << endl;
}

void print_event(const Event & ev, const bool cart,
                 const bool with_decay_products, std::ostream* ostr) {
  PseudoJet tot(0,0,0,0);
  *ostr << "weight = " << ev.weight << endl;
  for (unsigned i=0; i<ev.particles.size(); i++) {
    tot += ev.particles[i];
    *ostr << "particle " << i << " (pdg id="
         << ev.particles[i].flavour().flavour() << ") : ";
    if (cart) {
        ostr->width(9);
	*ostr << ev.particles[i].px() << "  "
	     << ev.particles[i].py() << "  ";
	*ostr << ev.particles[i].pz() << "  "
	     << ev.particles[i].E() << endl;
    } else {
        ostr->width(9);
        *ostr << ev.particles[i].perp() << "  "
             << ev.particles[i].rap()  << "  ";
        *ostr << ev.particles[i].phi()  << "  "
             << ev.particles[i].m()    <<  endl;
    }
    vector<fj::PseudoJet> dp = ev.particles[i].decay_products();
    int ndp = dp.size();
    fj::PseudoJet totdp(0,0,0,0);

    if (with_decay_products) {
      for (int jdp=0; jdp<ndp; jdp++) {
        totdp += dp[jdp];
        *ostr << "     decay product " << jdp+1 << ":  ";
        if (cart) {
          ostr->width(9);
          *ostr << dp[jdp].px() << "  "
               << dp[jdp].py() << "  ";
          *ostr << dp[jdp].pz() << "  "
               << dp[jdp].E() << endl;
        } else {
          ostr->width(9);
          *ostr << dp[jdp].perp() << "  "
               << dp[jdp].rap()  << "  ";
          *ostr << dp[jdp].phi()  << "  "
               << dp[jdp].m()    <<  endl;
        }
      }
      if (cart) {
        *ostr << "     sum of decay products: ";
        *ostr << totdp.px() <<"  "<< totdp.py() <<"  "<< totdp.pz();
        *ostr << "  " << totdp.E() << endl;
      }
    }
  }
  if (cart) {
      *ostr << "total momentum : ";
      *ostr << tot.px() <<"  "<< tot.py() <<"  "<< tot.pz() << "  " << tot.E();
  } else {
      *ostr << "pt_tot:  " << tot.perp() << "  HT:   " << ev.HT();
  }
  *ostr << endl;

}

}

