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
// $Id: TreeLevel.cc 947 2012-03-15 10:10:26Z sapeta $
//-------------------------------------------------------ENDHEADER----

#include "TreeLevel.hh"
#include "ReshufflePlugin.hh"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>

using namespace std;
using namespace fastjet;


void print_history(const ClusterSequence::history_element * hist) {
  cout << hist->parent1 << "  " << hist->parent2 << "  " << hist->child;
  cout << "  " << hist->jetp_index << "  " << hist->dij << "  " << endl;
}

//construct a pseudojet using M instead of E
PseudoJet PxPyPzM(const double px, const double py, 
		      const double pz, const double M) {
  double E = sqrt(px*px + py*py + pz*pz + M*M);
  return PseudoJet(px,py,pz,E);
}

//swap 2 pseudojet pointers
// GPS: how come not just use the stl swap??
void pointer_swap(const PseudoJet * p1, const PseudoJet * p2) {
  const PseudoJet * temp = p1;
  p1 = p2;
  p2 = temp;
}


TreeLevel::TreeLevel(const Event & ev, int nborn, double R, bool spread_virtuals) {
  _no_loop_event = ev;
  _nparticles = _no_loop_event.particles.size();
  _nborn = nborn;
  _spread_virtuals = spread_virtuals;

  //be sure that all flavours are defined
  //assert(_nparticles==int(_no_loop_event.flavours.size()));
  //assert(_nparticles==_no_loop_event.flavours.size());

  // make sure input event is OK
  verify_no_missing_pt(ev, "TreeLevel constructor (input event)");

  //find particles' status
  //place never virtual particles first in the event
  init(R);

  //resize in order to obtain all possible numbers of loops, 
  //from 0 to _nparticles
  _all_events.resize(_nparticles+1);

  //fill _all_events
  //Event my_event;
  //my_event.weight = ev.weight;
  //generate_all_loops(0,0,my_event);
  vector<int> virtuals;
  generate_all_relevant_loops(0,0,virtuals);

}


const vector<vector<Event> > & TreeLevel::extract_all_events(){
  return _all_events;
}

const vector<Event> & TreeLevel::extract_all_events(int iloops) {
  return _all_events[iloops];
}

const Event & TreeLevel::extract_original_event() {
  return _all_events[0][0];
}

const Event & TreeLevel::extract_event(int i, int j) {
  if (i>_nvirtuals || j>=int(_all_events[i].size())) {
    cout << "ERROR: indexes i and j out of range in ";
    cout << "TreeLevel::extract_event(i,j)!" << endl;
    exit(-1);
  }
  return _all_events[i][j];
}


void TreeLevel::print_all_events() {
  for (int i=0; i<_nvirtuals+1; i++) {
    cout << endl << "number of loops = " << i << endl << endl;
    for (unsigned j=0; j<_all_events[i].size(); j++) {
      cout << "all_events[" << i <<"][" << j << "]: " << endl;
      print_event(_all_events[i][j]);
      cout << endl;
    }
  }
}


void TreeLevel::init(double R) {
  _status.resize(_nparticles,undef);

  //init particles such that their user_index is equal to their place in the
  //particles' vector
  //Moreover, put some info on flavour (see Flavour.hh for more info)
  for (int i=0; i<_nparticles; i++) {
    // GPSToCheck: rather than added i+f_index, use translate(i)
    // [check this elsewhere, and make sure results are unchanged...]
    int user_index=_no_loop_event.particles[i].flavour().flavour_plus_index(i);
    _no_loop_event.particles[i].set_user_index(user_index);
  }

  //Cluster the event
  MEAFRecombiner rec;
  FlavourPlugin fp(R);
  JetDefinition jet_def(&fp);
  jet_def.set_recombiner(&rec);
  _cs.reset(new ClusterSequence(_no_loop_event.particles,jet_def));
  
  //sort the branchings into decreasing kt distance
  vector<AuxHist> aux_vector;
  for (unsigned i = 0; i < _cs->history().size(); i++) {
    AuxHist aux_element(*_cs,i);
    aux_vector.push_back(aux_element);
    //cout << i << " : " << "ktij = " << aux_element.aux_distance << "  ";
    //print_history(&_cs->history()[i]);
  }
  // sorts the aux vector into _decreasing_ kt-distance order
  // (cf definition of operator< below)
  sort(aux_vector.begin(), aux_vector.end());

  //find the status of all particles in the event
  //We go through all the clustering history until all particles have been
  //checked (condition: ip = _nparticles)
  int ip = 0, aux_index = 0;
  const ClusterSequence::history_element * hist;
  const PseudoJet * p1,* p2;
  while (ip != _nparticles) {
    int jet_hist_index =  aux_vector[aux_index].hist_index;
    // GPS: the following line is correct for cam, but not correct in
    // general (e.g. if we change the jet alg one day)
    //int jet_user_index = _cs->jets()[jet_hist_index].user_index();
    // ---------
    // GPS: the general alternative which accounts for the
    //      fact that jet indices and history indices do 
    //      not always coincide
    hist = &(_cs->history()[jet_hist_index]);

    int jet_user_index = _cs->jets()[hist->jetp_index].user_index();
    aux_index++;
    if (hist->parent1 < 0) {
      cout << "About to crash on negative hist->parent1 = "<< hist->parent1 << endl;
      cout << "_no_loop_event is:\n";
      print_event(_no_loop_event);
      throw fastjet::Error("Invalid hist->parent1 in TreeLevel::init");
    }

    //assert(hist->parent1 >= 0);
    // GPS: again care needed with history and jet indices
    p1 = &(_cs->jets()[_cs->history()[hist->parent1].jetp_index]);
    //if p1 clusters with the beam...
    if (hist->parent2<0) {
      fill_status(ip,jet_user_index,p1);
      continue;
    }
    //if p1 clusters with p2...
    // GPS: again care needed with history and jet indices
    p2 = &(_cs->jets()[_cs->history()[hist->parent2].jetp_index]);
    fill_status(ip,jet_user_index,p1,p2);
  }

  //assert(ip == _nparticles);

  //finish the initialisation according to the status found
  //place never virtual particles first in the _no_loop_event
  //store the status in particles' user_index.
  finish_init();

}


void TreeLevel::fill_status(int & ip, int jet_index, const PseudoJet * p1, 
			    const PseudoJet * p2) {

  //we delete info on flavour
  delete_flavour(jet_index);
  int index1 = p1->user_index();
  delete_flavour(index1);
  int status1 = _status[index1];
  //case where p1 clusters with the beam
  if (p2==NULL) {
    switch(status1) {
    case nc: return; //cannot become anything else
    case born: return; //already a born particle
    case undef: 
      if (ip<_nborn) _status[index1] = born; //becomes a born particle
      else _status[index1] = beam; //need to see later if it can become virtual
      
      // // pseudocode alternative to account for configurations like:
      // // (q1->q1+W) recoiling against (q2), which has the problem that
      // // when nborn = 1 the (q1+W) object is "the born" particle. Then
      // // q2 can get looped leaving a diagram in which q1 recoils
      // // against W, which is a radical and unphysical change of
      // // topology. To work around this, one option is to look at the
      // // flavour of born 1, and if it's a parton declare that there has to be
      // // a second born particle. 
      // //
      // // [Revisit this is we ever go to single top production!!!!]
      // if (ip<_nborn) {
      // 	_status[index1] = born; //becomes a born particle
      // 	if (_nborn == 1 && flavour_is_parton(index1)) _nborn = 2;
      // } else {
      // 	_status[index1] = beam; //need to see later if it can become virtual
      // }


      ip++;
      return;
    default:
    //status1 should be nc,born, or undef (ie <-1)  because it cannot cluster 
    //twice with the beam, or cluster with another harder particle later
      cout << "it is not possible for status1 to be " << status1;
      cout << " in function fill_status(...)!" << endl;
      cout << "There is a problem!" << endl;
      exit(-1);
    }
  }

  //case where p1 clusters with p2.
  int index2 = p2->user_index(); 
  delete_flavour(index2); //delete info on flavour
  int status2 = _status[index2];

  //We want index1 = jet_index so that the first particle is the one
  //whose index propagates
  if (jet_index == index2) {
    swap(status1,status2);
    swap(index1,index2);
  }
  //necessarily, p1 cannot become virtual.
  switch(status1) {
  case undef:
    if (ip<_nborn) _status[index1] = born;
    else _status[index1] = nc;
    ip++;
    break;
  case born: break;
  case nc: break;
  default: 
    //p1 was previously declared to cluster with the beam or another particle
    _status[index1] = nc;
    break;
  }
  //necessarily, p2 is declared as clustering with p1, except if 
  //it can be (or already is) a born particle, 
  //or if it was previously declared nc.
  switch(status2) {
  case undef:
    if (ip<_nborn) _status[index2] = born;
    else _status[index2] = index1;
    ip++;
    break;
  case born: break;
  case nc: break;
  default:
    cout << "Problem! status2 = " << status2 << " in function fill_status()";
    cout << " whereas it should not!" << endl;
    exit(-1);
  }

}


void TreeLevel::finish_init() {

  //change particles ordering in the event: place never virtual particles first
  Event ev;
  ev.weight = _no_loop_event.weight;
  vector<int> position(_nparticles); //position of real particles in the new event
  int pos = 0;
  for (int i=0; i<_nparticles; i++) {
    if (_status[i]>=-1) continue;
    ev.particles.push_back(_no_loop_event.particles[i]);
    ev.particles[pos].set_user_index(_status[i]);
    //ev.flavours.push_back(_no_loop_event.flavours[i]);
    position[i]=pos;
    pos++;
  }
  _nreals = pos;


  //fill _pseudo_virtuals and finish filling ev.particles
  for (int i=0; i<_nparticles; i++) {
    if (_status[i]<-1) continue; //particle i will never become virtual
    _pseudo_virtuals.push_back(pos);
    //as particles position changed, you have to reinit the status
    //and the flavours
    int status = _status[i]==beam ? beam : position[_status[i]];
    ev.particles.push_back(_no_loop_event.particles[i]);
    ev.particles[pos].set_user_index(status);
    //ev.flavours.push_back(_no_loop_event.flavours[i]);
    position[i] = pos;
    pos++;
  }
  _nvirtuals = _pseudo_virtuals.size();

  // GPS now create a plugin jet definition that will allow us simply
  // to reshuffle the existing cluster sequence
  _reshuffle_plugin.reset(new ReshufflePlugin(_cs.get(), position));
  _new_jet_def = _reshuffle_plugin.get();

  // get the new, reshuffled cluster sequence
  _new_cs.reset(new ClusterSequence(ev.particles, _new_jet_def));

  // and set info that allows the event particles to be referenced
  // within the cluster sequence
  for (int i=0; i<_nparticles; i++) {ev.particles[i].set_cluster_hist_index(i);}
  
  // finally, copy the reordered event across to the _no_loop_event
  for (int i=0; i<_nparticles; i++) {
    _no_loop_event.particles[i] = ev.particles[i];
    //_no_loop_event.flavours[i] = ev.flavours[i];
  }

  /*
  // reinit _status vector, but I think this step will be cancelled in the
  // final version, because _status vector is now useless
  // (all relevant information is stored in particles' user_index...)
  for (int i=0; i<_nparticles; i++) {
    _status[i] = ev.particles[i].user_index();
  }
  */

}



void TreeLevel::print_status() {
  /*
  cout << "particles' old status, with user's ordering: " << endl;
  for (int i=0; i<_nparticles; i++) {
    cout << "particle " << i << " : ";
    switch(_status[i]) {
    case undef:
      cout << "undef" << endl;
      break;
    case born:
      cout << "born" << endl;
      break;
    case nc:
      cout << "nc" << endl;
      break;
    case beam:
      cout << "beam" << endl;
      break;
    default:
      cout << _status[i] << endl;
    }
  }
  cout << endl;
  */

  //  cout << "new status after reordering: " << endl;
  for (int i=0; i<_nparticles; i++) {
    cout << "particle " << i << " : ";
    cout << _no_loop_event.particles[i].flavour().flavour() << " : ";
    switch(_no_loop_event.particles[i].user_index()) {
    case undef:
      cout << "undef" << endl;
      break;
    case born:
      cout << "born" << endl;
      break;
    case nc:
      cout << "nc" << endl;
      break;
    case beam:
      cout << "beam" << endl;
      break;
    default:
      cout << _no_loop_event.particles[i].user_index() << endl;
    }
  }
  cout << endl;
}


bool operator<(const AuxHist & aux1, const AuxHist & aux2) {
  return aux1.aux_distance > aux2.aux_distance;
}



void TreeLevel::generate_all_relevant_loops(int iloops, int ip,
				 vector<int> virtuals) {

  if (ip==_nvirtuals) {
    if (_spread_virtuals) {
      _all_events[iloops].push_back(make_virtual_spread(virtuals));
    } else {
      _all_events[iloops].push_back(make_virtual(virtuals));
    }
    return;
  }

  ip++;


  //case where particle number ip is real: we just continue...
  generate_all_relevant_loops(iloops,ip,virtuals);

  //case where particle number ip is virtual: we put it in the virtuals vector
  virtuals.push_back(_pseudo_virtuals[ip-1]);
  generate_all_relevant_loops(iloops+1,ip,virtuals);

}


//In all this function, we want to make a copy of _no_loop_event, but
//filling it only with recombined particles
Event TreeLevel::make_virtual(const vector<int> & virtuals) {

  //cout << "Hey there " << virtuals.size() << endl;

  int iloops = virtuals.size();

  //first, copy never virtual particles
  Event evt;
  //evt.k = _no_loop_event.k;
  evt.weight = iloops%2==0 ? _no_loop_event.weight : -_no_loop_event.weight;
  //cout << "----- " << _nreals << " " << _no_loop_event.particles.size() << " " << virtuals.size() << endl;
  for (int i=0; i<_nreals; i++) {
    evt.particles.push_back(_no_loop_event.particles[i]);
    //evt.flavours.push_back(_no_loop_event.flavours[i]);
  }

  //Recombine virtual particles with the good particles and make
  //the sum of all particles whose status is "beam"
  int nbeam = 0; //number of particles which cluster with the beam
  PseudoJet vbeam(0,0,0,0);
  vector<bool> is_recombined(_nparticles,false);
  for (int i=0; i<iloops; i++) {
    const PseudoJet * pi = &(_no_loop_event.particles[virtuals[i]]);
    //cout << i << " " << pi->cluster_hist_index() << endl;
    is_recombined[virtuals[i]] = true;
    if (pi->user_index()==beam) {
      nbeam++;
      vbeam += *pi;
      continue;
    }
    evt.particles[pi->user_index()] += *pi;
    //cout << "                " << pi->user_index() << " " << pi->cluster_hist_index() << endl;
    Flavour virtual_flavour = _no_loop_event.particles[virtuals[i]].flavour();
    //the compatibility between the 2 flavours was already checked
    //during the clustering
    evt.particles[pi->user_index()].set_flavour(evt.particles[pi->user_index()].flavour().recombine(virtual_flavour));
    is_recombined[pi->user_index()] = true;
  }

  //Finish filling the event with all the other real particles,
  //i.e. particles that are real for this event but which are
  //virtual in other events generated from _no_loop_event.
  for (int i=_nreals; i<_nparticles; i++) {
    if (!is_recombined[i]) {
      evt.particles.push_back(_no_loop_event.particles[i]);
      //evt.flavours.push_back(_no_loop_event.flavours[i]);
    }
  }

  /*
  // below, we have the original recoil code; one can also call
  // the do_recoils routine which should be equvialent modulo
  // rounding errors
  do_recoils(evt, vbeam);
  return evt; */
  

  //if there's at least one beam particle, then we will have to make a pt boost,
  //such that we get a purely transverse event, all particles having their
  //true mass. I store rapidities in order to restablish them after the boost
  //
  //If there's no beam particle, I restore true masses to recombined particles
  int ireals = _nparticles-iloops;
  //int ireals = evt.particles.size();
  double rapidities[ireals];
  double Etot = 0; //energy for the boost
  if (nbeam!=0) {
    for (int i=0; i<ireals; i++) {
      PseudoJet * pi = &(evt.particles[i]);
      rapidities[i] = pi->rap();
      *pi = PtYPhiM(pi->perp(),0,pi->phi(),evt.particles[i].flavour().m());
      Etot += pi->E();
    }
  }
  //If there's no beam particle, then I just have to put masses of recombined
  //particles to their true mass (keeping (pt,y,phi) fixed).
  //
  // GPS after consulting with MR: this procedure is equivalent to the
  // procedure that we have when nbeam != 0, but just uses fewer
  // manipulations.
  else {
    for (int i=0; i<_nreals; i++) {
      if (is_recombined[i]) {
	PseudoJet * pi = &(evt.particles[i]);
	*pi = PtYPhiM(pi->perp(),pi->rap(),pi->phi(),evt.particles[i].flavour().m());	
      }
    }
  }
  //to finish, we make the boost
  if (nbeam!=0) {
    // GPS: by mom.cons, vbeam.px, py should be exactly - the px,py of
    // of the sum of all the other particles
    PseudoJet k(vbeam.px(),vbeam.py(),0,Etot);
    //(evt.k).push_back(k);
    for (int i=0; i<ireals; i++) {
      PseudoJet * pi = &(evt.particles[i]);
      double Mi = evt.particles[i].flavour().m();
      //We should write our own unboost function, because we do not need
      //to recalculate k.m() for each unboosted particle...
      if (ireals!=1) pi->boost(k);
      else { //k.m() can be 0
	//*pi += k;
	*pi = PtYPhiM(0,rapidities[i],0,Mi);
	break;
      }
      //restore rapidities
      *pi = PtYPhiM(pi->perp(),rapidities[i],pi->phi(),Mi);
    }
  }
  //else { //no beam particles
  //  evt.k = fastjet::PseudoJet(0,0,0,1.0);
  //}

  // make sure event is OK
  verify_no_missing_pt(evt, "make_virtual",1);

  return evt;
}


//----------------------------------------------------------------------
/// a version fo the make_virtual code that spreads the looped particle's
/// momentum among all of the parents children that were clustered at
/// smaller angle
Event TreeLevel::make_virtual_spread(const vector<int> & virtuals) {

  //cout << "Hey there " << virtuals.size() << endl;

  int iloops = virtuals.size();

  // first, create the basics of the event
  Event evt;
  //evt.k = _no_loop_event.k;
  evt.weight = iloops%2==0 ? _no_loop_event.weight : -_no_loop_event.weight;

  // establish which particles will be recombined
  vector<bool> is_virtual(_nparticles,false);
  for (int i=0; i<iloops; i++) {
    is_virtual[virtuals[i]] = true;
  }

  // now go about recombining them, storing the result in a temporary event
  PseudoJet vbeam(0,0,0,0);
  Event tmpevt;
  tmpevt.particles = _no_loop_event.particles;
  //tmpevt.flavours  = _no_loop_event.flavours;
  for (int i=0; i<iloops; i++) {
    int ivirt = virtuals[i];
    const LSParticle & pvirt = _no_loop_event.particles[ivirt];
    if (pvirt.user_index()==beam) {
      vbeam += pvirt;
    } else {
      // first figure out flavour, which goes just onto the "trunk" particle
      int itrunk = pvirt.user_index();
      tmpevt.particles[itrunk].set_flavour(tmpevt.particles[itrunk].flavour().recombine(tmpevt.particles[ivirt].flavour()));

      /*// original method: momentum entirely onto trunk particle
	tmpevt.particles[itrunk] += pvirt;*/

      // next sort out the momentum, which gets shared across all real
      // children of the current particle's partner. We make use of the
      // fact that entries in _no_loop_event have been given
      // cluster_hist_index values that allow them to be located
      // within the _new_cs; also implicitly make use of the fact that
      // particles with "children" never get made virtual, so that a
      // virtual particle's first clustering is always with "trunk+X".
      PseudoJet partner;
      assert(_new_cs->has_partner(pvirt, partner));
      vector<fj::PseudoJet> ps_children = _new_cs->constituents(partner);
      unsigned int chsize = ps_children.size();
      vector<LSParticle> children(chsize);
      for (unsigned int ich=0; ich<chsize; ich++) {
	children[ich].reset(ps_children[ich]);
	//children[ich].set_user_index(ps_children[ich].user_index());
	//children[ich].set_cluster_hist_index(ps_children[ich].cluster_hist_index());
	children[ich].set_flavour(_no_loop_event.particles[ps_children[ich].cluster_hist_index()].flavour());
	children[ich].set_decay_products(_no_loop_event.particles[ps_children[ich].cluster_hist_index()].decay_products());
      }

      // get the ptsum of the real children, using the cluster_hist_index
      // to identify the position in the original sequence
      double ptsum = 0;
      for (unsigned ich = 0; ich<chsize; ich++) {
      	if (!is_virtual[children[ich].cluster_hist_index()]) ptsum += children[ich].perp();
      }
      
      // and assign the momentum of the virtual particle to all real
      // children in proportion to their original momenta
      // (identical in _new_cs and _no_loop_event).
      for (unsigned ich = 0; ich<children.size(); ich++) {
      	int child_index = children[ich].cluster_hist_index();
      	if (!is_virtual[child_index]) {
	  // remember the momentum before reshuffling in case there are some
	  // decay products to recoil
	  LSParticle old = tmpevt.particles[child_index];
      	  tmpevt.particles[child_index] += pvirt * (children[ich].perp()/ptsum);
	  assert(old.flavour().flavour()==tmpevt.particles[child_index].flavour().flavour());
	}
      }
    }
  }
  
  // now take the temporary event and decant just the part that is not virtual
  // MR: virtuals.size()==iloops, tmpevt.particles.size()==_nparticles
  evt.particles.reserve(tmpevt.particles.size() - virtuals.size());
  //evt.flavours.reserve(tmpevt.particles.size() - virtuals.size());
  for (unsigned i = 0; i < tmpevt.particles.size(); i++) {
    if (!is_virtual[i]) {
      evt.particles.push_back(tmpevt.particles[i]);
      //evt.flavours. push_back(tmpevt.flavours[i]);
    }
  }
  assert(evt.particles.size() == tmpevt.particles.size() - virtuals.size());

  // sort out recoils
  do_recoils(evt, vbeam);

  // make sure event is OK
  verify_no_missing_pt(evt, "make_virtual_spread",1);

  return evt; 
}

//----------------------------------------------------------------------
void TreeLevel::do_recoils(Event & evt, const PseudoJet & vbeam) {
  //if there's at least one beam particle, then we will have to make a pt boost,
  //such that we get a purely transverse event, all particles having their
  //true mass. I store rapidities in order to restablish them after the boost
  //
  //If there's no beam particle, I restore true masses to recombined particles
  double rapidities[evt.particles.size()];
  double Etot = 0; //energy for the boost
  
  bool no_beam_recoil = (vbeam.px() == 0 && vbeam.py() == 0);

  
  if (no_beam_recoil) {
    // just put all particles on mass shell, while retaining their
    // pt, rapidity and azimuth; note that for some particles, this
    // may be a "null" operation
    for (unsigned i = 0; i < evt.particles.size(); i++) {
      LSParticle & pi = evt.particles[i];
      evt.particles[i].reset_4vector(PtYPhiM(pi.perp(),pi.rap(),pi.phi(),
					     pi.flavour().m()));
      evt.particles[i].recoil_decay_products();
    }
  } else {
    // figure out and record the transverse boost vector that we'll need
    // (GPS: by mom.cons, vbeam.px, py should be exactly - the px,py of
    // of the sum of all the other particles)
    //
    // get the k vector right
    PseudoJet k(vbeam.px(),vbeam.py(),0,Etot);
    for (unsigned i = 0; i < evt.particles.size(); i++) {
        LSParticle * pi = &(evt.particles[i]);
        rapidities[i] = pi->rap();
        //*pi = PtYPhiM(pi->perp(),0,pi->phi(),evt.flavours[i].m());
	// new version (work in progress)
	//pi->reset_4vector(PtYPhiM(pi->perp(),0,pi->phi(),pi->flavour().m()));
        //pi->recoil_decay_products();
        Etot += pi->E();
    }
    k += PseudoJet(0,0,0,Etot);
    if (evt.particles.size() == 1) {
      // careful treatment needed for special case of just a single particle
      // to avoid possible issues with zero-energy boosts, etc.
      LSParticle pi = evt.particles[0];
      evt.particles[0].reset_4vector(PtYPhiM(0,rapidities[0],0,pi.flavour().m()));
      evt.particles[0].recoil_decay_products();
    } else {
      // put all particles at central rapidity
      // NB: if all particles that are present have zero mass and 
      //     identical phi(), then the resulting boost vector
      //     will be ill-defined (zero mass).
      for (unsigned i = 0; i < evt.particles.size(); i++) {
        LSParticle * pi = &(evt.particles[i]);
	Flavour fi = pi->flavour();
        double Mi = fi.m();	
        //We should write our own unboost function, because we do not need
        //to recalculate k.m() for each unboosted particle...
	pi->boost(k);
	//restore rapidities (do it before recoiling the decay products!)
	pi->reset_4vector(PtYPhiM(pi->perp(),rapidities[i],pi->phi(),Mi));
	pi->recoil_decay_products();
	assert(fi.flavour()==pi->flavour().flavour());
        
        
      }
    }
    // record it for posterity...
    //(evt.k).push_back(k);
  }
}


void TreeLevel::verify_no_missing_pt(const Event & evt, const char * label, bool print_orig) const {
  // check consistency

  if (evt.sum().perp()/evt.HT() > _missing_pt_tolerance) {
    if (!(evt.particles.size()==1 && abs(evt.sum().perp()) < _zero)) {
       cerr << "WARNING: evt.sum().pt() seems non-zero = " << evt.sum().perp();
       if (label) cerr << " in " << label;
       cerr << endl;
       print_event(evt,false,true,&std::cerr);
       if (print_orig) {
         cerr << "Clustered with R = " << _cs->jet_def().R() 
	      << " and original event was " << endl;
         print_event(_no_loop_event,false, true,&std::cerr);
       }
    }
  }
}


int TreeLevel::number_of_events(int iloops) {
  if (iloops<0) {
    int res = 2;
    for (int i=1; i<_nvirtuals; i++) {
      res *= 2;
    }
    return res;
  }
  if (iloops>_nvirtuals) return -1;
  return _all_events[iloops].size();
}


/// returns the user_index associated with the recombination of pa and pb.
///
/// The user index has its flavour part set according to the
/// Flavour::recombine function. Its jet index part is equal to that
/// of the parent that is "closer" to the child, as determined by a
/// combination of flavour considerations (flavour::is_closer) and
/// kinematics.
int MEAFRecombiner::recombined_user_index(const fastjet::PseudoJet & pa,
					  const fastjet::PseudoJet & pb) const {

  Flavour fA,fB;
  int iA = pa.user_index();
  int iB = pb.user_index();
  fA.set_flavour_from_user_index(iA);
  fB.set_flavour_from_user_index(iB);
  //if not compatible, then dAB should have been set infinity in
  //function FlavourBriefJet::distance(*)
  assert(fA.can_be_recombined(fB));

  Flavour fAB = fA.recombine(fB);
  switch(fAB.is_closer(fA,fB)) {
  case 0: //none is closer: only transverse energy can decide
    if (pa.mperp2()<pb.mperp2()) return fAB.flavour_plus_index(iB);
    else return fAB.flavour_plus_index(iA);
  case 1: //fA is closer
    return fAB.flavour_plus_index(iA);
  case -1: //fB is closer
    return fAB.flavour_plus_index(iB);
  default: //not possible
    cout << "ERROR: Problem in return value for function "; 
    cout << "Flavour::is_closer(*)!" << endl;
    exit(-1);
  }

}
