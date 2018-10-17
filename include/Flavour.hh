#ifndef _FLAVOUR_HH_
#define _FLAVOUR_HH_

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
// $Id: Flavour.hh 2255 2014-11-01 11:50:12Z sapeta $
//-------------------------------------------------------ENDHEADER----


#include <vector>

//useful for conversions between flavours and fastjet's user_indexes
const int flavour_annihilator = 0x0000FFFF;

/// class storing flavour information
/// the treatment of heavy flavours is incomplete and not tested
class Flavour {
public:
  /// Default constructor
  Flavour(){}

  /// flav is just pdg id of particle.
  /// For partons use the value 81
  Flavour(int pdg_id);
  /// legacy ctor that takes one of q,b,t,Z,W (not all correctly treated)
  Flavour(char flav);

  /// In order to use the generated mass and not the default one
  Flavour(int pdg_id, double m);
  /// legacy ctor that takes one of q,b,t,Z,W (not all correctly treated)
  Flavour(char flav, double m);

  /// access routines
  int flavour() const {return _pdg_id;}
  int pdg_id()  const {return _pdg_id;}
  bool is_a_parton() const {return _is_a_parton;}
  bool is_a_vector_boson() const {return _is_a_Z || _is_a_Wp || _is_a_Wm || 
                                         _is_a_gamma || _is_a_H;}
  bool is_a_Z() const {return _is_a_Z;}
  bool is_a_W() const {return _is_a_Wp || _is_a_Wm;}
  bool is_a_Wp() const {return _is_a_Wp;}
  bool is_a_Wm() const {return _is_a_Wm;}
  bool is_a_gamma() const {return _is_a_gamma;}
  bool is_a_H() const {return _is_a_H;}
  double m() const {return _mass;}

  /// decide if 2 flavours can be recombined
  /// should always be called to test if the recombine function can be used
  bool can_be_recombined (const Flavour & other_flav) const;
  /// recombine 2 flavours
  Flavour recombine(const Flavour & other_flav) const;

  /// translate a fastjet user_index into this flavour
  /// user_index = flavour_index+jet_index
  /// jet_index is in [0x00000000, 0x0000FFFF]
  /// flavour_index is in [0xFFFF0000, 0xFFFFFFFF] and it is just pdg id of
  /// a particle <<-shifted by 16 bits
  void set_flavour_from_user_index(int user_index);

  /// returns an integer in which the flavour information has been
  /// added (top two bytes) to the non-flavour part of the supplied index.
  ///
  /// This is also what should be taken as the user index for particle
  /// jet_index;
  /// Warning: removes previous flavours information.
  int flavour_plus_index(int jet_index) const;

  /// returns true if this flavour and the other one are both a parton or
  /// are both a vector boson
  bool is_of_the_same_kind(const Flavour & other_flav) const {
    return (this->is_a_parton())==(other_flav.is_a_parton());
  }

  /// decides which of fA and fB is "closer" to this flavour. For instance:
  /// -> t is closer to b than W
  /// -> W is closer to Z than q
  /// -> None of W+ and W- are closer to Z
  /// [NB: Useful to decide which user_index to propagate in TreeLevel]
  ///
  /// return:
  /// 0 if cannot be decided
  /// 1 if fA is the nearest
  /// -1 if fB is the nearest
  int is_closer(const Flavour & fA, const Flavour & fB);

private:
  int _pdg_id;
  bool _is_a_parton, _is_a_Z, _is_a_Wp, _is_a_Wm, _is_a_gamma, _is_a_H;
  double _mass;
  void _set_all_flags_to_false();
  void _set_flavour(int pdg_id);
  void _set_flavour(int pdg_id, double m);
};


/// delete info on flavour in user_index
void delete_flavour(int & user_index);


#endif //_FLAVOUR_HH_
