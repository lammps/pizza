# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# pair tool

oneline = "Compute LAMMPS pairwise energies"

docstr = """
p = pair("lj/charmm/coul/charmm")   create pair object for specific pair style

  available styles: lj/cut, lj/cut/coul/cut, lj/charmm/coul/charmm

p.coeff(d)			    extract pairwise coeffs from data object
p.init(cut1,cut2,...)		    setup based on coeffs and cutoffs

  init args are specific to pair style:
    lj/cut = cutlj
    lj/cut/coul/cut = cutlj,cut_coul (cut_coul optional)
    lj/charmm/coul/charmm = cutlj_inner,cutlj,cutcoul_inner,cut_coul
      (last 2 optional)
      
e_vdwl,e_coul = p.single(rsq,itype,jtype,q1,q2,...)   compute LJ/Coul energy

  pairwise energy between 2 atoms at distance rsq with their attributes
  args are specific to pair style:
    lj/cut = rsq,itype,jtype
    lj/cut/coul/cut = rsq,itype,jtype,q1,q2
    lj/charmm/coul/charmm = rsq,itype,jtype,q1,q2
"""

# History
#   8/05, Steve Plimpton and Paul Crozier (SNL): original version
#   9/05, Paul Crozier (SNL): added lj/cut and lj/cut/coul/cut

# ToDo list

# Variables

# Imports and external programs

from math import sqrt

# Class definition

class pair:

  # --------------------------------------------------------------------

  def __init__(self,style):

    # map generic methods to style-specific methods

    if style == "lj/cut":
      self.coeff_func = self.coeff_lj_cut
      self.init_func = self.init_lj_cut
      self.single_func = self.single_lj_cut
    elif style == "lj/cut/coul/cut":
      self.coeff_func = self.coeff_lj_cut_coul_cut
      self.init_func = self.init_lj_cut_coul_cut
      self.single_func = self.single_lj_cut_coul_cut
    elif style == "lj/charmm/coul/charmm":
      self.coeff_func = self.coeff_lj_charmm_coul_charmm
      self.init_func = self.init_lj_charmm_coul_charmm
      self.single_func = self.single_lj_charmm_coul_charmm
    else:
      raise StandardError, "this pair style not yet supported"

  # --------------------------------------------------------------------
  # generic coeff method

  def coeff(self,data):
    self.coeff_func(data)

  # --------------------------------------------------------------------
  # generic init method, as many args as needed

  def init(self,*list):
    self.init_func(list)

  # --------------------------------------------------------------------
  # generic single method, as many args as needed

  def single(self,*list):
    return self.single_func(list)

  # --------------------------------------------------------------------
  # --------------------------------------------------------------------
  # lj/cut methods

  def coeff_lj_cut(self,data):
    epsilon = data.get("Pair Coeffs",2)
    sigma = data.get("Pair Coeffs",3)
    ntypes = len(epsilon)
    
    self.lj3 = []
    self.lj4 = []     
    for i in xrange(ntypes):
      self.lj3.append(ntypes * [0])
      self.lj4.append(ntypes * [0])
      for j in xrange(ntypes):
        epsilon_ij = sqrt(epsilon[i]*epsilon[j])
        sigma_ij = sqrt(sigma[i]*sigma[j])
        self.lj3[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,12.0);
        self.lj4[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,6.0);
  
  # --------------------------------------------------------------------
  # args = cutlj

  def init_lj_cut(self,list):
    cut_lj = list[0]
    self.cut_ljsq = cut_lj*cut_lj

  # --------------------------------------------------------------------
  # args = rsq,itype,jtype
  
  def single_lj_cut(self,list):
    rsq = list[0]
    itype = list[1]
    jtype = list[2]
    
    r2inv = 1.0/rsq
    
    if rsq < self.cut_ljsq:
      r6inv = r2inv*r2inv*r2inv
      eng_vdwl = r6inv*(self.lj3[itype][jtype]*r6inv-self.lj4[itype][jtype])
    else: eng_vdwl = 0.0
      
    return eng_vdwl 

  # --------------------------------------------------------------------
  # --------------------------------------------------------------------
  # lj/cut/coul/cut methods

  def coeff_lj_cut_coul_cut(self,data):
    epsilon = data.get("Pair Coeffs",2)
    sigma = data.get("Pair Coeffs",3)
    ntypes = len(epsilon)
    
    self.lj3 = []
    self.lj4 = []     
    for i in xrange(ntypes):
      self.lj3.append(ntypes * [0])
      self.lj4.append(ntypes * [0])
      for j in xrange(ntypes):
        epsilon_ij = sqrt(epsilon[i]*epsilon[j])
        sigma_ij = sqrt(sigma[i]*sigma[j])
        self.lj3[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,12.0);
        self.lj4[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,6.0);
  
  # --------------------------------------------------------------------
  # args = cutlj, cut_coul (cut_coul optional)

  def init_lj_cut_coul_cut(self,list):
    self.qqr2e = 332.0636  			# convert energy to kcal/mol
    cut_lj = list[0]
    self.cut_ljsq = cut_lj*cut_lj
    
    if len(list) == 1: cut_coul = cut_lj
    else: cut_coul = list[1]
    self.cut_coulsq = cut_coul*cut_coul

  # --------------------------------------------------------------------
  # args = rsq,itype,jtype,q1,q2
  
  def single_lj_cut_coul_cut(self,list):
    rsq = list[0]
    itype = list[1]
    jtype = list[2]
    q1 = list[3]
    q2 = list[4]
    
    r2inv = 1.0/rsq
    
    if rsq < self.cut_coulsq: eng_coul = self.qqr2e * q1*q2*sqrt(r2inv)
    else: eng_coul = 0.0
    
    if rsq < self.cut_ljsq:
      r6inv = r2inv*r2inv*r2inv
      eng_vdwl = r6inv*(self.lj3[itype][jtype]*r6inv-self.lj4[itype][jtype])
    else: eng_vdwl = 0.0
      
    return eng_coul,eng_vdwl   
  
  # --------------------------------------------------------------------
  # --------------------------------------------------------------------
  # lj/charmm/coul/charmm methods

  def coeff_lj_charmm_coul_charmm(self,data):
    epsilon = data.get("Pair Coeffs",2)
    sigma = data.get("Pair Coeffs",3)
    ntypes = len(epsilon)
    
    self.lj3 = []
    self.lj4 = []     
    for i in xrange(ntypes):
      self.lj3.append(ntypes * [0])
      self.lj4.append(ntypes * [0])
      for j in xrange(ntypes):
        epsilon_ij = sqrt(epsilon[i]*epsilon[j])
        sigma_ij = 0.5 * (sigma[i] + sigma[j])
        self.lj3[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,12.0);
        self.lj4[i][j] = 4.0 * epsilon_ij * pow(sigma_ij,6.0);
  
  # --------------------------------------------------------------------
  # args = cutlj_inner,cutlj,cutcoul_inner,cut_coul (last 2 optional)

  def init_lj_charmm_coul_charmm(self,list):
    self.qqr2e = 332.0636  			# convert energy to kcal/mol
    cut_lj_inner = list[0]
    cut_lj = list[1]

    self.cut_lj_innersq = cut_lj_inner*cut_lj_inner
    self.cut_ljsq = cut_lj*cut_lj
    
    if len(list) == 2:
      cut_coul_inner = cut_lj_inner
      cut_coul = cut_lj
    else:
      cut_coul_inner = list[2]
      cut_coul = list[3]
    self.cut_coul_innersq = cut_coul_inner*cut_coul_inner
    self.cut_coulsq = cut_coul*cut_coul

    self.denom_lj = (self.cut_ljsq-self.cut_lj_innersq) *  \
      (self.cut_ljsq-self.cut_lj_innersq) *  \
      (self.cut_ljsq-self.cut_lj_innersq);
    self.denom_coul = (self.cut_coulsq-self.cut_coul_innersq) *  \
      (self.cut_coulsq-self.cut_coul_innersq) *  \
      (self.cut_coulsq-self.cut_coul_innersq);

  # --------------------------------------------------------------------
  # args = rsq,itype,jtype,q1,q2
  
  def single_lj_charmm_coul_charmm(self,list):
    rsq = list[0]
    itype = list[1]
    jtype = list[2]
    q1 = list[3]
    q2 = list[4]
    
    r2inv = 1.0/rsq
    
    if rsq < self.cut_coulsq:
      eng_coul = self.qqr2e * q1*q2*sqrt(r2inv)
      if rsq > self.cut_coul_innersq:
        switch1 = (self.cut_coulsq-rsq) * (self.cut_coulsq-rsq) *  \
                  (self.cut_coulsq + 2.0*rsq - 3.0*self.cut_coul_innersq) /  \
                  self.denom_coul
        eng_coul *= switch1
    else: eng_coul = 0.0
    
    if rsq < self.cut_ljsq:
      r6inv = r2inv*r2inv*r2inv
      eng_vdwl = r6inv*(self.lj3[itype][jtype]*r6inv-self.lj4[itype][jtype])
      if rsq > self.cut_lj_innersq:
        switch1 = (self.cut_ljsq-rsq) * (self.cut_ljsq-rsq) *  \
                  (self.cut_ljsq + 2.0*rsq - 3.0*self.cut_lj_innersq) /  \
                  self.denom_lj
        eng_vdwl *= switch1
    else: eng_vdwl = 0.0
      
    return eng_coul,eng_vdwl
