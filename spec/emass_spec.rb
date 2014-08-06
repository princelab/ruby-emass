require 'spec_helper'

require 'emass'
require 'mspire/mf'
require 'mspire/isotope'

#describe 'making superatomdata from mspire isotope table' do
  #obj = Emass::SuperAtomData.create_from_isotope_data
  #obj.class.should == Emass::SuperAtomData
  #obj.size.should == 84
   #Argon
  #obj[:Ar].first.map(&:mass).should == [35.967545106, nil, 37.9627324, nil, 39.9623831225]
  #obj[:Ar].first.map(&:rel_area).should == [0.003365, 0, 0.000632, 0, 0.996003]
#end

describe 'can calculate isotope dists' do
  before do
    myhash = Mspire::Isotope::NIST::BY_ELEMENT.dup
    myhash[:C][0].relative_abundance = 0.988930
    myhash[:C][1].atomic_mass = 13.0033554
    myhash[:C][1].relative_abundance = 0.011070

    h = myhash[:H]
    h[0].atomic_mass = 1.0078246
    h[1].atomic_mass = 2.0141021
    h[0].relative_abundance = 0.99985 
    h[1].relative_abundance = 0.00015

    @sad = Emass::SuperAtomData.create_from_isotope_data(myhash)
  end
  
  it 'works' do
    # C2H4
    # formula: C2H4 charge : 0 limit: 0.000000e+00
    # 28.031298 100.000000
    # 29.034730 2.298792
    # 30.038298 0.013887
    # 31.044401 0.000008

    puts Emass::Pattern.calculate(Mspire::MF['C2H4'], @sad).to_s
  end

#  it 'works for Carbon' do
    #Emass::Pattern.calculate(Mspire::MF['C'], @sad).should == Emass::Pattern.from_doublets([12.0, 0.98893], [13.0033554, 0.01107])
  #end
end

#describe 'pruning a pattern' do
  #subject{ 
    #Emass::Pattern.new( [[11,0.01], [12,0.6], [13,0.2], [14,0.6], [15,0.01]].map {|data| Emass::Peak.new(*data) } ) 
  #}
  #it 'does not prune if limit is too low' do
    #reply = subject.prune!(0.001)
    #reply.should == subject
  #end

  #it 'prunes' do
    #reply = subject.prune!(0.01)
    #reply.should ==  Emass::Pattern.from_doublets([12,0.6], [13,0.2], [14,0.6])
  #end

  #it 'prunes ends while protecting low internal values' do
    #reply = subject.prune!(0.3)
    #reply.should ==  Emass::Pattern.from_doublets([12,0.6], [13,0.2], [14,0.6])
  #end
#end

#describe 'convoluting patterns' do
  ##alpha = Pattern.new( [[12,0.9],[13,0.1]].map {|data| Peak.new(*data) } )
  ##beta = Pattern.new( [[1,0.99],[2,0.1]].map {|data| Peak.new(*data) } )

  #alpha = Emass::Pattern.new( [[12,0.9],[13,0.1],[14,0.3]].map {|data| Emass::Peak.new(*data) } )
  #beta = Emass::Pattern.new( [[1,0.99],[2,0.1],[3,0.4]].map {|data| Emass::Peak.new(*data) } )

  #p alpha.convolute(beta)
#end
