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

    puts Emass::Pattern.calculate(Mspire::MF['CH2'], @sad).to_s
    # 14.015649 100.000000
    # 15.019081 1.149396
    # 16.025302 0.000338
    
    #puts Emass::Pattern.calculate(Mspire::MF['C2H4'], @sad).to_s
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

describe 'convoluting patterns' do
  # takes an array of doublets
  def pattern_from_array(*array)
    Emass::Pattern.new( array.map {|data| Emass::Peak.new(*data) } )
  end

  def close_patterns(pattern1, pattern2)
    pattern2.each_with_index do |p2,i|
      p2.zip(pattern1[i]) do |lilp2, lilp1|
        expect(lilp1).to be_within(1e-10).of(lilp2)
      end
    end
  end

  it 'should give the same convolution as emass' do
    # these were taken directly from emass by printing the structures before
    # and after convolute
    alpha = pattern_from_array [0,1]
    beta = pattern_from_array [1.00782460,1.0], [2.01410210, 0.0001500225]
    alpha.convolute(beta).should == pattern_from_array([1.00782460, 1.0], [2.01410210, 0.0001500225])

    alpha = pattern_from_array [2.01564920, 1.0], [3.02192670, 0.0003000450], [4.02820420, 0.0000000225]
    beta = pattern_from_array [24.00000000, 1.0], [25.00335540, 0.0223878333], [26.00671080, 0.0001253038]
    close_patterns alpha.convolute(beta), [[26.01564920, 1.0],[27.01904324449037, 0.0226878783],[28.02250964963153, 0.0001320436], [29.02867613304374, 0.0000000381], [30.034915, 2.8193355000000002e-12]]
  end
end
