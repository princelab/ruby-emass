require 'spec_helper'

require 'emass'
require 'mspire/mf'
require 'mspire/isotope'


describe 'emass' do
  def close_patterns(pattern1, pattern2)
    pattern2.each_with_index do |p2,i|
      p2.zip(pattern1[i]) do |lilp2, lilp1|
        expect(lilp1).to be_within(1e-10).of(lilp2)
      end
    end
  end

  describe 'making superatomdata from mspire isotope table' do
    obj = Emass::SuperAtomData.create_from_isotope_data
    obj.class.should == Emass::SuperAtomData
    obj.size.should == 84
    obj[:Ar].first.map(&:mass).should == [35.967545106, nil, 37.9627324, nil, 39.9623831225]
    obj[:Ar].first.map(&:rel_area).should == [0.003365, 0, 0.000632, 0, 0.996003]
  end

  describe 'can calculate isotope dists' do
    before do
      # setting up an environment to exactly mimic emass so I can be sure we
      # are getting exactly the same answers
      myhash = Mspire::Isotope::NIST::BY_ELEMENT.dup
      c = myhash[:C]
      c[0].atomic_mass = 12.0  # <- should already be this
      c[0].relative_abundance = 0.988930
      c[1].atomic_mass = 13.0033554
      c[1].relative_abundance = 0.011070

      h = myhash[:H]
      h[0].atomic_mass = 1.0078246
      h[0].relative_abundance = 0.99985 
      h[1].atomic_mass = 2.0141021
      h[1].relative_abundance = 0.00015

      s = myhash[:S]
      s[0].atomic_mass = 31.972070
      s[0].relative_abundance = 0.9502
      s[1].atomic_mass = 32.971456  
      s[1].relative_abundance = 0.0075
      s[2].atomic_mass = 33.967866
      s[2].relative_abundance = 0.0421
      s[3].atomic_mass = 35.967080  
      s[3].relative_abundance = 0.0002

      @sad = Emass::SuperAtomData.create_from_isotope_data(myhash)
    end

    # these are the unnormalized results (as expected), but the normalized
    # result is exactly correct (verified with emass, hence these must be
    # exactly correct).

    it 'works with super atoms (e.g. H2)' do
      expected_isotope_pattern = [[2.0156492, 0.9997000225], 
                                  [3.0219267, 0.00029995499999999997], 
                                  [4.0282042, 2.2499999999999996e-08]]
      # this is the unnormalized result, but the normalized result is exactly
      # correct (verified with emass)
      Emass::Pattern.calculate(Mspire::MF['H2'], @sad).should == Emass::Pattern.from_doublets( *expected_isotope_pattern )
    end

    it 'works with more complex molecules (e.g., CH2)' do
      expected_isotope_pattern = [[14.0156492, 0.988633343250925], 
                                  [15.019080880184314, 0.011363313747224999],
                                  [16.025301550863496, 3.3427527749999995e-06],
                                  [17.0315596, 2.4907499999999994e-10]]
      Emass::Pattern.calculate(Mspire::MF['CH2'], @sad).should == Emass::Pattern.from_doublets( *expected_isotope_pattern )
    end

    it 'works for bigger molecules with gaps in isotopes' do
      actual_emass_output = %Q{
        779.597892 100.000000
        780.598983 27.280112
        781.593936 92.208971
        782.594910 23.779345
        783.590001 40.739147
        784.590876 9.945235
        785.586094 11.486218
        786.586885 2.657481
        787.582219 2.322238
        788.582942 0.509720
        789.578384 0.358657
        790.579050 0.074751
        791.574594 0.044019
        792.575216 0.008718
        793.570856 0.004410
        794.571445 0.000830
        795.567178 0.000368
        796.567742 0.000066
        797.563565 0.000026
        798.564114 0.000004
        799.560027 0.000002
      }.gsub(/^\s+/,'')

      Emass::Pattern.calculate(Mspire::MF['C10H20S20'], @sad).emass_normalize_truncate_percentize_stringify(6).should == actual_emass_output
    end

    it 'works on very large molecules in reasonable amount of time' do
      ksr1_aaseq = %Q{      
        MDRAALRAAAMGEKKEGGGGGDAAEGGAGAAASRALQQCGQLQKLIDISIGSLRGLRTKC
        AVSNDLTQQEIRTLEAKLVRYICKQRQCKLSVAPGERTPELNSYPRFSDWLYTFNVRPEV
        VQEIPRDLTLDALLEMNEAKVKETLRRCGASGDECGRLQYALTCLRKVTGLGGEHKEDSS
        WSSLDARRESGSGPSTDTLSAASLPWPPGSSQLGRAGNSAQGPRSISVSALPASDSPTPS
        FSEGLSDTCIPLHASGRLTPRALHSFITPPTTPQLRRHTKLKPPRTPPPPSRKVFQLLPS
        FPTLTRSKSHESQLGNRIDDVSSMRFDLSHGSPQMVRRDIGLSVTHRFSTKSWLSQVCHV
        CQKSMIFGVKCKHCRLKCHNKCTKEAPACRISFLPLTRLRRTESVPSDINNPVDRAAEPH
        FGTLPKALTKKEHPPAMNHLDSSSNPSSTTSSTPSSPAPFPTSSNPSSATTPPNPSPGQR
        DSRFNFPAAYFIHHRQQFIFPVPSAGHCWKCLLIAESLKENAFNISAFAHAAPLPEAADG
        TRLDDQPKADVLEAHEAEAEEPEAGKSEAEDDEDEVDDLPSSRRPWRGPISRKASQTSVY
        LQEWDIPFEQVELGEPIGQGRWGRVHRGRWHGEVAIRLLEMDGHNQDHLKLFKKEVMNYR
        QTRHENVVLFMGACMNPPHLAIITSFCKGRTLHSFVRDPKTSLDINKTRQIAQEIIKGMG
        YLHAKGIVHKDLKSKNVFYDNGKVVITDFGLFGISGVVREGRRENQLKLSHDWLCYLAPE
        IVREMTPGKDEDQLPFSKAADVYAFGTVWYELQARDWPLKNQAAEASIWQIGSGEGMKRV
        LTSVSLGKEVSEILSACWAFDLQERPSFSLLMDMLEKLPKLNRRLSHPGHFWKSAEINSS
        KVVPRFERFGLGVLESSNPKM
      }.each_line.map(&:strip).join
      ksr1_mf = Mspire::MF.from_aaseq(ksr1_aaseq)

      before = Time.now
      Emass::Pattern.calculate(ksr1_mf, @sad, 1e-30)
      puts "\t -- KSR1, with pruning (1e-30) took #{Time.now - before} sec"

      before = Time.now
      Emass::Pattern.calculate(ksr1_mf, @sad)
      puts "\t -- KSR1, with no pruning took #{Time.now - before} sec"

    end

  end

  describe 'pruning a pattern' do
    subject{ 
      Emass::Pattern.new( [[11,0.01], [12,0.6], [13,0.2], [14,0.6], [15,0.01]].map {|data| Emass::Peak.new(*data) } ) 
    }
    it 'does not prune if limit is too low' do
      reply = subject.prune!(0.001)
      reply.should == subject
    end

    it 'prunes' do
      reply = subject.prune!(0.01)
      reply.should ==  Emass::Pattern.from_doublets([12,0.6], [13,0.2], [14,0.6])
    end

    it 'prunes ends while protecting low internal values' do
      reply = subject.prune!(0.3)
      reply.should ==  Emass::Pattern.from_doublets([12,0.6], [13,0.2], [14,0.6])
    end
  end

  describe 'convoluting patterns' do

    it 'should give the same convolution as emass' do
      # these were taken directly from emass by printing the structures before
      # and after convolute
      alpha = Emass::Pattern.from_doublets [0,1]
      beta = Emass::Pattern.from_doublets [1.00782460,1.0], [2.01410210, 0.0001500225]
      alpha.convolute(beta).should == Emass::Pattern.from_doublets([1.00782460, 1.0], [2.01410210, 0.0001500225])

      alpha = Emass::Pattern.from_doublets [2.01564920, 1.0], [3.02192670, 0.0003000450], [4.02820420, 0.0000000225]
      beta = Emass::Pattern.from_doublets [24.00000000, 1.0], [25.00335540, 0.0223878333], [26.00671080, 0.0001253038]
      close_patterns alpha.convolute(beta), [[26.01564920, 1.0],[27.01904324449037, 0.0226878783],[28.02250964963153, 0.0001320436], [29.02867613304374, 0.0000000381], [30.034915, 2.8193355000000002e-12]]
    end
  end
end
