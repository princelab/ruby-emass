require "emass/version"
require 'mspire/molecular_formula'

# isotope file format:
# Pb  4
# 203.973020  0.014
# 205.974440  0.241

#def init_data(isotope_file)
  #sad = Emass::SuperAtomData.new
  #new_super_atom = nil
  #IO.foreach(isotope_file) do |line|
    #line.chomp!
    #if line[/\A[[:alpha:]]/]
      
    #if line.size > 0

      ## SuperAtomList
    #else
      #sad << new_super_atom

    #if 
  #end
#end

module Emass

  # ElemMap // map from element abbreviation to index in the elements table
  # FormMap // map from element index to the count of occurences in the formula

  # holds superatomlists. Indexed by element (element_number in original)
  class SuperAtomData < Hash
    # returns a new array with a nil for each missing mass_number within the
    # array of isotopes
    def self.insert_nils(elemental_isotopes)
      mass_num_cnt = elemental_isotopes.first.mass_number
      is_tmp = [elemental_isotopes.first]
      elemental_isotopes[1..-1].each do |isotope|
        loop do
          if isotope.mass_number - mass_num_cnt > 1
            is_tmp << nil
            mass_num_cnt += 1
          else
            mass_num_cnt += 1
            break
          end
        end
        is_tmp << isotope
      end
      is_tmp
    end

    # returns a superatomlist
    def self.create_from_isotope_data(isotopes_by_element=Mspire::Isotope::BY_ELEMENT)
      sad = self.new
      isotopes_by_element.each do |el, elemental_isotopes|
        with_nils = insert_nils(elemental_isotopes)
        # fill in missing isotopes with a nil
        observe = (elemental_isotopes.map(&:mass_number).first == 32)
        pattern = Pattern.new( 
          with_nils.map do |iso|
            (mass, rel_abund) = iso ? [iso.atomic_mass, iso.relative_abundance] : [nil,0]
            Peak.new(mass, rel_abund)
          end
        )
        sad[el] = SuperAtomList.new([pattern])
      end
      sad
    end
  end

  # holds patterns.  Indexed by bit_number (??) may hold up to 8 superatoms
  class SuperAtomList < Array
  end

  Peak = Struct.new(:mass, :rel_area)

  # holds peaks. indexed by peak number
  class Pattern < Array
    DIGITS_TO_PRINT = 6
    LIMIT = 0.0   # 1e-30 usually good enough, too

    # takes an array of doublets, casts them into Emass::Peaks
    # and returns a Pattern object
    def self.from_doublets(*array)
      self.new( array.map {|ar| Peak.new(*ar) } )
    end

    def self.calculate(formula, super_atom_data, limit=LIMIT)
      puts "STARTING CALCULATE"
      result = Pattern[Peak.new(0.0, 1.0)] # <- 0 mass, 1.0 area

      # for(FormMap::iterator i = fm.begin(); i != fm.end(); i++) {
      p formula
      formula.delete(:C)
      formula[:C] = 1
      p formula

      formula.each do |el,cnt|
        sal = super_atom_data[el]
        n = cnt
        j = 0
        while n > 0
          sz = sal.size
          puts "sz: #{sz}"
          puts "n: #{n}"
          if j == sz
            puts "j == sz"
            #sal << Pattern.new
            puts "SAL j-1:"
            puts sal[j-1].to_s
            sal[j] = sal[j-1].convolute(sal[j-1])
            puts "SAL j:"
            puts sal[j].to_s
            sal[j].prune!(limit)
          end
          if (n & 1) != 0
            puts "n & 1"
            result = result.convolute(sal[j])
            result.prune!(limit)
          end
          n = (n >> 1)
          j += 1
        end
      end
      puts "ENDING CALCULATE"
      result

      ## take charge into account, if any

    end

    # merge two patterns into one.  (convolute_basic in emass). Returns a new Pattern
    # reverse other, then convolve in a special way
    def convolute(other)
      #puts "SELF:"
      #p self
      #puts "SELF:"
      #p other
      o_sz = other.size
      other.reverse!
      pad = Array.new(o_sz-1)
      self.unshift(*pad) ; self.push(*pad)
      newpeaks = self.each_cons(o_sz).map do |mini_self|
        summass = 0.0
        sumweight = 0.0
        peaks = mini_self.zip(other) do |a,b|
          next unless a && b
          sumweight += (weight = a.rel_area * b.rel_area)
          summass += (a.mass+b.mass) * weight
        end
        Peak.new( 
                 sumweight == 0 ? nil : (summass / sumweight), 
                 sumweight
                )
      end
      self.pop(o_sz-1) ; self.shift(o_sz-1)
      other.reverse!
      Pattern.new(newpeaks)
    end

    def prune!(limit)
      self.slice! 0...(self.index {|p| p.rel_area > limit })
      self.slice! (self.rindex {|p| p.rel_area > limit }+1)..-1
      self
    end

    # similar to print_pattern
    def to_s(digits=DIGITS_TO_PRINT)
      max_area = self.map(&:rel_area).max
      return '' if max_area == 0

      # wcout.setf(ios::fixed);  <- format fixed (scientific??)
      # wcout.precision(digits);  <- num digits to round
      print_limit = (10.0**-digits) / 2

      self.map do |peak|
        mass = peak.mass
        rel_area = peak.rel_area
        val_perc = rel_area / max_area * 100
        if mass && val_perc >= print_limit
          sprintf "%.#{digits}f %.#{digits}f\n", mass, val_perc
        end
      end.join
    end

  end
end
