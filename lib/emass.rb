require "emass/version"

module Emass

  # holds patterns
  class SuperAtomList < Array
  end

  # holds superatomlists
  class SuperAtomData < Array
  end

  Peak = Struct.new(:mass, :rel_area)

  # holds peaks
  class Pattern < Array

    # merge two patterns into one.  (convolute_basic in emass). Returns a new Pattern
    # reverse other, then convolve in a special way
    def convolute(other)
      o_sz = other.size
      other.reverse!
      pad = Array.new(o_sz-1)
      self.unshift(*pad) ; self.push(*pad)
      newpeaks = self.each_cons(o_sz).map do |mini_self|
        summass = 0.0
        sumweight = 0.0
        peaks = mini_self.zip(other) do |a,b|
          next unless a
          sumweight += (weight = a.rel_area * b.rel_area)
          summass += (a.mass+b.mass) * weight
        end
        Peak.new(sumweight == 0 ? nil : summass / sumweight, sumweight)
      end
      self.pop(o_sz-1) ; self.shift(o_sz-1)
      other.reverse!
      Pattern.new(newpeaks)
    end

    def prune!(limit)
      self.slice! 0...(self.index {|p| p.rel_area > limit })
      self.slice! self.rindex {|p| p.rel_area > limit }..-1
      self
    end

  end
end

include Emass

#alpha = Pattern.new( [[12,0.9],[13,0.1]].map {|data| Peak.new(*data) } )
#beta = Pattern.new( [[1,0.99],[2,0.1]].map {|data| Peak.new(*data) } )

alpha = Pattern.new( [[12,0.9],[13,0.1],[14,0.3]].map {|data| Peak.new(*data) } )
beta = Pattern.new( [[1,0.99],[2,0.1],[3,0.4]].map {|data| Peak.new(*data) } )

p alpha.convolute(beta)

abort 'here'

=begin
  // [g0,g1,g2]
  // [f0,f1,f2,f3]
  //
  // [g0,g1,g2][f0,f1,f2,f3]
  // [0  1  2   3  4  5  6
  h.clear();
  size_t g_n = g.size(); //3
  size_t f_n = f.size(); //4
  if(g_n == 0 || f_n == 0)
     return;
  // 0...(3+4-1 [i.e. 0..5])
  for(size_t k = 0; k < g_n + f_n - 1; k++) {  // k=0,1,2,3,4,5,6
    double sumweight = 0, summass = 0;
    size_t start = k < (f_n - 1) ? 0 : k - f_n + 1; // max(0, k-f_n+1)  // s=0,0,0,0,1,2,3
    size_t end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)  // e=0,1,2,2,2,2,2
    for(size_t i = start; i <= end; i++) {  // single it, i=0;0..1;0..2;0..2,1..2,2..2,3..2
      double weight = g[i].rel_area * f[k - i].rel_area; // g0*f0;g0*f1"+"g1*f0;g0*f2"+"g1*f1"+"g2*f0;g0*f3+g1*f2+g2*f1;g1*f3+g2*f2;g2*f3;
      double mass = g[i].mass + f[k - i].mass; // g0+f0;g0+f1"+"g1+f0;g0+f2"+"g1+f1"+"g2+g2+f0;
      sumweight += weight;  // ^
      summass += weight * mass; // w*m
    }
    peak p;
    if(sumweight == 0)
      p.mass = DUMMY_MASS;
    else
      p.mass = summass / sumweight;
    p.rel_area = sumweight;
    h.push_back(p);


// Merge two patterns to one.
void convolute_basic(Pattern & h, const Pattern & g, const Pattern & f)
{
  // [g0,g1,g2]
  // [f0,f1,f2]
  //
  // [g0,g1,g2][f0,f1,f2]
  // [0  1  2   3  4  5 ]
  h.clear();
  size_t g_n = g.size(); //3
  size_t f_n = f.size(); //3
  if(g_n == 0 || f_n == 0)
     return;
  // k=>0..4
  for(size_t k = 0; k < g_n + f_n - 1; k++) {  // k=0,1,2,3,4,5
    double sumweight = 0, summass = 0;                                  //       |
    size_t start = k < (f_n - 1) ? 0 : k - f_n + 1; // max(0, k-f_n+1)  // s=0,0,0,1,2,3 <- 0 for (0...f_n) then 0..g_n
    size_t end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)  // e=0,1,2,2,2,2
    for(size_t i = start; i <= end; i++) {  // single it, i=0;0..1;0..2;0..2,1..2,2..2,3..2
      double weight = g[i].rel_area * f[k - i].rel_area; // g0*f0;g0*f1"+"g1*f0;g0*f2"+"g1*f1"+"g2*f0;g1*f3+g1*f2+g2*f1;g1*f3+g2*f2;
      double mass = g[i].mass + f[k - i].mass;
      sumweight += weight;  // ^
      summass += weight * mass; // w*m
    }
    peak p;
    if(sumweight == 0)
      p.mass = DUMMY_MASS;
    else
      p.mass = summass / sumweight;
    p.rel_area = sumweight;
    h.push_back(p);
  }
}

=end
