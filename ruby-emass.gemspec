# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'emass/version'

Gem::Specification.new do |spec|
  spec.name          = "ruby-emass"
  spec.version       = Emass::VERSION
  spec.authors       = ["John T. Prince"]
  spec.email         = ["jtprince@gmail.com"]
  spec.summary       = %q{ruby implementation of Haimi's Rockwood and Haimi's emass code.}
  spec.description   = %q{A pure ruby implementation of the Perttu Haimi's implementation of the Rockwood & Haimi 2006 JASMS "Efficient Calculation of Accurate Masses of Isotopic Peaks" method.}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  [
    #["nokogiri", "~> 1.6.1"],
    ["mspire-molecular_formula", "~> 0.1.0"],
  ].each do |args|
    spec.add_dependency(*args)
  end

  [
    ["bundler", "~> 1.6.2"],
    ["rake"],
    ["rspec", "~> 2.14.1"], 
    ["rdoc", "~> 4.1.1"], 
    ["simplecov", "~> 0.8.2"],
    ["fftw3"], # <- remove when can
  ].each do |args|
    spec.add_development_dependency(*args)
  end

end
