# cl-buhlmann
Common Lisp implementation of the Buhlmann decompression model just for fun.
I used [this document written by Paul Chapman](https://aquatec.wordpress.com/wp-content/uploads/2011/03/decompression-theory.pdf) to understand the theory and write most if not all of the code.
As such there is probably stuff I don't fully understand and it goes without saying that this shouldn't be used to actually plan any dives.

## TODO
[Wikipedia](https://en.wikipedia.org/wiki/B%C3%BChlmann_decompression_algorithm#Tissue_inert_gas_limits) says that the Buhlmann model also specifies how the constants for multiple inert gases when both nitrogen and helium are present in the tissue, using:
$$ a = a_{\text{N}_2} (1-R) + a_{\text{He}} R, $$
and
$$ b = b_{\text{N}_2} (1-R) + b_{\text{He}} R. $$
