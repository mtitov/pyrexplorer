 0.4.1 (2015-06-04)
  * performance enhancement

 0.4.0 (2015-05-29)
  * current implementation of SPADE will produce only maximal sequential patterns (no sub-sequences of frequent sequences)
  * name of master class is changed from SPADE to SPADEm (where "m" means "maximal")
  * removed option "--maximal" (it is applied by default)

 0.3.2 (2015-05-27)
  * fixed algorithm implementation (fixed incorrect logic and added corresponding methods, e.g. methods to deal with equivalence classes)
  * improved memory usage (for long sequences, with length more than 10^3, it is recommended to request maximal frequent sequences only)
  * option "--k" is changed to the option "--length" (maximum length of frequent sequences)
  * added option "--maximal" to get only maximal frequent sequences (no sub-sequences of frequent sequences)
  * added option "--sort" to sort output sequences (base on sequence length)

 0.2.1 (2015-04-27)
  * light enhancements

 0.2.0 (2015-04-20)
  * added option to get top-n longest sequential patterns (option "--top n")

 0.1.2
  * enhanced calculating of difference between sequences
  * improved procedure for processing of 1- and 2-sequences

 0.1.1 (2015-03-31)
  * doc-strings updates
  * moved sequence length checking procedure out of loop in method enumerate_frequent_sequences

 0.1.0 (2015-03-24)
  * first release
