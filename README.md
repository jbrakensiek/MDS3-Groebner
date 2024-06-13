# MDS3-Groebner
Groebner basis code for verifying MDS(3) constructions. This is supplementary to the paper

Improved field size bounds for higher order MDS codes
J Brakensiek, M Dhar, S Gopi - 2023 IEEE International Symposium on Information Theory
<https://arxiv.org/abs/2212.11262>

* `claim-5-2.jl` verifies Claim 5.2
* `mds3-groebner.jl` verifies Claim 5.4
* `mds3-groebner-char2.jl` gives an alternative construction of an MDS(3) over characteristic 2

The scrip `file.jl` can be run with `julia file.jl`. Note that the [Oscar.jl](https://github.com/oscar-system/Oscar.jl) library is required.

The intended output of `file.jl` is available at `file-out.txt`.
