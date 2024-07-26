### Self-Aware Transport of Economic Agents


Contains Julia programs for computing macroeconomic equlibrium in models of Krusell-Smith or Aiyagari-Bewley-Huggett type. The meaning and the purpose of these programs is exaplained in the paper <a href="http://arxiv.org/abs/2303.12567">Self-Aware Transport of Economic Agents</a>;  &nbsp;[SATHA]&nbsp; by Andrew Lyasoff (last revised on July 22, 2024), which this repository supplements.

The code in &nbsp;&#8902;-RMT-&#8902;.jl&nbsp; implements the method described in Ch 18 of the book &ldquo;<i>Recursive Macroeconomic Theory</i>&rdquo;  &nbsp;[RMT]&nbsp; by L. Ljungvist and T. Sargent for calculating the equilibrium in Huggett's benchmark pure-exchange economy with no shared risk. The model parameters are borrowed from &nbsp;[RMT]&nbsp; and the output is illustrated.

The code in &nbsp;&#8902;-SATHA-0&#8902;.jl&nbsp; implements the alternative algorithm described in Sec. 3 of &nbsp;[SATHA]&nbsp; and yields a numerically verifiable equilibrium for the pure-echnage economy with no aggregate risk borrowed from &nbsp;[RMT]&nbsp; and covered by &nbsp;&#8902;-RMT-&#8902;.jl&nbsp;. The need for developing this alternative algorithm is explained in the paper [SATHA] (see Section 1).

The code in &nbsp;&#8902;-SATHA-&#8902;.jl&nbsp; implements the algorithm described in Sec. 4 of &nbsp;[SATHA]&nbsp; and yields an equilibrium for the benchmark economy with aggregate risk and production described in the paper &nbsp;&ldquo;<i>Income and wealth heterogeneity in the macroeconomy</i>&rdquo; &nbsp; by P. Krusell and A. Smith, Journal of Political Economy 106(1998) pp. 867-896. The reason for for developing this alternative algorithm is explained in the paper [SATHA] (see Section 1).

The code in &nbsp;functions-&#8902;.jl&nbsp; contains the actual programs. The code in &nbsp;ini-setup-&#8902;.jl&nbsp; defines the model parameters. The code in &nbsp;run-&#8902;.jl&nbsp; generates and stores the output in &nbsp;output-&#8902;.jls&nbsp;, which are included in this repository. The code in &nbsp;post-run-&#8902;.jl&nbsp; can read the respective output-&#8902;.jls and generates the plots included in &nbsp;[SATHA]&nbsp;, in addition to other illustrations. 

All reported output is produced with generic Linux-on-x86 binaries (glibc) v1.9.0 1.9.0 (2023-05-07 retreived from https://julialang.org/downloads/
