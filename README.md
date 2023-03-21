### Incomplete-Market Models with Large Number of Heterogeneous Economic Agents



Contains Julia programs for computing macroeconomic equlibrium in certain incomplete-market models with large number of heterogeneous agents. The meaning and the purpose of these programs is exaplained in the paper &ldquo;<i>Dynamic transportation of economic agents</i>&rdquo;  &nbsp;[DTEA]&nbsp; by Andrew Lyasoff (last revised on Mar 21, 2023), which this repository supplements.

The code in &nbsp;&#8902;-RMT-&#8902;.jl&nbsp; implements the method described in Ch 18 of the book &ldquo;<i>Recursive Macroeconomic Theory</i>&rdquo;  &nbsp;[RMT]&nbsp; by L. Ljungvist and T. Sargent for calculating the equilibrium in Huggett's benchmark pure-exchange economy with no shared risk. The model parameters are borrowed from &nbsp;[RMT]&nbsp; and the output is illustrated.

The code in &nbsp;&#8902;-DIW-0&#8902;.jl&nbsp; implements the alternative algorithm described in Sec. 3 of &nbsp;[DIW]&nbsp; and yields a numerically verifiable equilibrium for the pure-echnage economy with no aggregate risk borrowed from &nbsp;[RMT]&nbsp; and covered by &nbsp;&#8902;-RMT-&#8902;.jl&nbsp;.

The code in &nbsp;&#8902;-DIW-1&#8902;.jl&nbsp; implements the algorithm described in Sec. 4 of &nbsp;[DIW]&nbsp; and yields an equilibrium for the benchmark economy with aggregate risk and production described in the paper &nbsp;&ldquo;<i>Income and wealth heterogeneity in the macroeconomy</i>&rdquo; &nbsp; by P. Krusell and A. Smith, Journal of Political Economy 106(1998) pp. 867-896. 

The code in &nbsp;functions-&#8902;.jl&nbsp; contains the actual programs. The code in &nbsp;ini-setup-&#8902;.jl&nbsp; defines the model parameters. The code in &nbsp;run-&#8902;.jl&nbsp; generates and stores the output in &nbsp;output-&#8902;.jls&nbsp;, which are included in this repository. The code in &nbsp;post-run-&#8902;.jl&nbsp; can read the respective output-&#8902;.jls and generates the plots included in &nbsp;[DIW]&nbsp;, in addition to other illustrations. 

All reported output is produced with generic Linux-on-x86 binaries (glibc) v1.8.4 (December 23, 2022) retreived from https://julialang.org/downloads/
