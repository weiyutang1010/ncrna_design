## Evaluation Data

Contains the designed solution for our work in ``Sampling-based Continuous Optimization for RNA Design''

- pyx_projection.csv: This Work (projection) optimizing for $p(y^\star \mid x)$ (Boltzmann probability)
- pyx_softmax.csv: This Work (softmax) optimizing for $p(y^\star \mid x)$ (Boltzmann probability)
- ned_softmax.csv: This Work (softmax) optimizing for NED$(x, y^\star)$ (normalized ensemble defect)
- dist_softmax.csv: This Work (softmax) optimizing for $d(\mathrm{MFE}(x), y^\star)$ (structural distance)
- ddg_softmax.csv: This Work (softmax) optimizing for $\Delta \Delta G^{\circ}(x, y^\star)$ (free energy gap)

### CSV Columns
- ID: Eterna puzzle ID
- Name: puzzle name
- Length: structure length
- Structure: secondary structure $y^\star$ (Eterna puzzle)
- Unpaired: number of unpaired nucleotides $\lvert \mathit{unpaired}(y) \rvert$
- Pairs: number of pairs nucleotides $\lvert \mathit{pairs}(y) \rvert$
- Total: number of valid sequences $\lvert \mathcal{X}(y^\star) \rvert$
- p(y | x): best $p(y^\star \mid x)$ found
- p(y|x): best $p(y^\star \mid x)$ sequence found
- NED: best NED$(x, y^\star)$ found
- NED seq: best NED$(x, y^\star)$ sequence found
- dist: best $d(\mathrm{MFE}(x), y^\star)$ found
- dist seq: best $d(\mathrm{MFE}(x), y^\star)$ sequence found
- DDG: best $\Delta \Delta G^{\circ}(x, y^\star)$ found
- DDG seq: best $\Delta \Delta G^{\circ}(x, y^\star)$ sequence found
- is_mfe: whether an mfe solution is found
- mfe_seq: one of the mfe solutions found
- is_umfe: whether an umfe solution is found
- umfe_seq: one of the umfe solution found
