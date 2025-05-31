## Evaluation Data

Designed sequences used for evaluation in ''Sampling-based Continuous Optimization for RNA Design''

- pyx_softmax.csv: SamplingDesign (softmax) optimizing for $p(y^\star \mid x)$ (Boltzmann probability)
- ned_softmax.csv: SamplingDesign (softmax) optimizing for $\mathrm{NED}(x, y^\star)$ (normalized ensemble defect)
- dist_softmax.csv: SamplingDesign (softmax) optimizing for $d(\mathrm{MFE}(x), y^\star)$ (structural distance)
- ddg_softmax.csv: SamplingDesign (softmax) optimizing for $\Delta \Delta G^{\circ}(x, y^\star)$ (free energy gap)
- pyx_projection.csv: SamplingDesign (projection) optimizing for $p(y^\star \mid x)$ (Boltzmann probability)

### CSV Columns
- ID: Eterna puzzle ID
- Name: puzzle name
- Length: puzzle length
- Structure: secondary structure $y^\star$ (Eterna puzzle)
- Unpaired: number of unpaired nucleotides $\lvert \mathit{unpaired}(y) \rvert$
- Pairs: number of base pairs $\lvert \mathit{pairs}(y) \rvert$
- Designability: "designable"/"undesignable"/"unknown"
- p(y|x): best $p(y^\star \mid x)$ value found
- p(y|x) seq: best $p(y^\star \mid x)$ sequence found
- NED: best $\mathrm{NED}(x, y^\star)$ value found
- NED seq: best $\mathrm{NED}(x, y^\star)$ sequence found
- dist: best $d(\mathrm{MFE}(x), y^\star)$ value found
- dist seq: best $d(\mathrm{MFE}(x), y^\star)$ sequence found
- DDG: best $\Delta \Delta G^{\circ}(x, y^\star)$ value found
- DDG seq: best $\Delta \Delta G^{\circ}(x, y^\star)$ sequence found
- is_mfe: whether an mfe solution is found
- mfe_seq: one of the mfe solutions found
- is_umfe: whether an umfe solution is found
- umfe_seq: one of the umfe solution found
