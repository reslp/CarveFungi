# CarveFungi

CarveFungi is a genome-scale metabolic model reconstruction pipeline able to create a compartmentalized metabolic model of any fungal specie from its protein sequences. It is implemented using python.

## Open-source reconstruction (SCIP)

This fork replaces CPLEX with **SCIP/PySCIPOpt**. The original MILP formulation,
reaction scoring, and downstream ensemble workflow are retained. COBRApy uses
GLPK for FBA; no commercial solver or license is needed.

SCIP's stored feasible solutions are filtered to within 0.1% of the best
incumbent score, then the existing code selects up to five per reconstruction.
This is **not** CPLEX's active solution-pool search: fewer alternatives (possibly
one) and different ensembles are expected. Non-optimal incumbents produce a
warning with the MIP gap; no feasible solution raises an error. The existing
`feast` and `opti` arguments now set SCIP's feasibility and dual-feasibility
tolerances, respectively. Biological scoring assumptions are unchanged.

### Inputs and scope

The container supports **localization prediction from protein FASTA** and
**reconstruction** from:

- Legacy eggNOG-mapper annotations (four header lines; query ID in column 1,
  score in column 4, EC in column 8).
- CarveFungi localization CSV (`.loc_pred`) with `Ids`, `E`, `M`, `P`, `O` columns.
  Gene IDs must match the annotation file.

It includes the universal reaction database, TensorFlow/Keras, and Sandra
Castillo's pretrained [localization model](https://doi.org/10.5281/zenodo.18953301)
(CC BY 4.0; downloaded and checksum-verified during the build). No model download
or network access is needed at runtime. eggNOG-mapper is not included; PSIPRED
is not needed for this newer predictor. Current eggNOG-mapper output needs
conversion before reconstruction. The bundled yeast annotation/localization
files use different gene identifiers and are not a paired example.

### Docker

From this directory, build for a typical x86-64 HPC node (use `linux/arm64`
for native ARM). TensorFlow on x86 requires AVX. When cross-building x86 on an
ARM Mac without AVX emulation, add `--build-arg SKIP_MODEL_TEST=1` to the build
command; all other tests still run. Model inference then needs an AVX-capable
x86 host and is not validated by that cross-build.

```bash
docker build --platform linux/amd64 -t carvefungi:scip .
# Predict from the same protein FASTA used for eggNOG-mapper:
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" --entrypoint python carvefungi:scip \
  /opt/carvefungi/bin/predict_localization.py \
  --fasta /work/proteins.faa --output /work/input.loc_pred --batch-size 8
# Reconstruct:
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" carvefungi:scip \
  --id fungus --eggnog /work/input.emapper.annotations \
  --localization /work/input.loc_pred --output /work/output
```

The predictor writes a new CSV file (never overwrites), retaining FASTA IDs.
It preserves the author's preprocessing: first 750 residues, zero padding,
unknown token = 0, and `O = cytoplasm + Golgi` (not normalized). Its output
columns are `E,M,P,O,Ex,N,Ids`; reconstruction ignores `Ex` and `N`. Inference
uses bounded batches; decrease `--batch-size` if memory is limited. `--model`
can override the bundled model path. A failed prediction may leave a partial
CSV, which must not be used for reconstruction.

For reconstruction, `--output` must be new or empty. The CLI uses a temporary workspace to support
read-only container installations; outputs are written directly to the mounted
directory. The original resume logic is not used. Original SBML ensemble and
intermediate annotation outputs are retained.

### Apptainer / Singularity

Save the Docker image and transfer the archive to the cluster, then build a SIF
on a node with Apptainer (no registry required):

```bash
docker save -o carvefungi.tar carvefungi:scip
# On the cluster:
apptainer build carvefungi.sif docker-archive://carvefungi.tar
apptainer exec --cleanenv --bind "$PWD:/work" carvefungi.sif \
  python /opt/carvefungi/bin/predict_localization.py \
  --fasta /work/proteins.faa --output /work/input.loc_pred
apptainer exec --cleanenv --bind "$PWD:/work" carvefungi.sif \
  python /opt/carvefungi/bin/carvefungi.py \
  --id fungus --eggnog /work/input.emapper.annotations \
  --localization /work/input.loc_pred --output /work/output
```

Use the same `exec` command in your scheduler job script. `singularity` can
replace `apptainer`. Both stages run on CPU; no GPU is needed. Solver search uses
SCIP's defaults; OpenMP/BLAS threads default to one in the image. Runtime and
memory requirements for real proteomes have not been benchmarked.

### Local installation and tests

```bash
python3.11 -m venv .venv
. .venv/bin/activate
pip install -r requirements-reconstruction.txt
# Also needed for localization prediction:
pip install -r requirements-reconstruction.txt -r requirements-prediction.txt
python -m unittest discover -s tests -v
python bin/carvefungi.py --help
```

Tests cover synthetic-network constraints, scores, infeasibility, solution
filtering, COBRApy SBML round-tripping, workspace isolation, and construction
of the bundled universal model's MILP. Prediction tests check preprocessing,
batching, score mapping, and inference against three published example predictions.
Tests run during Docker builds; model inference can be explicitly skipped for
cross-builds and is skipped locally if the model has not been downloaded. They do not validate biological accuracy, end-to-end
reconstruction, or equivalence to CPLEX.

## Upstream protein-to-model workflow

Requirements for the upstream annotation/localization stages:
- Tensorflow
- keras
- biopython
- pandas
- numpy
- cobrapy
- EggNog annotation software: https://github.com/eggnogdb/eggnog-mapper
- Download the deep learning models from " https://doi.org/10.5281/zenodo.18953301".


 CarveFungi uses a deep learning model to predict the cellular localization of the proteins and this information is used to score the reactions. 
![Deep neural network](/images/CNN.png)

CarveMe (https://doi.org/10.1093/nar/gky537,https://github.com/cdanielmachado/carveme) uses the scoring of the reactions to obtain a functional metabolic model that is able to produce biomass. 

