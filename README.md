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

The container performs **reconstruction only**, from precomputed:

- Legacy eggNOG-mapper annotations (four header lines; query ID in column 1,
  score in column 4, EC in column 8).
- CarveFungi localization CSV (`.loc_pred`) with `Ids`, `E`, `M`, `P`, `O` columns.
  Gene IDs must match the annotation file.

It includes the universal reaction database, but not eggNOG-mapper, PSIPRED,
TensorFlow, or pretrained localization models. Current eggNOG-mapper output
needs conversion before use. The bundled yeast annotation/localization files
use different gene identifiers and cannot directly serve as a paired example.

### Docker

From this directory, build for a typical x86-64 HPC node (use `linux/arm64`
instead on ARM; emulated builds/runs on a different architecture can be slower):

```bash
docker build --platform linux/amd64 -t carvefungi:scip .
docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" carvefungi:scip \
  --id fungus --eggnog /work/input.emapper.annotations \
  --localization /work/input.loc_pred --output /work/output
```

`--output` must be new or empty. The CLI uses a temporary workspace to support
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
  python /opt/carvefungi/bin/carvefungi.py \
  --id fungus --eggnog /work/input.emapper.annotations \
  --localization /work/input.loc_pred --output /work/output
```

Use the same `exec` command in your scheduler job script. `singularity` can
replace `apptainer`. No GPU is needed for reconstruction. Solver search uses
SCIP's defaults; OpenMP/BLAS threads default to one in the image. Runtime and
memory requirements for real proteomes have not been benchmarked.

### Local installation and tests

```bash
python3.11 -m venv .venv
. .venv/bin/activate
pip install -r requirements-reconstruction.txt
python -m unittest discover -s tests -v
python bin/carvefungi.py --help
```

Tests cover synthetic-network constraints, scores, infeasibility, solution
filtering, COBRApy SBML round-tripping, workspace isolation, and construction
of the bundled universal model's MILP. They run during Docker builds. They do
not validate end-to-end biological reconstruction or equivalence to CPLEX.

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

