# Pretrained localization model

Sandra Castillo, *Protein localization prediction model for CarveFungi* (2026).
Source: https://doi.org/10.5281/zenodo.18953301
License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Docker downloads the unmodified `model_binary_output10_18.keras` and verifies
SHA-256 `2765ce4f7cbf3a0dac9bd67e00c013df2a792bf820f2a06fb5a534d2d46f0ca0`.
`bin/predict_localization.py` adapts the accompanying notebook to batched CLI
inference, preserving its feature values, 750-residue limit, and score mapping.
