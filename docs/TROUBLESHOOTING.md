# Troubleshooting

## Common Errors

**`skipped_no_samples` for a cell type**
Zero samples remained after metadata filtering, or the cell type name
doesn't match between your data and config. Check that:
1. Sample IDs in metadata match column names in count matrices (hyphens
   are automatically converted to dots)
2. `sample_filters` in config don't eliminate all samples
3. Cell type names match between `cell_mtx/` filenames and metadata

**`skipped_min_cells` after deduplication**
Fewer than 3 unique donors remain after dedup. This is common for rare
cell types in small cohorts. Either reduce filtering stringency or exclude
this cell type from analysis.

**`skipped_preflight` with "too few genes"**
Fewer than 10 genes passed the `min_gene_counts` / `min_gene_proportion`
filter. The cell type likely has very low expression depth. Try lowering
`min_gene_counts` or `min_gene_proportion` in the config.

**DESeq2 crashes with "every gene contains at least one zero"**
The count matrix is extremely sparse. This typically happens with very
rare cell types. The preflight check (minimum 3 samples, 10 genes) should
catch most cases, but borderline situations can still trigger this error.
Check the intermediate `counts_filtered.csv` to verify data quality.

**fGSEA returns 0 significant pathways**
Not necessarily an error. Check:
1. The ranking statistic produces a reasonable spread (inspect
   `results_deseq.tsv` stat column)
2. Gene symbols in the GMT file match your gene names (human symbols
   expected by default)
3. `fgsea.min_gene_set_size` is not too high for your gene set

**`no count matrix found` / `Could not find paired count matrices`**
The pipeline cannot locate interest + other matrices for this cell type.
Verify your directory layout matches one of the supported formats:
- `cell_mtx/{CellType}_persample_RNA_counts.tsv` + `allbut_mtx/...`
- `{CellType}_interest.counts.csv` + `{CellType}_other.counts.csv`

**Snakemake shows "Nothing to be done"**
All output files already exist. Delete `results/` (or the specific cell
type directory) to force a re-run.

**Sample name mismatches between metadata and count matrices**
Common when matrices use hyphens (`SRR-001`) but metadata uses dots
(`SRR.001`). Step 02 automatically converts hyphens to dots in both
metadata and matrix column names. If you still see mismatches, check for
other naming inconsistencies.

---

## Reading Log Files

Steps 02--04 write structured logs to `logs/{celltype}.log`:

```
[2026-04-03 19:43:38] [INFO ] celltype=Beta | step=start | ...
[2026-04-03 19:43:38] [INFO ] celltype=Beta | step=load_metadata | n_samples=191
[2026-04-03 19:43:38] [SKIP ] celltype=X | step=no_samples | reason=no_samples | ...
```

To see all warnings and errors after a run:

```bash
grep -rh "WARN\|ERROR\|SKIP" logs/
```

Step 05 also collects all flagged lines into
`results/summary/log_summary.tsv` for easy review.

---

## Resource Files

All resources are included in the repository -- no additional downloads
needed.

### Gene exclusion lists

Used by step 04 to remove ribosomal and mitochondrial genes before fGSEA
ranking:

```
resources/gsea_files/rpl_genes.csv     ribosomal proteins, large subunit
resources/gsea_files/rps_genes.csv     ribosomal proteins, small subunit
resources/gsea_files/mtr_genes.csv     mitochondrially encoded genes
```

### Pathway GMT files

Merged Reactome + KEGG pathways from MSigDB:

```
resources/gsea_files/reactome_kegg.gmt.txt     merged (default)
```

To use a different collection, set `fgsea.gene_sets_file` in your config.

**Updating pathways:** Download new GMT files from
[MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp) and point
`fgsea.gene_sets_file` to the new file.
