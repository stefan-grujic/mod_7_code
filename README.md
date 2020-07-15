#Script run order:

1. prokka_to_filestructure.py
2. Xstring_to_newick.R

>>> Run SCOARY using:

	scoary -g /dir/SCOARY/INPUT/cog_table.csv -t /dir/SCOARY/INPUT/cog_metadata.csv -n /dir/SCOARY/INPUT/tree_tmp.newick - e 50 -o dir/SCOARY/OUTPUT

3. sig_cog_table_gen_scoary.py
4. tSNE_scoary.R
5. logistic_regression_models.py
6. sensitivity_analysis.py
