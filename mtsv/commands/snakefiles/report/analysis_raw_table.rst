**Raw Analysis Results Table**

===========================

The analysis table shows the results of the `equivalence hypothesis test <https://en.wikipedia.org/wiki/Equivalence_test>`_ 
(TOST: two-one-sided t-tests) for each candidate taxid (at the species and genus level).
Candidate taxids are determined using the user-provided **filter_params**. :math:`^*`

The equivalence test compares the difference (Δ) between the **expected** and **observed** ratio of **unique signature hits** to **total unique hits**.
The **expected** values are calculated by simulating and running queries (:math:`N =` **sample_n_kmers**) using sequences in the provided sequence database for each candidate taxid.
Approximately 10% of the queries for each species is used to update the expected value at the genus level (a simple bloom filter is used to avoid double counting the same query). 


The TOST procedure defines an upper (ΔU) and lower (–ΔL) equivalence bound that is based on the smallest effect size
of interest (based on the user-provided Cohen's **h** parameter). Two composite **null hypotheses** are tested: 
:math:`H_{01}`: Δ ≤ –ΔL and :math:`H_{02}`: Δ ≥ ΔU. When both of these one-sided tests can be statistically rejected,
we can conclude that –ΔL < Δ < ΔU.

| :math:`\Delta U = p_exp - \sin{\frac{h - (2 \times \sin^{-1}(\sqrt{p_exp})}{-2}}^2`
| :math:`\Delta L = p_exp - \sin{\frac{h + (2 \times \sin^{-1}(\sqrt{p_exp})}{2}}^2`
| Only taxids with at least one signature hit are included in the table.

- **Total_Hits** is the total number of queries that aligned to a taxid multiplied by the number of times the query appeared in the reads.
- **Unique_Hits** is the number of *unique* queries that aligned to a taxid.
- **Signature_Hits** is the number of queries that aligned *only* to this taxid multiplied by the number of times the queries appeared in the reads.
- **Unique_Signature_Hits** is the number of *unique* queries that aligned *only* to this taxid.
- **Weighted_Support** is the number of queries that aligned to a taxid, inversely weighted by the total number of taxids that the query aligned to. This is multiplied by the number of times the queries appeared in the reads (:math:`\sum 1/N_taxa \cdot N_copies`).
- **Unique_Weighted_Support** is the number of *unique* queries that aligned to a taxid, inversely weighted by the total number of taxids that query aligned to (:math:`\sum 1/N_taxa`).
- **Exp_Unique_Hits** is the number of *unique simulated* queries that aligned to the taxid.
- **Exp_Unique_Signature_Hits** is the number of *unique simulated* queries that aligned to *only* this taxid.
- **Exp_Prop(USH/UH)** is the ratio of Exp_Unique_Signature_Hits/Exp_Unique_Hits.
- **Obs_Prop(USH/UH)** is the observed Unique_Signature_Hits/Unique_Hits from the sample.
- **P_value** is the p-value calculated in the equivalence TOST.
- **Significance** indicates whether the p-value is statistically significant based on user specified **alpha**.
- **User_Filter** shows whether the taxid passed the user specified **filter_params**.

------------

:math:`^*` The **filter_params** ensure that any taxid that passes the filter will be included in the analysis;
however, if expected values are available (in the datastore from previous runs) for taxids that did not pass the filter, they will also be included in the analysis.
