**Raw Summary Table** 

========================

Raw counts for number of query hits for taxa at taxonomic levels: species, genus, family, order, class, phylum, and superkingdom. :math:`^*` All taxids with query hits are included in the table.

- **Total** is the total number of queries that aligned to a taxid multiplied by the number of times the query appeared in the reads.
- **Unique** is the number of *unique* queries that aligned to a taxid.
- **Signature** is the number of queries that aligned *only* to this taxid multiplied by the number of times the queries appeared in the reads.
- **Unique_Signature** is the number of *unique* queries that aligned *only* to this taxid.
- **Weighted_Support** is the number of queries that aligned to a taxid, inversely weighted by the total number of taxids that the query aligned to. This is multiplied by the number of times the queries appeared in the reads (:math:`\sum 1/N_taxa \cdot N_copies`).
- **Unique_Weighted_Support** is the number of *unique* queries that aligned to a taxid, inversely weighted by the total number of taxids that query aligned to (:math:`\sum 1/N_taxa`).

------------

:math:`^*` Queries that have more than **max_taxa_per_query** hits are not included in the totals.
