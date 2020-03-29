
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phyloR

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.0.1.9000-blue.svg)](https://github.com/cparsania/tidywrappers)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

``` r
library(magrittr)
library(dplyr)
```

The goal of `phyloR` is to pre-process the NCBI blast output for
downstream phylogenetic analysis.

## Assign column names to blast tabular output

``` r

blast_out_file <- system.file("extdata" ,"blast_output_01.txt", package = "phyloR")
blast_out_tbl <- readr::read_delim(blast_out_file , delim = "\t" , comment = "#" ,col_names = F)

colnames(blast_out_tbl) <- phyloR::get_blast_outformat_7_colnames()
colnames(blast_out_tbl)
#>  [1] "query_acc_ver"    "subject_acc_ver"  "identity"         "alignment_length"
#>  [5] "mismatches"       "gap_opens"        "q_start"          "q_end"           
#>  [9] "s_start"          "s_end"            "evalue"           "bit_score"       
#> [13] "positives"
```

## Filter blast hits

``` r

## evalue <= 1e-5

blast_out_tbl %>% phyloR::filter_blast_hits(evalue = 1e-5)
#> # A tibble: 545 x 13
#>    query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>    <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#>  1 KAE8371401.1  KAE8371401.1       100                253          0         0
#>  2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#>  3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#>  4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#>  5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#>  6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#>  7 KAE8371401.1  GES63448.1          39.3              229        135         3
#>  8 KAE8371401.1  XP_001214927.1      39.3              229        135         3
#>  9 KAE8371401.1  OQD86157.1          38.6              228        132         3
#> 10 KAE8371401.1  XP_002375911.1      39.1              230        129         5
#> # … with 535 more rows, and 7 more variables: q_start <dbl>, q_end <dbl>,
#> #   s_start <dbl>, s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>

##  Filter blast hits : bitscore >= 200

blast_out_tbl %>% phyloR::filter_blast_hits(bit_score = 200)
#> # A tibble: 6 x 13
#>   query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>   <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#> 1 KAE8371401.1  KAE8371401.1       100                253          0         0
#> 2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#> 3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#> 4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#> 5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#> 6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#> # … with 7 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>,
#> #   s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>

##  Filter blast hits : query coverage >= 90

blast_out_tbl %>% phyloR::filter_blast_hits(query_cov = 90, query_length = 253)
#> # A tibble: 89 x 13
#>    query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>    <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#>  1 KAE8371401.1  KAE8371401.1       100                253          0         0
#>  2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#>  3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#>  4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#>  5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#>  6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#>  7 KAE8371401.1  GES63448.1          39.3              229        135         3
#>  8 KAE8371401.1  XP_001214927.1      39.3              229        135         3
#>  9 KAE8371401.1  XP_002375911.1      39.1              230        129         5
#> 10 KAE8371401.1  OQD72420.1          40                230        127         5
#> # … with 79 more rows, and 7 more variables: q_start <dbl>, q_end <dbl>,
#> #   s_start <dbl>, s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>

##  Filter blast hits : percent identity >= 90

blast_out_tbl %>% phyloR::filter_blast_hits(identity = 90)
#> # A tibble: 1 x 13
#>   query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>   <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#> 1 KAE8371401.1  KAE8371401.1         100              253          0         0
#> # … with 7 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>,
#> #   s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>

##  Filter blast hits : evalue <= 1e-5, bitscore >= 200, query coverage >= 90, percent identity >= 60, 

blast_out_tbl %>% phyloR::filter_blast_hits(evalue = 1e-5, bit_score = 200, query_cov = 90, query_length = 253, identity = 60)
#> # A tibble: 6 x 13
#>   query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>   <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#> 1 KAE8371401.1  KAE8371401.1       100                253          0         0
#> 2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#> 3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#> 4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#> 5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#> 6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#> # … with 7 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>,
#> #   s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>
```

## Format fasta headers

``` r

fa_file <- system.file("extdata" ,"blast_output_01.fasta", package = "phyloR")
fa <- Biostrings::readBStringSet(fa_file)

## existing headers 
names(fa) %>% head()
#> [1] "KAE8371401.1:1-253 trypsin-like cysteine/serine peptidase domain-containing protein [Aspergillus bertholletius]" 
#> [2] "XP_031924425.1:1-247 trypsin-like cysteine/serine peptidase domain-containing protein [Aspergillus caelatus]"    
#> [3] "XP_031940570.1:1-247 trypsin-like cysteine/serine peptidase domain-containing protein [Aspergillus pseudonomius]"
#> [4] "KAB8261773.1:1-247 trypsin-like cysteine/serine peptidase domain-containing protein [Aspergillus pseudonomius]"  
#> [5] "XP_022392557.1:1-228 hypothetical protein ABOM_003016 [Aspergillus bombycis]"                                    
#> [6] "OQD81666.1:1-248 hypothetical protein PENANT_c026G03216 [Penicillium antarcticum]"

## new headers 
fa_file %>% phyloR::format_fasta_headers() %>% head()
#>   A BStringSet instance of length 6
#>     width seq                                               names               
#> [1]   253 MKSLQVLGNLLGVAFLMTQPVNA...NGRPDVNNDPVVVTSWINTLRAD KAE8371401_1__1_2...
#> [2]   247 MKSVKALGNLPYCILLIAQSVSA...YCEPHSNGRPDVNNDPLSASDWI XP_031924425_1__1...
#> [3]   247 MKPMRSVGRLLFFISFIALPVNA...YCEPSSNGRPDVNNDALSANSWI XP_031940570_1__1...
#> [4]   247 MKPMRSVGRLLFFISFIALPVNA...YCEPSSNGRPDVNNDALSANSWI KAB8261773_1__1_2...
#> [5]   228 MKPTKALGRLLFFTSIIASPVNA...YCQPSSNGRPDVNNDPLSADDWI XP_022392557_1__1...
#> [6]   248 MKPTQGFSNLLCLASLIAQPVNA...CDKASNGRPDINSDTLPAKDWID OQD81666_1__1_248...
```

## Remove redundancy from blast tabular ouput

``` r
blast_out_tbl %>% phyloR::remove_redundant_hits()
#> # A tibble: 716 x 13
#> # Groups:   subject_acc_ver [716]
#>    query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>    <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#>  1 KAE8371401.1  1FN8_A              31.6              209        127         6
#>  2 KAE8371401.1  1GDU_A              31.9              210        127         6
#>  3 KAE8371401.1  1PPZ_A              31.9              210        127         6
#>  4 KAE8371401.1  AAQ04074.1          27.5              218        135         9
#>  5 KAE8371401.1  AAR91718.1          30.8              117         73         2
#>  6 KAE8371401.1  AAR91719.1          28.7              157         94         6
#>  7 KAE8371401.1  AAW31593.1          27.2              228        147         6
#>  8 KAE8371401.1  ADK37838.1          29.5              227        136         9
#>  9 KAE8371401.1  ADY16698.1          35.5              242        151         4
#> 10 KAE8371401.1  AER28315.1          30.5              236        148         5
#> # … with 706 more rows, and 7 more variables: q_start <dbl>, q_end <dbl>,
#> #   s_start <dbl>, s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>
```

## Add taxonomy columns to blast tabular output

``` r

## add kingdom

with_kingdom <- blast_out_tbl %>% slice(1:10) %>% phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver", 
                                                                       taxonomy_level = "kingdom")
#> ✓ done.  Time taken 5.2
#> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ● Rank search begins...
#> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ✓ done.  Time taken 0.190485954284668

with_kingdom
#> # A tibble: 10 x 15
#>    query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>    <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#>  1 KAE8371401.1  KAE8371401.1       100                253          0         0
#>  2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#>  3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#>  4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#>  5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#>  6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#>  7 KAE8371401.1  GES63448.1          39.3              229        135         3
#>  8 KAE8371401.1  XP_001214927.1      39.3              229        135         3
#>  9 KAE8371401.1  OQD86157.1          38.6              228        132         3
#> 10 KAE8371401.1  XP_002375911.1      39.1              230        129         5
#> # … with 9 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>,
#> #   s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>, taxid <chr>,
#> #   kingdom <chr>

## add species 

with_kingdom_and_species <- with_kingdom %>% phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver" ,
                                                                  taxonomy_level = "species")
#> ── WARNING ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ── WARNING ENDS ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ● Rank search begins...
#> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ✓ done.  Time taken 0.0854787826538086

with_kingdom_and_species
#> # A tibble: 10 x 16
#>    query_acc_ver subject_acc_ver identity alignment_length mismatches gap_opens
#>    <chr>         <chr>              <dbl>            <dbl>      <dbl>     <dbl>
#>  1 KAE8371401.1  KAE8371401.1       100                253          0         0
#>  2 KAE8371401.1  XP_031924425.1      71.3              247         71         0
#>  3 KAE8371401.1  XP_031940570.1      66.4              247         83         0
#>  4 KAE8371401.1  KAB8261773.1        66.0              247         84         0
#>  5 KAE8371401.1  XP_022392557.1      64.4              247         69         1
#>  6 KAE8371401.1  OQD81666.1          61.7              248         95         0
#>  7 KAE8371401.1  GES63448.1          39.3              229        135         3
#>  8 KAE8371401.1  XP_001214927.1      39.3              229        135         3
#>  9 KAE8371401.1  OQD86157.1          38.6              228        132         3
#> 10 KAE8371401.1  XP_002375911.1      39.1              230        129         5
#> # … with 10 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>,
#> #   s_end <dbl>, evalue <dbl>, bit_score <dbl>, positives <dbl>, taxid <chr>,
#> #   kingdom <chr>, species <chr>
```
