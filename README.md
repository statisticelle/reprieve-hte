# reprieve-hte
This repository contains R code used for causal survival forest analysis of 5-year incident MACE and diabetes in the REPRIEVE trial. It is intended to provide more details on analysis decisions and model hyperparameters. 

Note that data is sensitive and is not made publicly available by this repository. De-identified data may be downloaded or requested elsewhere.

------

**Corresponding Paper:** Davies Smith et al. Characterizing the Heterogeneity of Statin Effects to Personalize Cardiovascular Disease Prevention in People with HIV. *Pending submission.* 

------

Causal survival forests are fit by `makedata/mk_crf_mace.R` and `makedata/mk_crf_diab.R`, and hyperparameter values can be found therein. We fit causal survival forests (CSFs) separately for MACE and diabetes. A 5-year time horizon was considered to balance the number of events and inverse probability of censoring weight (IPCW) stability. All risk factors were used to construct CSFs and estimate IPCWs. Ten-thousand honest trees were fit for each CSF, with 50% sample splitting for each tree and an additional 50% sample splitting for honest estimation. To ensure tree depth and sensitivity with few events, each tree node considered all risk factors, a minimum of 30 participants, and a minimum of 1% of observed events.

------

**Correspondence:** Emma Davies Smith (esmith at sdac.harvard.edu)
