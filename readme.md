# comparative-genomics-project

[todo list for this project](TODO.md)

## Important notes

> [!CAUTION]
> never at any cost merge two branches together (at least discuss first), i.e. the command `git merge...` should never be used in this project. Always work in one branch, if something is needed from another branch, first discuss it, will be made available in a proper way (risk of deleting files or messing up the code otherwise)

The main directories to be aware of in master branch
```text
.
├── TODO.md             # -- general todo list
├── analysis
│   ├── Ks
│   ├── code.Rproj
│   └── duplicated_genes
├── data                # -- large data_files (hidden)
├── docs                # -- documentation files (some notes if you wanna add something)
├── output              # -- output files and results of pipelines and codes
├── pipeline            # -- overall pipeline structure (merge of all sub-pipelines)
├── pipeline_1          # -- pipeline from raw to duplicated genes
├── pipeline_2          # -- next pipeline stage: from dup genes to Ks computation
├── requirements.txt    # -- python packages installed for this project
└── scripts             # -- scripts used for seperate tasks in the project
```
<!-- 
> To mark a change you've added or done to the project, please commit with a clear message on what was changed and run  
> ``` ./dev/version_tracker.sh ```, follow the instructions (choose an option - default z unless you did a big change, write each step you've done in a seperate line and the press `ENTER` to get out of the prompt => version will be updated automatically)
 -->

At the end all of our works will be merged in the `main` branch, this is why preserving these main directory names is essential to avoid conflicts.

Each specific branch will have a specific readme.md with an `examples` folder containing example data files to test the code and pipelines.

## Transposable Elements (TEs)

For this branch there are 2 types of codes:

Analysis, annotation and that kind of stuff regarding TEs: `analysis/TE_analysis` (or rename/change)

Scripts related to TEs and covergae computation (heavy computation/metrics type of codes on terminal): add to `scripts/`)