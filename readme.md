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
│   ├── family_sizes    # -- analysis of this part here
│   ├── duplicated_genes
│   └── TE_analysis
├── figures/
│   └── family_sizes    # -- figures from family sizes analysis
├── data                # -- large data_files (hidden)
├── docs                # -- documentation files (some notes if you wanna add something)
└── output              # -- output files and results of pipelines and codes
    ├── gene_lists      # -- here can put .txt files with lists of genes (small_families.txt, large_families.txt...)
    └── family_sizes    # -- any output files from family sizes analysis
```
<!-- 
> To mark a change you've added or done to the project, please commit with a clear message on what was changed and run  
> ``` ./dev/version_tracker.sh ```, follow the instructions (choose an option - default z unless you did a big change, write each step you've done in a seperate line and the press `ENTER` to get out of the prompt => version will be updated automatically)
 -->

## Accessing branches

If you already clones the repository locally, you can access any branch by running first the `fetch` command to get the latest changes from remote repository:
```bash
git fetch origin
```

Then you can checkout the branch you want to work on:
```bash
git checkout -b <branch_name> origin/<branch_name>
```
