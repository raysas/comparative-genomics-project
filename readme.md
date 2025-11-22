# comparative-genomics-project

[todo list for this project](TODO.md)

You can find the duplicated genes (gene families output, from pipeline_1) in [`output/clusters/`](output/clusters/) in 2 formats:
- `.txt` : mcl output format (space-separated gene names per line, line=family)
- `.tsv` : tab-separated format (geneName \t familyID, mapping of each gene to its family)

*preferable to use tsv format, easier and clearer to handle*

## Important notes

> [!CAUTION]
> never at any cost merge two branches together (at least discuss first), i.e. the command `git merge...` should never be used in this project. Always work in one branch, if something is needed from another branch, first discuss it, will be made available in a proper way (risk of deleting files or messing up the code otherwise)

The main directories to be aware of in master branch
```text
.
├── TODO.md             # -- general todo list
├── analysis
│   ├── Ks
│   ├── duplicated_genes
│   └── TE_analysis
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

At the end all of our works will be merged, this is why preserving these main directory names is essential to avoid conflicts.
Each specific branch will have a specific readme.md with an `examples` folder containing example data files to test the code and pipelines (`pipeline_2` and `general_analysis` branches already have it).

These files are important, change them in your respective branches if needed, will merge later on:

* requirements.txt : contains the list of python packages needed for the project (if any new package is needed, please add it to the list)
* scripts/SETUP.sh : contains installation commands for softwares needed for the project (if any new software is needed, please add it to the script)
* TODO.md : contains the todo list

> [!NOTE]
> Each branch has a different todo list tailored to its specific tasks, in `master` branch you will find the general todo list for the whole project


## Branches in this repository

* TE
* dev (try moving in and out there before doing stuff on master)
* general_analysis
* master
* pipeline_1
* pipeline_2
* ppi

## Accessing branches

If you already clones the repository locally, you can access any branch by running first the `fetch` command to get the latest changes from remote repository:
```bash
git fetch origin
```

Then you can checkout the branch you want to work on:
```bash
git checkout -b <branch_name> origin/<branch_name>
```

### for pipeline 2

```bash
# -- first time accessing the branch
git fetch origin
git checkout -b pipeline_2 origin/pipeline_2
# -- making sure you're on the right branch
git branch # this should highlight pipeline_2

# -- then every time you wanna push your changes on github:
git add .
git commit -m "your message here"
git push origin pipeline_2
```

### for general_analysis

```bash
# -- first time accessing the branch
git fetch origin
git checkout -b general_analysis origin/general_analysis
# -- making sure you're on the right branch
git branch # this should highlight general_analysis

# -- then every time you wanna push your changes on github:
git add .
git commit -m "your message here"
git push origin general_analysis
```