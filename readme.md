# comparative-genomics-project

Main directory structure:
```text
.
├── analysis/   # -- analysis code, scripts/notebooks (after pipeline run)
├── data/
├── dev/       # -- development scripts, tools, installation scripts (reproducibility and logging)
├── docs/      # -- documentation files
├── logs/      # -- log files and version
├── readme.md  # -- project overview and instructions
├── output/    # -- output files and results
└── scripts/   # -- pipeline scripts
```

> [!IMPORTANT]  
> To mark a change you've added or done to the project, please commit with a clear message on what was changed and run  
> ``` ./dev/version_tracker.sh ```, follow the instructions (choose an option - default z unless you did a big change, write each step you've done in a seperate line and the press `ENTER` to get out of the prompt => version will be updated automatically)