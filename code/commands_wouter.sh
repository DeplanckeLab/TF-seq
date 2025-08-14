python code/convert_notebook.py code/ZZ_Endogenous/2-create_figures.py code/ZZ_Endogenous/README.md


python code/convert_notebook.py code/ZZ_Dose_comparison/dose_comparison.py code/ZZ_Dose_comparison/README.md

jupyter nbconvert code/ZZ_Power/power.ipynb   --to markdown   --output README   --output-dir code/ZZ_Power   --ExtractOutputPreprocessor.dir=".figures"

jupyter nbconvert code/ZZ_Scaling/scaling.ipynb   --to markdown   --output README   --output-dir code/ZZ_Scaling   --ExtractOutputPreprocessor.dir=".figures"