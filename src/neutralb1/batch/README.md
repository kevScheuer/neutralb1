This directory handles any python scripts related to batch job submission. Right now some old cfg files and the run_fit.sh script still lives here
and will need to be reorganized, along with `submit.py`, due to the massive project restructuring that has occurred. In the future this should be part of 
a script at the project root, that simply uses something like a YAML file to pass the arguments to the submitter, allowing for some sort of history to know
what config was submitted to the farm. Could be simply included as part of a jupyter notebook.