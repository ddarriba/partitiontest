# Logging and Checkpointing #

Partitiontest generates the log and checkpoint files in the desired output folder. By default, this folder is partest\_ALIGN\_FILE, but it can be set trough the configuration file.


# Details #

Output folder contains:
  * _models_ file, with the result of the model selection for each partition.
  * _schemes_ file, with the result of the schemes selection on each step of the algorithm.
  * _results_ file, with a XML summary of the best-fitting scheme.
  * _ckpfiles_ directory, with the checkpointing files.

# Checkpointing #

Each partitioning scheme identifier is coded in base64, and each partition is located through its unique ID. This ID is used as the checkpoint filename. However, length for file names is restricted on each Operating System, thus for checkpointing this ID is divided into "filename", shorter than this limit, and a "hash" code, for locating each single partition inside the checkpointing file.