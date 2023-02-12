# cutflow1 = Cutflow.from_file("/path/to/cutflow1.cflow")
# cutflow2 = Cutflow.from_file("/path/to/cutflow2.cflow")
# ...

# Then, you can add up all of the yields in the cutflow by doing

# merged_cutflow = cutflow1 + cutflow2

# Or, you can make a "collection" of cutflows to display their yields side-by-side:

# CutflowCollection(
#     {"Group1": cutflow1, "Group2": cutflow2 + cutflow3, ...}
# )

#!/usr/bin/env python3

import os
import glob
import argparse
from tqdm import tqdm
from multiprocessing import Pool
from subprocess import Popen, PIPE 
from tools.cutflow import Cutflow, CutflowCollection

