"""move here:
- all strictly plotting-related functions

may want to break them into individual classes or modules later
"""

import matplotlib.pyplot as plt

from neutralb1 import utils

# use this in class, or in the big fit result class as an optional bool
WORKSPACE_DIR = utils.get_workspace_dir()
plt.style.use(f"{WORKSPACE_DIR}/config/neutralb1.mplstyle")
