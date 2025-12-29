__version__ = '1.3.3'

import multiprocessing as mp

from .index import index
from .fetch import fetch
from .profile import profile
from .rebuild import rebuild

try:
    mp.set_start_method('fork', force=True)
except RuntimeError:
    pass
