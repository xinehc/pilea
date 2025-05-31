__version__ = '1.1.1'

import resource

from .index import index
from .fetch import fetch
from .profile import profile
from .rebuild import rebuild

## fix macOS limit
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
if soft_limit < 1024:
	resource.setrlimit(resource.RLIMIT_NOFILE, (1024, hard_limit))
