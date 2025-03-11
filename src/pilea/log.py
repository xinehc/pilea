import os
import sys
import logging

from tqdm import tqdm
from functools import partialmethod

logging.basicConfig(
    level='INFO',
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

if sys.stderr.isatty():
    logging.addLevelName(logging.INFO, f'\033[1m{logging.getLevelName(logging.INFO)}\033[1;0m')
    logging.addLevelName(logging.WARNING, f'\033[1m\x1b[33;20m{logging.getLevelName(logging.WARNING)}\033[1;0m')
    logging.addLevelName(logging.CRITICAL, f'\033[1m\x1b[31;20m{logging.getLevelName(logging.CRITICAL)}\033[1;0m')
else:
    tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

log = logging.getLogger(__name__)