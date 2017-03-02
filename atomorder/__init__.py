
# done this way to ensure the loading in the right order
from .globals import settings

from .definitions import Constants
# constants defined in the Constants class
constants = Constants()

from .structure import Reaction
from .utils import parse_args

