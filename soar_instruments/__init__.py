# Import the modules under this package, to trigger any class
# registering that may be neededfrom astrodata import version

from astrodata import version
__version__ = version()

import astrodata

from . import soar

from . import sami
from . import goodman