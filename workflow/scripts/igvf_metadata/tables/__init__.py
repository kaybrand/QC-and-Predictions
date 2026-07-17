"""Importing this package registers every table module with igvf_metadata.registry.
Add one `from . import <module>` line per new table as it's designed."""

from . import prediction_tabular_files  # noqa: F401
from . import signal_files  # noqa: F401
from . import bedpe_index_file  # noqa: F401
from . import prediction_set  # noqa: F401
from . import principal_pseudobulk_set  # noqa: F401
from . import documents  # noqa: F401
