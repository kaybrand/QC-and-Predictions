"""Idempotent uploader backbone for the IGVF Data Portal metadata tables.

Runs alongside (not instead of) the Synapse upload system in
manage_synapse_manifest.py -- separate destination, separate state, same
UPLOAD_ELIGIBLE_CLUSTERS feed from resolve_exclusions.py. See registry.py,
state.py, and orchestrator.py for the actual mechanism; tables/ holds one
module per IGVF metadata table.
"""
