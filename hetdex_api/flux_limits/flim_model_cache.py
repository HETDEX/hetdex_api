"""

Module to cache the flux limit model
to avoid having to repeatedly initialise it
when looping over sensitivity cubes

"""

cached_sim_interp = None
cached_model = None


