# config.py
# 
# James S Blachly MD
# 2014-10-04

class Config(object):
    """Object to hold configuration options passed from cmdline or parsed from JSON config file"""

    def __init__(self):
            # run_name = ''   # As named in JSON
            outputDir = '' # As used in mucor.py -- TO DO: need consistency between config options and var names to facilitate programmatic introspection etc.
            gtf = ''
            union = False
            fast = False
            feature = ''
            filters = []
            outputFormats = []
            samples = []
            inputFiles = []    # derived from samples
            database = []
            regions = []
