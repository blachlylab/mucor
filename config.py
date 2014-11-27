# config.py
# 
# James S Blachly MD
# 2014-10-04

class Config(object):
    """Object to hold configuration options passed from cmdline or parsed from JSON config file"""

    def __init__(self):
            self.outputDir = ''
            self.gff = ''
            self.union = False
            self.fast = False
            self.featureType = ''
            self.filters = []
            self.outputFormats = []
            self.samples = []
            self.inputFiles = []    # derived from samples
            self.database = []
            self.regions = []

    def __str__(self):
        string_rep =    "No. samples: {0}\n\n" + \
                        "JSON Options:\n" + \
                        "\tFeature Type: {1}\n" + \
                        "\tOutput Dir: {2}\n" + \
                        "\tUnion: {3}\n" + \
                        "\tFast: {4}\n" + \
                        "\tAnnotation: {5}\n" + \
                        "\tDatabase(s): \n\t\t" + \
                        "\n\t\t".join(self.database) + "\n" + \
                        "\tFilter(s): \n\t\t" + \
                        "\n\t\t".join(self.filters) + "\n" + \
                        "\tInput Files: \n\t\t" + \
                        "\n\t\t".join(self.inputFiles) + "\n" + \
                        "\tOutput Formats: \n\t\t" + \
                        "\n\t\t".join(self.outputFormats) + "\n"
        string_rep = string_rep.format(str(len(self.inputFiles)), self.featureType, \
                                      self.outputDir, self.union, self.fast, self.gff)
        return string_rep
