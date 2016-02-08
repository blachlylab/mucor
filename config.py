#    Copyright 2013-2015 James S Blachly, MD and The Ohio State University
#
#    This file is part of Mucor.
#
#    Mucor is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Mucor is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Mucor.  If not, see <http://www.gnu.org/licenses/>.

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
            self.databases = {}
            self.regions = []
            self.filename2samples = {}
            self.source = {}
            self.format = {}

    def __str__(self):
        string_rep =    "No. samples: {0}\n\n" + \
                        "JSON Options:\n" + \
                        "\tFeature Type: {1}\n" + \
                        "\tOutput Dir: {2}\n" + \
                        "\tUnion: {3}\n" + \
                        "\tFast: {4}\n" + \
                        "\tAnnotation: {5}\n" + \
                        "\tDatabase(s): \n\t\t" + \
                        "\n\t\t".join(self.databases) + "\n" + \
                        "\tFilter(s): \n\t\t" + \
                        "\n\t\t".join(self.filters) + "\n" + \
                        "\tInput Files: \n\t\t" + \
                        "\n\t\t".join(self.inputFiles) + "\n" + \
                        "\tOutput Formats: \n\t\t" + \
                        "\n\t\t".join(self.outputFormats) + "\n"
        string_rep = string_rep.format(str(len(self.inputFiles)), self.featureType, \
                                      self.outputDir, self.union, self.fast, self.gff)
        return string_rep
