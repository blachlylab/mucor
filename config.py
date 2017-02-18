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

import os
import sys
import json

# mucor modules
import info
from databases import checkAndOpen
import output

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
        self.total = 0

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

    def parseJSON(self, json_config):
        '''Import the JSON config file from mucor_config.py.
        Reads the config file into a dictionary, then writes each dictionary entry into the respective Config class position.
        '''

        try:
            JD = json.load(open(json_config, 'r'))
        except ValueError as json_error:
            info.throwWarning(json_error.message)
            info.abortWithMessage("Could not load the given JSON config file. See the example config for proper formatting.")

        # write dictionary values into a more human-friendly config class
        self.featureType = JD['feature']
        self.outputDir = os.path.expanduser(JD['outputDir'])
        self.union = JD['union']
        self.fast = JD['fast']
        self.gff = JD['gff']
        self.outputFormats = list(set(JD['outputFormats'])) # 'set' prevents repeated formats from being written multiple times
        if JD['databases']:
            if 'tabix' in sys.modules: # make sure tabix is imported 
                for name,db in JD['databases'].items():
                    dbPointer = checkAndOpen(db)
                    if dbPointer:
                        #check for non-null pointers
                        self.databases[name] = dbPointer
            else: # the user supplied databases but did not successfully import tabix
                info.throwWarning("tabix module not found; database features disabled")
        if str(JD['regions']):
            self.regions = JD['regions']
        else:
            self.regions = []
        # comma separated list of acceptable VCF filter column values
        self.filters = JD['filters']

        self.inputFiles = []
        self.samples = []
        for i in JD['samples']:
            self.samples.append(i['id'])
            for j in i['files']:
                if j['type'] == "bam":
                    # skip the bam files
                    continue
                filename = os.path.basename(j['path'])
                self.filename2samples[filename] = i['id']
                self.source[filename] = j['source']
                self.inputFiles.append(j['path'])
        return

    def validate(self):
        '''
        determine whether we should proceed with execution
        '''
        if not os.path.exists(self.gff):
            info.abortWithMessage("Could not find GFF file {0}".format(self.gff))
        if os.path.exists(self.outputDir) and [x for x in os.listdir(self.outputDir) if x in output.Writer().file_names.values()]:
            info.abortWithMessage("The directory {0} already exists and contains output. Will not overwrite.".format(self.outputDir))
        elif not os.path.exists(self.outputDir):
            try:
                os.makedirs(self.outputDir)
            except OSError:
                info.abortWithMessage("Error when creating output directory {0}".format(self.outputDir))

        # check that all specified variant files exist
        for fn in self.inputFiles:
            if not os.path.exists(fn):
                info.abortWithMessage("Could not find variant file: {0}".format(fn))

        # Total is used in frequency calculations. 
        #   supplying config.samples will use sample count as the denominator [canonical operation, ie: comparing samples]
        #   or, using samples.inputFiles will use file count [non-canonical operation, ie: comparing tools, or otherwise comparing many vcf files with no regard for sample ID]
        # TODO: consider a flag to switch for this
        self.total = len(set(self.samples))
