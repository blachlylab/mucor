import unittest
import io
import os
from pdb import set_trace as stop
import HTSeq

# mucor modules
import config
import mucorfilters as mf

class TestMucorFilterMethods(unittest.TestCase):
    def setUp(self):
        self.mFilters = mf.MucorFilters()

    def tearDown(self):
        self.mFilters = None

    def test_parseRegionBed(self):
        f_bed = os.path.join("tests", "test.bed")
        self.mFilters.parseRegionBed(f_bed)
        # check region dict keys
        self.assertIn("chr2", self.mFilters.regionDict)
        self.assertIn("chr3", self.mFilters.regionDict)
        self.assertIn("chr4", self.mFilters.regionDict)
        self.assertIn("chr7", self.mFilters.regionDict)
        # check region dict values
        for r_bed in [("127473530", "127474697", "Pos3"),
                      ("127471196", "127472363", "Pos1"),
                      ("127472363", "127473530", "Pos2")]:
            self.assertIn(r_bed, self.mFilters.regionDict["chr7"])
        self.assertIn(("127434530", "127436530", "Pos6"), self.mFilters.regionDict["chr4"])
        self.assertIn(("127477530", "127479530", "Pos5"), self.mFilters.regionDict["chr3"])
        self.assertIn(("127472530", "127474530", "Pos4"), self.mFilters.regionDict["chr2"])

    def test_generateRegionFilterKeys(self):
        f_bed = os.path.join("tests", "test.bed")
        regions = [f_bed, "chr1:1234567-1235467", "chr5"]
        self.mFilters.generateRegionFilter(regions)
        # check region dict keys
        self.assertIn("chr1", self.mFilters.regionDict)
        self.assertIn("chr2", self.mFilters.regionDict)
        self.assertIn("chr3", self.mFilters.regionDict)
        self.assertIn("chr4", self.mFilters.regionDict)
        self.assertIn("chr5", self.mFilters.regionDict)
        self.assertIn("chr7", self.mFilters.regionDict)

    def test_generateRegionFilterValues(self):
        f_bed = os.path.join("tests", "test.bed")
        regions = [f_bed, "chr1:1234567-1235467", "chr5"]
        self.mFilters.generateRegionFilter(regions)
        # check region dict values
        '''
        ("chr7", 127471196, 127472363),
        ("chr7", 127472363, 127473530),
        ("chr7", 127473530, 127474697),
        ("chr4", 127434530, 127436530),
        ("chr3", 127477530, 127479530),
        ("chr2", 127472530, 127474530),
        ("chr1", 1234567, 1235467),
        ("chr5", 0, 0)
        '''
        for r_bed in [("127473530", "127474697", "Pos3"),
                      ("127471196", "127472363", "Pos1"),
                      ("127472363", "127473530", "Pos2")]:
            self.assertIn(r_bed, self.mFilters.regionDict["chr7"])
        self.assertIn(("127434530", "127436530", "Pos6"), self.mFilters.regionDict["chr4"])
        self.assertIn(("127477530", "127479530", "Pos5"), self.mFilters.regionDict["chr3"])
        self.assertIn(("127472530", "127474530", "Pos4"), self.mFilters.regionDict["chr2"])
        self.assertIn(("1234567", "1235467"), self.mFilters.regionDict["chr1"])
        self.assertIn((0, 0), self.mFilters.regionDict["chr5"])

    def test_filterLocPass(self):
        f_bed = os.path.join("tests", "test.bed")
        regions = [f_bed, "chr1:1234567-1235467", "chr5"]
        self.mFilters.generateRegionFilter(regions)
        for roi in [("chr1", 1234568, 1234569),
                    ("chr2", 127472531, 127472533),
                    ("chr3", 127479520, 127479528),
                    ("chr4", 127434533, 127434535),
                    ("chr5", 1, 2),
                    ("chr7", 127473530, 127474696),
                    ("chr7", 127471198, 127471198),
                    ("chr7", 127472366, 127472366)]:
            self.assertFalse(self.mFilters.filterLoc(roi[0], roi[1], roi[2]))

    def test_filterLocFail(self):
        f_bed = os.path.join("tests", "test.bed")
        regions = [f_bed, "chr1:1234567-1235467", "chr5"]
        self.mFilters.generateRegionFilter(regions)
        for roi in [("chr7", 127471190, 127471195), # left -window
                    ("chr7", 127474700, 127474710), # right -window
                    ("chr7", 127474690, 127474710), # straddles -window
                    ("chr4", 127434528, 127434528), # left -piont
                    ("chr3", 127479532, 127479532), # right -point
                    ("chr2", 127472529, 127472529), # ledge -piont
                    ("chr1", 1235468, 1235468)]:    # redge -point
            self.assertTrue(self.mFilters.filterLoc(roi[0], roi[1], roi[2]))

    def test_filterVCFRowPass(self):
        self.mFilters.vcfFilters = ["PASS",".","aD","lowBQ"]
        row = HTSeq.VariantCall()
        row.pos = HTSeq.GenomicPosition("chr1", 1234568)
        # vcf tests
        row.filter = "PASS"
        self.assertFalse(self.mFilters.filterVCFRow(row, "vcf"))

        row.filter = "aD"
        self.assertFalse(self.mFilters.filterVCFRow(row, "vcf.gz"))

        row.filter = "aD;lowBQ"
        self.assertFalse(self.mFilters.filterVCFRow(row, "vcf"))

        # mutect.out tests
        row = {}
        fieldId = {"judgement":1}
        row[fieldId['judgement']] = "lowBQ"
        self.assertFalse(self.mFilters.filterVCFRow(row, "out", fieldId))

    def test_filterVCFRowFail(self):
        self.mFilters.vcfFilters = ["PASS",".","aD"]
        row = HTSeq.VariantCall()
        row.pos = HTSeq.GenomicPosition("chr1", 1234568)
        # vcf tests
        row.filter = "lowBQ"
        self.assertTrue(self.mFilters.filterVCFRow(row, "vcf.gz"))

        row.filter = "aD;lowBQ"
        self.assertTrue(self.mFilters.filterVCFRow(row, "vcf"))

        # mutect.out tests
        row = {}
        fieldId = {"judgement":1}
        row[fieldId['judgement']] = "lowBQ"
        self.assertTrue(self.mFilters.filterVCFRow(row, "out", fieldId))

        # location out of range
        f_bed = os.path.join("tests", "test.bed")
        regions = [f_bed, "chr1:1234567-1235467", "chr5"]
        self.mFilters.generateRegionFilter(regions)
        row = HTSeq.VariantCall()
        row.pos = HTSeq.GenomicPosition("chr1", 1234566)
        row.filter = "PASS"
        self.assertTrue(self.mFilters.filterVCFRow(row, "vcf"))

class TestConfigMethods(unittest.TestCase):
    def setUp(self):
        self.config = config.Config()

    def tearDown(self):
        self.config = None

    def test_parseJSON(self):
        f_json = os.path.join("tests","test.json")
        self.config.parseJSON(f_json)

        # output formats are de-duped
        self.assertEqual(len(self.config.outputFormats), 11)
        # bam files are ignored, vcf and mutect are accepted
        self.assertTrue("/home/path/r1.bam" not in self.config.inputFiles)
        self.assertTrue("/home/path/sample2.out" in self.config.inputFiles)
        self.assertTrue("/home/path/sample2.out.vcf" in self.config.inputFiles)
        # sample name lookup dict works - not relevant for multi-sample VCFs
        self.assertEqual(self.config.filename2samples["sample1.out.vcf"], "Sample1")
        self.assertEqual(self.config.filename2samples["sample2.out.vcf"], "Sample2")
        self.assertEqual(self.config.filename2samples["sample2.out"], "Sample2")


    # def test_upper(self):
    #     self.assertEqual('foo'.upper(), 'FOO')

    # def test_isupper(self):
    #     self.assertTrue('FOO'.isupper())
    #     self.assertFalse('Foo'.isupper())

    # def test_split(self):
    #     s = 'hello world'
    #     self.assertEqual(s.split(), ['hello', 'world'])
    #     # check that s.split fails when the separator is not a string
    #     with self.assertRaises(TypeError):
    #         s.split(2)

if __name__ == '__main__':
    unittest.main()