#! /usr/bin/env python2

import vcf
import vcf.utils
import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
from bioutils.assemblies import make_name_ac_map

__author__ = 'Claudia Juno'


class Assignment3:
    
    def __init__(self, f, m, s):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        ## Check if hgvs is installed
        print("HGVS version: %s" % hgvs.__version__)

        self.father = f
        self.mother = m
        self.son = s

    def get_total_number_of_variants_mother(self):
        vcf_reader = vcf.Reader(open(self.mother, 'r'))
        z = 0
        for r in vcf_reader:
            z += 1
        print (z)
        
    def get_total_number_of_variants_father(self):
        vcf_reader = vcf.Reader(open(self.father, 'r'))
        z = 0
        for r in vcf_reader:
            z += 1
        print (z)
       
        
    def get_variants_shared_by_father_and_son(self):
        vcf_reader_f = vcf.Reader(open(self.father, 'r'))
        vcf_reader_s = vcf.Reader(open(self.son, 'r'))
        z = 0
        for r in vcf.utils.walk_together(vcf_reader_f, vcf_reader_s):
            if (r[0] == r[1]) and not (r[0] == None) and not (r[1] == None):
                # falls neben der Anzahl auch die Varianten gefragt sind, folgende Zeile einkommentieren.
                #print (r[0])
                z += 1
        print (z)


    def get_variants_shared_by_mother_and_son(self):
        vcf_reader_m = vcf.Reader(open(self.mother, 'r'))
        vcf_reader_s = vcf.Reader(open(self.son, 'r'))
        z = 0
        for r in vcf.utils.walk_together(vcf_reader_m, vcf_reader_s):
            if (r[0] == r[1]) and not (r[0] == None) and not (r[1] == None):
                # falls neben der Anzahl auch die Varianten gefragt sind, folgende Zeile einkommentieren.
                # print (r[0])
                z += 1
        print (z)


    def get_variants_shared_by_trio(self):
        vcf_reader_f = vcf.Reader(open(self.father, 'r'))
        vcf_reader_m = vcf.Reader(open(self.mother, 'r'))
        vcf_reader_s = vcf.Reader(open(self.son, 'r'))
        z = 0
        for r in vcf.utils.walk_together(vcf_reader_f, vcf_reader_m, vcf_reader_s):
            if (r[0] == r[1] == r[2]) and not (r[0] == None) and not (r[1] == None) and not (r[2] == None):
                # falls neben der Anzahl auch die Varianten gefragt sind, folgende Zeile einkommentieren.
                # print (r[0])

                z += 1
        print (z)
        

    def merge_mother_father_son_into_one_vcf(self):
        vcf_reader_f = vcf.Reader(open(self.father, 'r'))
        vcf_reader_m = vcf.Reader(open(self.mother, 'r'))
        vcf_reader_s = vcf.Reader(open(self.son, 'r'))

        mf = open('merge_mother_father_son.vcf', 'w')
        writer = vcf.Writer(mf, vcf_reader_s)
        for r in vcf.utils.walk_together(vcf_reader_f, vcf_reader_m, vcf_reader_s):
            if (r[0] == r[1] == r[2]) and not (r[0] == None) and not (r[1] == None) and not (r[2] == None):
                writer.write_record(r[0])
        mf.close()
        print ("find result in \"merge_mother_father_son.vcf\"")


    def convert_first_variants_of_son_into_HGVS(self):

        z = 0      # zaehler fuer 100 variants
        z_ok = 0   # zaehler fuer erfolgreiche conversions
        z_exceptions = 0    # zaehler fuer exceptions

        ## Connect to UTA
        hdp = hgvs.dataproviders.uta.connect()

        assembly_mapper = hgvs.assemblymapper.AssemblyMapper(hdp, normalize=False)  # EasyVariantMapper before, normalize=False, um Warning zu beseitigen
        ## Used for parsing
        hgvsparser = hgvs.parser.Parser()  # Parser

        vcf_reader_s = vcf.Reader(open(self.son, 'r'))    # reader wie oben

        for r in vcf_reader_s:
            if z < 100:
                refseq_nc_number = make_name_ac_map("GRCh37.p13")[r.CHROM[3:]]
                genome_hgvs = "%s:g.%s%s>%s" % (refseq_nc_number, str(r.POS), str(r.REF), str(r.ALT[0]))

                try:
                    g = hgvsparser.parse_hgvs_variant(genome_hgvs)
                    for t in assembly_mapper.relevant_transcripts(g):
                        try:
                            c = assembly_mapper.g_to_c(g, t)  # c: coding DNA reference sequence, g: genomic reference sequence
                            z_ok += 1
                            print ("%s\t%s" % (g, c))
                        except hgvs.exceptions.HGVSUsageError:
                            n = assembly_mapper.g_to_n(g, t)  # n: non-coding RNA reference sequence (gene producing an RNA transcript but not a protein)
                            z_ok += 1
                            print ("%s\t%s" % (g, n))
                        except:
                            z_exceptions += 1
                except Exception:
                    z_exceptions += 1

            else:
                break

            z += 1

        # Summary ausgeben
        print ("Successful conversions: %s" % (z_ok))
        print ("Exceptions occurred: %s" % (z_exceptions))




        
    
    def print_summary(self):
        print ("-----------------------")
        print ("total variants number mother:")
        self.get_total_number_of_variants_mother()
        print("-----------------------")
        print ("total variants number father:")
        self.get_total_number_of_variants_father()
        print("-----------------------")
        print ("total variants shared by father & son:")
        self.get_variants_shared_by_father_and_son()
        print("-----------------------")
        print ("total variants shared by mother & son:")
        self.get_variants_shared_by_mother_and_son()
        print("-----------------------")
        print ("total variants shared by Trio:")
        self.get_variants_shared_by_trio()
        print("-----------------------")
        print ("merge vcf of father, mother & son in one file:")
        self.merge_mother_father_son_into_one_vcf()
        print("-----------------------")
        print("convert first 100 variants of son into HGVS:")
        self.convert_first_variants_of_son_into_HGVS()
        
if __name__ == '__main__':
    print("Assignment 3")
    print(__author__)
    assignment1 = Assignment3("AmpliseqExome.20141120.NA24149.vcf", "AmpliseqExome.20141120.NA24143.vcf", "AmpliseqExome.20141120.NA24385.vcf" )
    assignment1.print_summary()
    
    

