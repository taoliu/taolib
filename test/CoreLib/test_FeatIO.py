#!/usr/bin/env python

from Cistrome.CoreLib.FeatIO import SOFT

import unittest

class TestSOFT(unittest.TestCase):
    def setUp(self):
        self.softstring = """
GEO SOFT Filed Name	 Protocaol Name	 Protocol Content
 The Experiment is AB1791_H3_GLP1TS_AD_1_EVERETT
!Sample_amplification_protocol_ch1	Worm_LM-PCR_Amplification_for_ChIP-chip_v1	Worm_LM-PCR_Amplification_for_ChIP-chip_v1. ChIP DNA was amplified with a modified ligation mediated PCR (LM-PCR) protocol derived from Ren, R et al (2000) Science 290, 2306-9.
!Sample_amplification_protocol_ch2	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1. 1/3 of ChIP and 10ng of input are blunted (T4 polymerase), ligated (concentrated T4 ligase) to annealed linkers and amplified by PCR using longer oligonucleotide as a primer. Generally two rounds of amplification are used to get the amount needed for microarray. Amplified DNA is tested by q-PCR and DNA gel is run to ensure small size and lack of degradation.
!Sample_amplification_protocol_comment_ch1		
!Sample_amplification_protocol_comment_ch2		
!Sample_biomaterial_provider_ch1		
!Sample_biomaterial_provider_ch2		
!Sample_characteristics_ch1		
!Sample_characteristics_ch2		
!Sample_data_processing	ChIP-chip_normalization_standard_MA2C_v1	ChIP-chip_normalization_standard_MA2C_v1. First, all the IP and INPUT log ratio values are read from the pairdata file. Secondly, we build GC bins for INPUT and IP based on the GC counts for every probe sequence, which means the INPUT or IP values for any probes who have the same GC counts will be put together. After that, for each GC bin, we calculate the mean for IP and INPUT data, and the covariance between this two channels. By default, the robust mean variance method is applied, which generalizes Tukeys theory of bi-weight estimation where the constant C is set to 2. At last, we adjust the log ratio values for each probe by using the mean and covariance values for their corresponding GC bins, then these values are further normalized by their mean and standard derivation. In case of replicates, when we calculate the MA2Cscore afterwards, we take the median as the score from all the replicates for all the probes within the sliding window defined by bandwidth parameter.
!Sample_description		
!Sample_extract_protocol_ch1	Worm_chromatin_immunoprecipitation_vCW1	Worm_chromatin_immunoprecipitation_vCW1. 2mg extract was used for each ChIP with 5% taken as input directly into elution buffer (1% SDS in TE, 0.1M NaHCO3). Antibody was added to each IP sample and incubated overnight at 4C. Immune complexes were incubated (2hrs at 4C) with 10 ul of protein A sepharose, and washed 5 minutes each with 1.4 mL of each of the following solutions: ChIP Buffer, ChIP Buffer+500mM NaCl, ChIP Buffer+1M NaCl, LiCl solution (10mM Tris-HCl pH 8.0, 250mM LiCl, 0.5% NP-40, 0.5% sodium deoxycholate, 1mM EDTA), and 1X TE (10mM Tris-HCl pH 8.0, 1mM EDTA) and treated with 20ug RNase A for 30 minutes. Samples were then washed once with 1X TE and eluted twice with 200uL elution buffer for 15minutes. 16ul 5M NaCl was added to each sample then transferred to 65C overnight to reverse crosslinks. DNA was cleaned up with a Zymo DNA clean up kit. For a detailed protocol see http://www.modencode.org/.
!Sample_extract_protocol_ch2		
!Sample_extract_protocol_comment_ch1		
!Sample_extract_protocol_comment_ch2		
!Sample_geo_accession		
!Sample_growth_protocol_ch1	Worm_L3_growth_and_harvest_vPK1	Worm_L3_growth_and_harvest_vPK1. About 2-7 million of worms are bleached and then hatched in M9 for 24-42 hrs. About 100 embryos are seeded onto the plate to test for contamination and hatching efficiency. Remaining hatched L1 larvae are inoculated in a proper volume of liquid culture. Next day when larvae reach the L3 stage they are cleaned by M9 washes and sucrose gradient and collected by freezing in liquid nitrogen. Just before collection DIC pictures are taken and about 50ul of worms are stained for DAPI to assess the stage.
!Sample_growth_protocol_ch2	Worm_embryo_growth_and_harvest_v1	Worm_embryo_growth_and_harvest_v1. Embryos were prepared by bleaching from gravid N2 adults grown in standard S-basal media liquid culture. Live embryos were cross-linked in M9 + 2% formaldehyde for 30 minutes at room temperature followed by quenching with 125mM glycine for 5 minutes. Embryos were then washed twice with M9 Buffer and once by FA buffer (50 mM HEPES/KOH pH 7.5, 1 mM EDTA, 1% Triton X-100, 0.1 % sodium deoxycholate; 150 mM NaCl). Pellets were frozen at -80C. For a detailed protocol see http://www.modencode.org/.
!Sample_growth_protocol_comment_ch1		
!Sample_growth_protocol_comment_ch2		
!Sample_hyb_protocol	ChIP-chip_label_hyb_nimblegen_v1	ChIP-chip_label_hyb_nimblegen_v1. DNA was labeled and hybridized to C. elegans tiling array by Roche NimbleGen according to the protocol described in chapter 3 and 4 of the NimbleGen Arrays User?s Guide ChIP-chip Analysis, Version 3.1, 27 May 2008. Briefly, Amplified IP or input DNA was either labeled with Cy5 or Cy3 in the presence of Klenow fragment. The reaction was stopped by the addition of EDTA. Labeled DNA was recovered by isopropanol precipitation, and dried. The labeled DNA was hybridized to C. elegans tiling array for 16 - 20 hours at 42?C.
!Sample_label_ch1	Cy3	Cy3
!Sample_label_ch2	Cy5	Cy5
!Sample_molecule_ch1	genomic DNA	genomic DNA
!Sample_molecule_ch2	mRNA	mRNA
!Sample_organism_ch1	Caenorhabditis elegans	Caenorhabditis elegans
!Sample_organism_ch2	Caenorhabditis elegans	Caenorhabditis elegans
!Sample_platform_id		
!Sample_scan_protocol	ChIP-chip_scanning_nimblegen_v1	ChIP-chip_scanning_nimblegen_v1. Array scanning and raw data extraction were performed at Roche NimbleGen, according to the protocol described in chapter 5 and 6 of the NimbleGen Arrays User?s Guide ChIP-chip Analysis, Version 3.1, 27 May 2008. Briefly, array signal was scanned by using a GenePix 4000B Scanner with associated software and saved as .tif files of the 532nm and 635nm images individually. Raw signal intensities of the images were extracted  and saved as .pair files by using NimbleScan software according to the NimbleScan v2.4 User?s Guide.
!Sample_source_name_ch1		
!Sample_source_name_ch2		
!Sample_title	AB1791_H3_GLP1TS_AD_1_EVERETT	
!Sample_treatment_protocol_ch1	Worm_L3_extraction_vPK1	Worm_L3_extraction_vPK1. Worms are frozen, ground, and crosslinked for 10 minutes in 1% formaldehyde. Later, washed pellets are resuspended in FA buffer and subjected to sonication in Bioruptor (14 pulses of 30 seconds with 1 minute rests in between). Extracts are then spun down and soluble fraction is stored for quality tests and future ChIP.
!Sample_treatment_protocol_ch2	Worm_embryo_extraction_v1	Worm_embryo_extraction_v1. Embryos were resuspended in FA buffer (50 mM HEPES/KOH pH 7.5, 1 mM EDTA, 1% Triton X-100, 0.1 % sodium deoxycholate; 150 mM NaCl) + protease inhibitors (Calbiochem Cat# 539131). Using a Branson sonifier microtip, samples were sonicated on ice at the following settings: 35% amplitude, 0.9 sec on, 0.1 sec off, 12 pulses, 7 times. Cell debris was removed by centrifuging at 13,000 g for 15 minutes at 4?C and taking the supernatant. Protein concentration was determined by Bradford Assay and extracts were aliquoted at stored at -80C. For a detailed protocol see http://www.modencode.org/.
!Sample_treatment_protocol_comment_ch1		
!Sample_treatment_protocol_comment_ch2		
^SAMPLE	AB1791_H3_GLP1TS_AD_1_EVERETT	
 The Experiment is AB1791_H3_N2_L3_1LM
!Sample_amplification_protocol_ch1	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1. 1/3 of ChIP and 10ng of input are blunted (T4 polymerase), ligated (concentrated T4 ligase) to annealed linkers and amplified by PCR using longer oligonucleotide as a primer. Generally two rounds of amplification are used to get the amount needed for microarray. Amplified DNA is tested by q-PCR and DNA gel is run to ensure small size and lack of degradation.
!Sample_amplification_protocol_ch2	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1	Worm_LM-PCR_Amplification_for_ChIP-chip_vPK1. 1/3 of ChIP and 10ng of input are blunted (T4 polymerase), ligated (concentrated T4 ligase) to annealed linkers and amplified by PCR using longer oligonucleotide as a primer. Generally two rounds of amplification are used to get the amount needed for microarray. Amplified DNA is tested by q-PCR and DNA gel is run to ensure small size and lack of degradation.
!Sample_amplification_protocol_comment_ch1		
!Sample_amplification_protocol_comment_ch2		
!Sample_antibody_name	AB1791_H3	An affinity purified rabbit polyclonal antibody to H3 obtained from Abcam (H3-AB1791);used for ChIP.
!Sample_biomaterial_provider_ch1		
!Sample_biomaterial_provider_ch2		
!Sample_characteristics_ch1		
!Sample_characteristics_ch2		
!Sample_data_processing	ChIP-chip_normalization_standard_MA2C_v1	ChIP-chip_normalization_standard_MA2C_v1. First, all the IP and INPUT log ratio values are read from the pairdata file. Secondly, we build GC bins for INPUT and IP based on the GC counts for every probe sequence, which means the INPUT or IP values for any probes who have the same GC counts will be put together. After that, for each GC bin, we calculate the mean for IP and INPUT data, and the covariance between this two channels. By default, the robust mean variance method is applied, which generalizes Tukeys theory of bi-weight estimation where the constant C is set to 2. At last, we adjust the log ratio values for each probe by using the mean and covariance values for their corresponding GC bins, then these values are further normalized by their mean and standard derivation. In case of replicates, when we calculate the MA2Cscore afterwards, we take the median as the score from all the replicates for all the probes within the sliding window defined by bandwidth parameter.
!Sample_description		
!Sample_extract_protocol_ch1	Worm_chromatin_immunoprecipitation_vPK1	Worm_chromatin_immunoprecipitation_vPK1. Appropriate amount of extract is incubated overnight with a proper amount of antibody (exceptional antibodies due to better results are incubated 2hrs). Afterwards, 40ul of equilibrated magnetic beads (either protein A or G, depending on antibody) are added and incubated for 2 hrs. Later, washes with FA, 500mM-salt FA, 1M salt FA, TEL, and TE buffer are performed and DNA is eluted in elution buffer (1% SDS in TE with 250 mM NaCl) ? two times with 57 ml volume each, at 65?C. Samples are treated with RNAse, proteinase K and then crosslinks are reversed overnight at 65?C. DNA is purified on qiagen PCR purification columns, tested by q-PCR for ChIP quality, and stored in -20?C for future applications and.
!Sample_extract_protocol_ch2	Worm_chromatin_immunoprecipitation_vPK1	Worm_chromatin_immunoprecipitation_vPK1. Appropriate amount of extract is incubated overnight with a proper amount of antibody (exceptional antibodies due to better results are incubated 2hrs). Afterwards, 40ul of equilibrated magnetic beads (either protein A or G, depending on antibody) are added and incubated for 2 hrs. Later, washes with FA, 500mM-salt FA, 1M salt FA, TEL, and TE buffer are performed and DNA is eluted in elution buffer (1% SDS in TE with 250 mM NaCl) ? two times with 57 ml volume each, at 65?C. Samples are treated with RNAse, proteinase K and then crosslinks are reversed overnight at 65?C. DNA is purified on qiagen PCR purification columns, tested by q-PCR for ChIP quality, and stored in -20?C for future applications and.
!Sample_extract_protocol_comment_ch1		
!Sample_extract_protocol_comment_ch2		
!Sample_geo_accession		
!Sample_growth_protocol_ch1	Worm_L3_growth_and_harvest_vPK1	Worm_L3_growth_and_harvest_vPK1. About 2-7 million of worms are bleached and then hatched in M9 for 24-42 hrs. About 100 embryos are seeded onto the plate to test for contamination and hatching efficiency. Remaining hatched L1 larvae are inoculated in a proper volume of liquid culture. Next day when larvae reach the L3 stage they are cleaned by M9 washes and sucrose gradient and collected by freezing in liquid nitrogen. Just before collection DIC pictures are taken and about 50ul of worms are stained for DAPI to assess the stage.
!Sample_growth_protocol_ch2	Worm_L3_growth_and_harvest_vPK1	Worm_L3_growth_and_harvest_vPK1. About 2-7 million of worms are bleached and then hatched in M9 for 24-42 hrs. About 100 embryos are seeded onto the plate to test for contamination and hatching efficiency. Remaining hatched L1 larvae are inoculated in a proper volume of liquid culture. Next day when larvae reach the L3 stage they are cleaned by M9 washes and sucrose gradient and collected by freezing in liquid nitrogen. Just before collection DIC pictures are taken and about 50ul of worms are stained for DAPI to assess the stage.
!Sample_growth_protocol_comment_ch1		
!Sample_growth_protocol_comment_ch2		
!Sample_hyb_protocol	ChIP-chip_label_hyb_nimblegen_v1	ChIP-chip_label_hyb_nimblegen_v1. DNA was labeled and hybridized to C. elegans tiling array by Roche NimbleGen according to the protocol described in chapter 3 and 4 of the NimbleGen Arrays User?s Guide ChIP-chip Analysis, Version 3.1, 27 May 2008. Briefly, Amplified IP or input DNA was either labeled with Cy5 or Cy3 in the presence of Klenow fragment. The reaction was stopped by the addition of EDTA. Labeled DNA was recovered by isopropanol precipitation, and dried. The labeled DNA was hybridized to C. elegans tiling array for 16 - 20 hours at 42?C.
!Sample_label_ch1	Cy3	Cy3
!Sample_label_ch2	Cy5	Cy5
!Sample_molecule_ch1	genomic DNA	genomic DNA
!Sample_molecule_ch2	genomic DNA	genomic DNA
!Sample_organism_ch1	Caenorhabditis elegans	Caenorhabditis elegans
!Sample_organism_ch2	Caenorhabditis elegans	Caenorhabditis elegans
!Sample_platform_id		
!Sample_scan_protocol	ChIP-chip_scanning_nimblegen_v1	ChIP-chip_scanning_nimblegen_v1. Array scanning and raw data extraction were performed at Roche NimbleGen, according to the protocol described in chapter 5 and 6 of the NimbleGen Arrays User?s Guide ChIP-chip Analysis, Version 3.1, 27 May 2008. Briefly, array signal was scanned by using a GenePix 4000B Scanner with associated software and saved as .tif files of the 532nm and 635nm images individually. Raw signal intensities of the images were extracted  and saved as .pair files by using NimbleScan software according to the NimbleScan v2.4 User?s Guide.
!Sample_source_name_ch1		
!Sample_source_name_ch2		
!Sample_title	AB1791_H3_N2_L3_1LM	
!Sample_treatment_protocol_ch1	Worm_L3_extraction_vPK1	Worm_L3_extraction_vPK1. Worms are frozen, ground, and crosslinked for 10 minutes in 1% formaldehyde. Later, washed pellets are resuspended in FA buffer and subjected to sonication in Bioruptor (14 pulses of 30 seconds with 1 minute rests in between). Extracts are then spun down and soluble fraction is stored for quality tests and future ChIP.
!Sample_treatment_protocol_ch2	Worm_L3_extraction_vPK1	Worm_L3_extraction_vPK1. Worms are frozen, ground, and crosslinked for 10 minutes in 1% formaldehyde. Later, washed pellets are resuspended in FA buffer and subjected to sonication in Bioruptor (14 pulses of 30 seconds with 1 minute rests in between). Extracts are then spun down and soluble fraction is stored for quality tests and future ChIP.
!Sample_treatment_protocol_comment_ch1		
!Sample_treatment_protocol_comment_ch2		
^SAMPLE	AB1791_H3_N2_L3_1LM
"""
    def test_parseSOFT (self):
        entities = SOFT.parsestring(self.softstring)
        self.assertEqual(entities[0].entity_dict["SAMPLE"],"AB1791_H3_GLP1TS_AD_1_EVERETT")
        self.assertEqual(entities[1].entity_dict["SAMPLE"],"AB1791_H3_N2_L3_1LM")


if __name__ == '__main__':
    unittest.main()
