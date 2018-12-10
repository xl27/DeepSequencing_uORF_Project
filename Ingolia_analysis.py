####################################
### Python Script for Annotation ###
####################################
#Import Relevent Packages
from __future__ import division
from plastid import GTF2_TranscriptAssembler, BAMGenomeArray, VariableFivePrimeMapFactory
import twobitreader as twobit
import numpy as np
import re
from collections import Counter, defaultdict
import itertools
from plastid.util.io.filters import CommentReader

#Other things we need
import pandas as pd
from plastid import CenterMapFactory

#Read in the data
test = pd.read_csv('Supp3.csv') #Ingolia Harringtonine Peak Annotations
mm9_transcriptome = list(GTF2_TranscriptAssembler('UCSC_mm9_old.gtf')) #MM9 Transcriptome including all 'old transcripts'
genome_twobit=twobit.TwoBitFile('mm9_genome.2bit') #MM9 genome as a .2bit file

#Create a dictionary of transcripts
transcriptome_dict=defaultdict()
for i in mm9_transcriptome:
	name = i.get_name()
	transcriptome_dict[name]=i

#This is actually really simple. If it is the canonical CDS he's annotated it as 'canonical' and the dist to CDS is zero
#If it is an uORF it is annotated as uORF. Don't need to do anything but extract those. 

#Filter Harringtonine Peaks to only contain peaks of interest
elements_of_intrest=['canonical', 'uorf', 'uorf-overlap']
uORF_CDS_peaks=test[np.isin(test.Product, elements_of_intrest)]
uORF_CDS_peaks=uORF_CDS_peaks[~np.isnan(uORF_CDS_peaks['CDS Length [codons]'])] #there were a few Na's that mess things up, removed.

#Set up a big loop over the list of unique transcripts 
unique_genes_toannotate=uORF_CDS_peaks.knownGene.unique() #4071 genes to annotate
sum(np.isin(uORF_CDS_peaks['knownGene'], [i for i in transcriptome_dict.iterkeys()])) #pre-flight check: all are in our annotations

logs=pd.DataFrame(index=range(0,len(unique_genes_toannotate)),columns=['transcript','peaks','error']) #Summary table
cds_list=list() #list of segmentchains
uORF_list=list() #list of segmentchains

#Run the Loop!
for index,transcript in enumerate(unique_genes_toannotate):
    print 'Working on transcript:'
    print transcript 
    print 'Index is: ' 
    print index
    #append into the log
    logs.transcript.iloc[index]=transcript
    #Subset the dataframe to look at that transcript. Get the annotated transcript model
    relevant_peaks=uORF_CDS_peaks[uORF_CDS_peaks.knownGene == transcript]
    transcript_model=transcriptome_dict[transcript]
    #Get the sequence of the transcript
    try:
        rna_seq=transcript_model.get_sequence(genome_twobit)
    #
    except KeyError:
        logs.peaks.iloc[index]='KeyError on Sequence. Check Chrom.'
        continue
    #
    #partition into uORFs and CDS
    annotated_cds=relevant_peaks[relevant_peaks['Dist to CDS [codons]']==0]
    annotated_uORF=relevant_peaks[relevant_peaks['Dist to CDS [codons]']<0]
    #Figure out whether cds and uORF exist
    cds_exists=annotated_cds.shape[0]>0
    uORF_exists=annotated_uORF.shape[0]>0
    #Classify as appropraite
    if (not uORF_exists and not cds_exists):
        logs.peaks.iloc[index]='No CDS or uORF'
        continue
    #
    #Finish our alternative possibilities
    elif (not uORF_exists and cds_exists):
        if relevant_peaks.shape[0]>1:
            logs.error.iloc[index]='More than One CDS'
            continue
        # First Test.
        putative_start=relevant_peaks['Init Codon [nt]'].iloc[0]
        #Test initiation codon
        if rna_seq[putative_start-3:putative_start+4]==relevant_peaks['Init Context [-3 to +4]'].iloc[0]:
            logs.peaks.iloc[index]='CDS only'
            cds_list.append(transcript_model.get_cds())
            continue
        #
        else: 
            logs.peaks.iloc[index]='CDS only'
            logs.error.iloc[index]='More than One CDS'
            continue
        #
    #
    elif (uORF_exists and not cds_exists):
        logs.peaks.iloc[index]='uORF only'
        cds_list.append(transcript_model.get_cds()) #annotate the CDS regardless
        #Loop through possible uORFs. 
        error_count=0
        for z in range(0,annotated_uORF.shape[0]):
            putative_start=annotated_uORF['Init Codon [nt]'].iloc[z]
            if rna_seq[putative_start-3:putative_start+4]==annotated_uORF['Init Context [-3 to +4]'].iloc[z]:
                putative_stop=annotated_uORF['CDS Length [codons]'].iloc[z]*3+3+putative_start #includes stop codon
                discovered_subchain = transcript_model.get_subchain(putative_start,putative_stop) #Cannot simply do my_transcript.get_subchain(start,stop,'ID'='meep') b/c hardcoded that ID of subchain = ID transcript + 'subchain'. This apparently cannot be overwritten
                discovered_subchain.attr['ID'] = transcript_model.get_name() + '_' + str(putative_start) + '_' + str(putative_stop) #overwrite the attribute to identify uORF uniquely
                uORF_list.append(discovered_subchain)
                continue
            #
            else:
                error_count=error_count+1
                continue
            #
        #
        if error_count>0:
            logs.error.iloc[index]='QC For uORFs Failed: ' + str(error_count)
        #
        continue
    #
    elif (uORF_exists and cds_exists):
        logs.peaks.iloc[index]='uORF and CDS'
        if annotated_cds.shape[0]>1:
            logs.error.iloc[index]='More than One CDS'
            continue
        # First Test.
        putative_start=annotated_cds['Init Codon [nt]'].iloc[0]
        #Test initiation codon
        if rna_seq[putative_start-3:putative_start+4]==annotated_cds['Init Context [-3 to +4]'].iloc[0]:
            cds_list.append(transcript_model.get_cds())
        #
        else:
            logs.error.iloc[index]='QC for CDS failed'
            continue
        #Now on to uORFs 
        error_count=0
        for z in range(0,annotated_uORF.shape[0]):
            putative_start=annotated_uORF['Init Codon [nt]'].iloc[z]
            if rna_seq[putative_start-3:putative_start+4]==annotated_uORF['Init Context [-3 to +4]'].iloc[z]:
                putative_stop=annotated_uORF['CDS Length [codons]'].iloc[z]*3+3+putative_start #includes stop codon
                discovered_subchain = transcript_model.get_subchain(putative_start,putative_stop) #Cannot simply do my_transcript.get_subchain(start,stop,'ID'='x') b/c hardcoded that ID of subchain = ID transcript + 'subchain'. This apparently cannot be overwritten
                discovered_subchain.attr['ID'] = transcript_model.get_name() + '_' + str(putative_start) + '_' + str(putative_stop) #overwrite the attribute to identify uORF uniquely
                uORF_list.append(discovered_subchain)
                continue
            #
            else:
                error_count=error_count+1
                continue
            #
        #
        if error_count>0:
            logs.error.iloc[index]='QC For uORFs Failed: ' + str(error_count)
        #
        continue
    #
    else:
        logs.peaks.iloc[index]='error'
    #
#

#Write Python .obj
import dill
dill.dump(uORF_list,open('./updated_quantitation/ingolia_uORFs.obj','wb'))

#Final To-Do: Create Bed and FASTA files
#Write uORFs to file
fout=open('./updated_quantitation/Ingolia_uORFs.bed','w')
for i in uORF_list:
    fout.write(i.as_bed())

fout.close()

#Write uORF FASTA
fout=open('./updated_quantitation/Ingolia_uORFs.fasta','w')
for i in uORF_list:
    fout.write(i.get_fasta(genome_twobit))

fout.close()

#Write CDS to file
fout=open('./updated_quantitation/Ingolia_CDS.bed','w')
for i in cds_list:
    fout.write(i.as_bed())

fout.close()

#Write CDS FASTA
fout=open('./updated_quantitation/Ingolia_CDS.fasta','w')
for i in cds_list:
    fout.write(i.get_fasta(genome_twobit))

fout.close()

#Write Summary Frame
logs.to_csv('./updated_quantitation/Ingolia_Transcripts_Summary.csv')

##Write the -50 to +50 sequences around uORFs to file as 
windows_list = []
for index, i in enumerate(uORF_list):
    print 'Working on ' + i.get_name()
    print 'Index is: ' + str(index)
    #
    roi_sequence = i.get_sequence(genome_twobit)
    full_transcript = transcriptome_dict[i.attr['transcript_id']]
    transcript_sequence = full_transcript.get_sequence(genome_twobit)
    roi_start_in_transcript = transcript_sequence.find(roi_sequence)
    if roi_start_in_transcript - 50 >= 0:
        lower_bound = roi_start_in_transcript - 50
    #
    else:
        lower_bound = 0 
    #
    if roi_start_in_transcript + 50 <= full_transcript.get_length():
        upper_bound = roi_start_in_transcript + 50
    #
    else:
        upper_bound = full_transcript.get_length()
    #
    subchain_window = full_transcript.get_subchain(lower_bound, upper_bound)
    subchain_window.attr['ID'] = i.get_name() + '_' + str('100bpwindow') + '_' + str(roi_start_in_transcript-lower_bound) #overwrite the attribute to identify uORF uniquely
    windows_list.append(subchain_window)

fout=open('./updated_quantitation/uORFs_100bpwindows_aroundstart.fasta','w')
for i in windows_list:
    fout.write(i.get_fasta(genome_twobit))

fout.close()

#Do the Same thing for CDS
cds_windows_list = []
errors = []
for index, i in enumerate(cds_list):
    print 'Working on ' + i.get_name()
    print 'Index is: ' + str(index)
    #
    roi_sequence = i.get_sequence(genome_twobit)
    full_transcript = transcriptome_dict[i.attr['transcript_id']]
    transcript_sequence = full_transcript.get_sequence(genome_twobit)
    roi_start_in_transcript = transcript_sequence.find(roi_sequence)
    #Just a Very Quick Extra Check
    if roi_start_in_transcript != full_transcript.cds_start:
        errors.append(index)
    #
    if roi_start_in_transcript - 50 >= 0:
        lower_bound = roi_start_in_transcript - 50
    #
    else:
        lower_bound = 0 
    #
    if roi_start_in_transcript + 50 <= full_transcript.get_length():
        upper_bound = roi_start_in_transcript + 50
    #
    else:
        upper_bound = full_transcript.get_length()
    #
    subchain_window = full_transcript.get_subchain(lower_bound, upper_bound)
    subchain_window.attr['ID'] = i.attr['ID'] + '_' + str('100bpwindow') + '_' + str(roi_start_in_transcript-lower_bound) #overwrite the attribute to identify uORF uniquely
    cds_windows_list.append(subchain_window)

fout=open('./updated_quantitation/CDS_100bpwindows_aroundstart.fasta','w')
for i in cds_windows_list:
    fout.write(i.get_fasta(genome_twobit))

fout.close()

#Folks wanted the CDS Counts and uORF counts for everything from the CHX data. Import the Data 
CHX_1 = BAMGenomeArray(['path.to.reads'])
CHX_1.set_mapping(VariableFivePrimeMapFactory.from_file(open('path.to.psite'))) #doesn't read string like it should, so my workaround is to give it the open file handle instead

CHX_2 = BAMGenomeArray(['path.to.reads'])
CHX_2.set_mapping(VariableFivePrimeMapFactory.from_file(open('path.to.psite'))) #doesn't read string like it should, so my workaround is to give it the open file handle instead

CHX_3 = BAMGenomeArray(['path.to.reads'])
CHX_3.set_mapping(VariableFivePrimeMapFactory.from_file(open('path.to.psite')))

#Count uORF reads across 3 CHX replicates
uORF_counts_table=pd.DataFrame(index=range(0,len(uORF_list)),columns=['transcript', 'ID', 'Length', 'CHX_1_ct', 'CHX_2_ct', 'CHX_3_ct']) 
for index, feature in enumerate(uORF_list):
    print 'Working on ' + feature.get_name()
    print 'Index is: ' + str(index)
    uORF_counts_table.transcript[index] = feature.attr['transcript_id']
    uORF_counts_table.ID[index] = feature.get_name()
    uORF_counts_table.Length[index] = feature.get_length()
    uORF_counts_table.CHX_1_ct[index] = np.nansum(feature.get_counts(CHX_1))
    uORF_counts_table.CHX_2_ct[index] = np.nansum(feature.get_counts(CHX_2))
    uORF_counts_table.CHX_3_ct[index] = np.nansum(feature.get_counts(CHX_3))

#Save uORF feature count to file
uORF_counts_table.to_csv(path_or_buf='./updated_quantitation/uORF_feature_counts_Ingolia.tsv', sep='\t',header=True, index=False, index_label=False)

#Count CDS reads across 3 CHX replicates
CDS_counts_table=pd.DataFrame(index=range(0,len(cds_list)),columns=['transcript', 'ID', 'Length', 'CHX_1_ct', 'CHX_2_ct', 'CHX_3_ct']) 
for index, feature in enumerate(cds_list):
    print 'Working on ' + feature.get_name()
    print 'Index is: ' + str(index)
    CDS_counts_table.transcript[index] = feature.attr['transcript_id']
    CDS_counts_table.ID[index] = feature.get_name()
    CDS_counts_table.Length[index] = feature.get_length()
    CDS_counts_table.CHX_1_ct[index] = np.nansum(feature.get_counts(CHX_1))
    CDS_counts_table.CHX_2_ct[index] = np.nansum(feature.get_counts(CHX_2))
    CDS_counts_table.CHX_3_ct[index] = np.nansum(feature.get_counts(CHX_3))

#Save CDS feature count to file
CDS_counts_table = CDS_counts_table.drop('ID',axis=1) #don't need the ID column really
CDS_counts_table.to_csv(path_or_buf='./updated_quantitation/CDS_feature_counts_Ingolia.tsv', sep='\t',header=True, index=False, index_label=False)

#Count RNA reads in CDS and across the transcript
rnaseq = BAMGenomeArray(['path.to.reads'])
rnaseq.set_mapping(CenterMapFactory()) #seems reasonable- scales the read across length 1/L (l=read length) at each position in read

mRNA_counts_table=pd.DataFrame(index=range(0,len(cds_list)),columns=['transcript', 'Length_CDS', 'CDS_ct_mRNA', 'Len_Transcript', 'Transcript_ct_mRNA']) 
for index, feature in enumerate(cds_list):
    print 'Working on ' + feature.get_name()
    print 'Index is: ' + str(index)
    mRNA_counts_table.transcript[index] = feature.get_name()
    mRNA_counts_table.Length_CDS[index] = feature.get_length()
    mRNA_counts_table.CDS_ct_mRNA[index] = np.nansum(feature.get_counts(rnaseq))
    full_transcript = transcriptome_dict[feature.get_name()]
    mRNA_counts_table.Len_Transcript[index] = full_transcript.get_length()
    mRNA_counts_table.Transcript_ct_mRNA[index] = np.nansum(full_transcript.get_counts(rnaseq))

#Save the RNA counts to file
mRNA_counts_table.to_csv(path_or_buf='./updated_quantitation/mRNA_feature_counts_Ingolia.tsv', sep='\t',header=True, index=False, index_label=False)

#Write summary statistics to file
summary_table_ingolia=pd.DataFrame(index=range(0,4),columns=['dataset', 'mapped_reads']) 
for index, element in enumerate([CHX_1, CHX_2, CHX_3, rnaseq]):
    summary_table_ingolia.mapped_reads[index] = element.sum()

summary_table_ingolia.dataset = np.asarray(['CHX_1', 'CHX_2', 'CHX_3', 'RNAseq'])

summary_table_ingolia.to_csv(path_or_buf='./updated_quantitation/summary_ingolia_CHX_mRNA_datasets.tsv', sep='\t',header=True, index=False, index_label=False)

#Folks wanted the distance of uORF to CDS
distance_table = pd.DataFrame(index=range(0,len(uORF_list)),columns=['uORF_ID', 'dist_uORFstart_to_CDS', 'dist_uORFend_to_CDS', 'Len_uORF'])
start_stop_regex = re.compile(".*?_([0-9]+)_([0-9]+)")
for index, feature in enumerate(uORF_list):
    print 'Working on: ' + feature.get_name()
    print 'Index is: ' + str(index)
    distance_table.uORF_ID[index] = feature.get_name()
    
    #Extract Start and Stop of uORF in Transcript Coordinates
    regex_matches = start_stop_regex.match(feature.get_name())
    start_location, stop_location = [int(i) for i in regex_matches.groups()]
    
    full_transcript = transcriptome_dict[feature.attr['transcript_id']]
    distance_table.dist_uORFstart_to_CDS[index] = full_transcript.cds_start - start_location
    distance_table.dist_uORFend_to_CDS[index] = full_transcript.cds_start - stop_location
    distance_table.Len_uORF[index] = feature.get_length()

distance_table.to_csv(path_or_buf='./updated_quantitation/uORFs_distance_toCDSstart.tsv', sep='\t',header=True, index=False, index_label=False)
