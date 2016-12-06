********************************************************************************
Motifx.py 

    Automate submitting jobs to Motif-x website.

    Input  : Plain text file listing excel files to process, one excel file name per line.

    Excel file format:
    
        Ppep    Group   Localized_Sequence  Motif_X_Input_Peptide
        YGL076C_T8_S11  Induced AAEKILtPEsQLKK  AAEKILT*PES*QLKK

        Column order is unimportant, column names must match above.

    usage: Motifx.py -f inputfiles 

        Optional Arguments:
        -o Minimum number of times each of your extracted motifs to occur in the data set (10)
        -s P-value threshold for the binomial probability (.000001)
        -u upload a newer version of SGD proteome (orf_trans.fasta) than used by motifx.
        
            If you choose to use the default from Motifx, then the follow on motifx matching
            may miss a few motifx, as Motifx uses an older version of SGD proteome.
    
        -w Number of total characters in motif, (13)

    Output : Table with amino acid sequence, motif and png logo images plus log files for each input Excel file.


    NOTE: Motifx.py automatically calls exactFastaMatch.py, so it does not need to be 
          run separately.

    To see Python Docs and get a better explaination of the program:

    Open python console and enter
    import sys
    sys.path.append('/home/GLBRCORG/mplace/scripts/motifx')
    import Motifx
    help(Motifx)

********************************************************************************
Alternate ORF fasta file from SGD R64-2-1 is located:    

    /home/GLBRCORG/mplace/scripts/motifx/orf_trans_all.20150113.fasta  
    
********************************************************************************
exactFastaMatch.py

    Purpose: Exact match a short sequence to a gene in fasta file  

    Required parameters:
    
        -f fasta file, expected to be an orf file.
           where the gene name is right after the >
        
            example:  >YAL001C 
            
        -g gene file, one gene name per line.
           Used to reduce the fasta search space.
           one gene name per line, like:
           YAL001C
           YAL002W

        -o output file name
        -s sequence file, short sequences one per line.
           one short sequence per line, expect a CSV file
           assumes file looks like Motifx.py results:

           PRARSSSVSNAAL,Induced,...R..S.S....
           LRERSRSNSSALA,Induced,...R..S.S....
           GSERTRSISFSKL,Induced,...R..S.S....
           SRSRSRSKSNANA,Induced,...R..S.S....

           
    output file: Matched gene names in text file
        YPL137C PRARSSSVSNAAL   ...R..S.S....   Induced
        YOL036W LRERSRSNSSALA   ...R..S.S....   Induced
        YNR047W GSERTRSISFSKL   ...R..S.S....   Induced
        YBL007C SRSRSRSKSNANA   ...R..S.S....   Induced

********************************************************************************

