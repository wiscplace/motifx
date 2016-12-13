********************************************************************************
  Full example in /home/GLBRCORG/mplace/scripts/motifx 
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
        -u upload your own version of SGD proteome (orf_trans.fasta).
           
           newer orf fasta located: /home/GLBRCORG/mplace/scripts/motifx/orf_trans_all.20150113.fasta
        -w Number of total characters in motif, (13)

    Output :
        A final text table named after the input file containing all the motifs matched to a gene.
        Given an input file named,  motifx_sample.xlsx the final results file
        will be: motifx_sample-Motifx-results.txt 

        The other results are put in a directory.
        For instance if your input file is called motifx_sample.xlsx

        3 directories will be created one for each central character:

           motifx_sample_T    motifx_sample_S    motifx_sample_Y

        These contain the LOGO pngs and the original html results page.

    To see Python Docs and get a better explaination of the program:

        Open python console and enter
        import sys
        sys.path.append('/home/GLBRCORG/mplace/scripts/motifx')
        import Motifx
        help(Motifx)

********************************************************************************
MotifxPreAlign.py

    Purpose: Pre-align peptides to an orf fasta file for use with Motifx website.

    Required parameters:

        -p peptide info file,

    Example, one per line of the following:

        YBR162W-A_T5,Induced,AVQtPR,AVQT*PR 
        YGL076C_T8_S11,Induced,AAEKILtPEsQLKK,AAEKILT*PES*QLKK

    Optional parameters:

        -f fasta file, expected to be an orf file.
        -w width of peptides, all peptides will be the same length.

    Outut: A list of peptides centered on the phosphorylated amino acid.
           List is meant to be used as input for the motifx website.

    To see Python Docs and get a better explaination of the program:

        Open python console and enter
        import sys
        sys.path.append('/home/GLBRCORG/mplace/scripts/motifx')
        import MotifxPreAlign
        help(MotifxPreAlign)

********************************************************************************

