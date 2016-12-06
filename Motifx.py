#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program: motifx.py

Purpose: Automate submitting jobs to Motif-x website. http://motif-x.med.harvard.edu/motif-x.html
         Each job will be run with 3 central characters S, T, Y.  This means each input file will
         have 3 different central characters used.  User also has the option of matching the 
         output peptides to the yeast proteome, exact match only.  This means if you are
         using anything other than S288C peptides for motifx, nothing will match.  
         
         NOTE: The default SGD proteome file used by motifx is different (older) 
         than the current SGD file.  This results in a few motifx which cannot be 
         mapped back to identify the gene.  If you don't care just use the default
         otherwise use -u to upload the latest SGD orf_trans fasta file.  

Input : One or more peptide excel files, column order is unimportant but column name must match, 
        each file looks like:

Ppep    Group    Localized_Sequence    Motif_X_Input_Peptide
YGL076C_T8_S11    Induced    AAEKILtPEsQLKK    AAEKILT*PES*QLKK
YNR047W_T428    Induced    AASEPNGLQLASATSPtSSSAR    AASEPNGLQLASATSPT*SSSAR

Output : A log file, results file, directories for each central character containing
         logo images and motif results page in a log file.

Mofif-x input parameters: ms/ms Extend from: SGD Yeast Proteome Central Character Occurances = 10

Dependencies: Python 3 
              Python modules: argparse, BeautifulSoup4, requests, xlrd, unicodecsv
              exactFastaMatch.py (matches motif to gene name)

author: Mike Place
Date:   8/19/2016
"""
import argparse	                # handle command line args
from bs4 import BeautifulSoup       # html parser
from collections import OrderedDict 
from collections import defaultdict
import exactFastaMatch as efm       # written by Mike Place to match motifx motifs to gene names
import os
import re                           # regex
import requests                     # HTTP library
import shutil                       # 
import sys
import time                         # sleep
from tqdm import tqdm               # progress bar
import unicodecsv                   # handle unicode chars from excel file
import xlrd                         # handle excel files

class Motifx ( object ):
    """
    Methods and data structures for Motifx class 
    """
    def __init__(self, file, occ, sig, prot, wid ):
        """
        Set up Motifx object
        text        : text version of input excel file
        pep         : list of peptides to submit to motifx website
        occurrence  : occurrences, default is 10
        sig         : significance threshold, default = 0.000001
        proteome    : Yeast ORF fasta file to use, default is the motifx website default
        width       : width of motif, default = 13       
        """
        self.occurance   = occ
        self.sig         = sig
        self.proteome    = prot
        self.width       = wid
        self.dir         = os.getcwd()
        self.text        = self.excelToText(file)
        self.centralRes  = ['S','T','Y']       # central character on motif-x web form required by motifx site
        self.fileName    = re.sub(r'.xlsx', '', file)     # filename 
        self.group       = ''
        self.logo        = []                  # list of locations for logo images
        self.pep         = defaultdict(list)   # dict of lists, key = YNR047W_T428, value = list with the peptide sequence [0], group [1] positions
        self.geneList    = set()               # unique list of the genes, used to pare down the yeast proteome for matching.
        self.getPep()
        self.peptideFile = 'pepFile.txt'       # temp peptide file name
        self.result      = OrderedDict()
        self.writeGeneList()                   # write gene list to file for downstream matching of motif to orf names

    def getPep(self):
        """
        Get the list of amino acid sequences to submit to Motifx website. 
        Load self.pep dictionary
        key = geneName_x   like YNR047W_T428
        value is a list , position[0] = peptide, position[1] = group 
        
        'YNL243W_S284': ['KREPS*VTPAR', 'Induced']        
        """
        with open(self.text, 'r') as f:
            for row in f:
                row = row.rstrip()
                if row.startswith('Ppep'):             # this should be the header
                    headers  = row.split(',')
                    pepCol   = headers.index("Motif_X_Input_Peptide")    # get the peptide column, makes column order unimportant
                    geneCol  = headers.index("Ppep")                     # get the gene name looks like: YNR047W_T428
                    groupCol = headers.index("Group")                    # get group column index 
                else:
                    data = row.split(',')
                    self.pep[data[geneCol]].append(data[pepCol])
                    self.pep[data[geneCol]].append(data[groupCol])
                    geneName = data[geneCol].split('_')                  # split gene name info, YER178W_Y309
                    self.geneList.add(geneName[0])                    # add gene names to list to limit proteome matching space
        
        # get group type i.e. induced/repressed
        self.group = next(iter (self.pep.values()))[1]
    
    def writeGeneList(self):
        """
        Write unique gene list to file for down stream matching of orf sequence
        """
        with open('GeneList.txt', 'w') as gOut:
            for gene in self.geneList:
                gOut.write('%s\n' %(gene))
    
    def pepFile(self):
        """
        Write peptide to file to use when submitting job to Motifx website.
        This is just a temporary file, will be overwritten if multiple files are processed.
        The last processed file's peptides will be in this file.
        """
        with open('pepFile.txt', 'w') as pf:
            for key in self.pep.keys():
                pf.write("%s\n" %(self.pep[key][0]))
        pf.close()

    def excelToText( self, file ):
        """
        Converts an Excel file to a CSV file.
        If the excel file has multiple worksheets, only the first worksheet is converted.
        Uses unicodecsv, so it will handle Unicode characters.
        Uses a recent version of xlrd, so it should handle old .xls and new .xlsx equally well.
        """        
        wb = xlrd.open_workbook(file)
        sh = wb.sheet_by_index(0)

        fh = open('Temp_file.txt',"wb")
        csv_out = unicodecsv.writer(fh, encoding='utf-8')
        
        # go through the file line by line
        for row_number in range (sh.nrows):
            csv_out.writerow(sh.row_values(row_number))

        fh.close()
        return 'Temp_file.txt'
        
    def submitMotifX( self, char ):
        """
        Submit protein list to the Motif-x website.
                
        url:  http://motif-x.med.harvard.edu/cgi-bin/multimotif-x.pl
        
        result url:
        http://motif-x.med.harvard.edu/  plus  /cgi-bin/jobres.pl?jobid=20160830-8155-04402686 (example only)
        which is found in tag:  <A HREF="/cgi-bin/jobres.pl?jobid=20160830-8155-04402686" TARGET="_blank">Check results</A>
        
        Inputs to the web form:     
            form text     -- object name
            upload file   -- fgfile 
            ms/ms         -- fgtype   
            extend from   -- fgextenddb
            central char  -- fgcentralres
            width         -- width
            occurances    -- occurrences
            background    -- bgdb
            select get_motifs button
            
            if using a user provided proteome fasta file:
            bgdb          -- background database
            bgtype        -- type fasta in this case
            bgcentralres  -- central Character
            bgextenedb    -- extention database (required, even though the user provide proteome is used.)
            bgfile        -- user provided proteome fasta            
            
        return the url of the result page
        """
        url       =  'http://motif-x.med.harvard.edu/cgi-bin/multimotif-x.pl'           # submission url
        
        # set values required for the form submission  'submit':'get
        # Here we set the values for using the default SGD proteome provided by motifx
        if self.proteome == 'default':
            form_data = {  'fgtype' : 'ms', 'fgextenddb': 'SGD_Yeast.fasta', 'fgcentralres' : char, 'width' : self.width,
                     'occurrences' : self.occurance, 'significance': self.sig, 'bgdb' : 'SGD_Yeast.fasta'  } 
            file     = { 'fgfile': open( self.dir + '/' + self.peptideFile, 'rb') }          # peptide upload file, one peptide per line
        # Set values for using the user provided proteome
        else:            
            form_data = {  'fgtype' : 'ms', 'fgcentralres' : char, 'fgextenddb': 'SGD_Yeast.fasta', 'width' : self.width,
                     'occurrences' : self.occurance, 'significance': self.sig, 'bgdb' : 'uploaded', 
                     'bgtype' : 'fasta', 'bgcentralres' : char, 'bgextenddb' : 'SGD_Yeast.fasta' }  
            # Now we have 2 files to submit to motifx
            file     = { 'fgfile': open( self.dir + '/' + self.peptideFile, 'rb'),           # peptide upload file, one peptide per line
                         'bgfile': open( self.dir + '/' + self.proteome, 'rb')   }           # user provided proteome fasta
                
        # send value to website, post will time out after 6 min seconds with no response 
        try:
            response  = requests.post( url, data = form_data, files = file, allow_redirects = True, timeout = 360 )
            soup      = BeautifulSoup( response.text, 'html.parser')     # store resulting webpage as text, need to parse to get result page
            tag       = str(soup.a)
            jobID     = re.findall('"([^"]*)"', tag )       # get the results page location jobID[0]
        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)
        
        time.sleep(120)                                  # give motifx time to run the job
  
        # return the results location
        return 'http://motif-x.med.harvard.edu/' + jobID[0] 
                  
    def parseResults( self, resultURL, char ):
        """
        Parse the results page from a Motif-x job, resultURL is the results webpage.
        Write the html results page to a file.
        Download all the logo images to current directory.
        """
        centralChar = char
        # make directory for each central character's results
        resultsDir = self.fileName + '_' + centralChar        
        os.mkdir( resultsDir )
        
        # get the results page
        try:
            response = requests.get(resultURL, timeout = 60)          # html results 
        except requests.exceptions.RequestException as e:             
            print(e)
            sys.exit(1)
            
            
        # write the html page to log file.
        with open( resultsDir + '/' + self.fileName + '-' + centralChar + '.log', 'w') as log:
            log.write(resultURL)
            log.write(response.text)
       
        # parse the returned result html
        tree = BeautifulSoup( response.text, "lxml" )
        # get the body of the page
        data = tree.body.find_all('font')[3]
        
        # find motif asscociated with peptide
        for m in data.findNextSiblings('a'):
            if m.text not in self.result:          # add key to dictionary, key is the motif ex: ".....T....S"
                self.result[m.text] = []
            for i in m.children:
                self.result[m.text].append(str(i.next))         # append a string will all the motifs combined
      
        # Get logo png files
        for logo in tree.find_all('a'):
            for img in logo.find_all('img'):
                path = 'http://motif-x.med.harvard.edu/' + img['src']      # get web address for logo image
                req = requests.get(path, stream = True)                    # download image
                if req.status_code == 200:       
                    imageFile = re.sub(r'/logos', '', img['src']) 
                    # good to go, now download
                    with open( resultsDir + imageFile, 'wb') as png:
                        req.raw.decode_content = True
                        shutil.copyfileobj(req.raw, png)
                else:                                                      # if unable to download logo image record path
                    with open("logo-image.log", 'a') as log:
                        log.write("Unable to download image: %s" %(path))              
            
    def writeResults( self ):
        """
        Write peptide table, 3 columns, comma separated.
        peptide,Group,motif
        """
        with open(self.fileName + '-Motifx-results.txt', 'w') as out:
            for k,v in self.result.items():
                motifs = v[0].split('\n')
                [ motifs.pop(0) for x in range(2)]
                [ motifs.pop() for x in range(3)]
                for row in motifs:
                    line = "".join(row)                           # this joins the peptides into a single string
                    line = line + "," + self.group + "," + k
                    out.write("%s\n" %(line))
                
def main():
    """
    Process command line arguments and run program.
    """
    cmdparser = argparse.ArgumentParser(description="Automate submitting jobs to Motif-x website.",
                                        usage='%(prog)s -f <File listing excel files>  ', prog='Motifx.py'  )                                  
    cmdparser.add_argument('-f', '--file', action='store', dest='FILE', help='file listing excel files to process expected', metavar='')
    cmdparser.add_argument('-o', '--occurrence', action='store', dest='OCC', help='occurrences, default is 10', metavar='')
    cmdparser.add_argument('-s', '--sig',  action='store', dest='SIG', help='significance, default = 0.000001', metavar='', type=float)
    cmdparser.add_argument('-u', '--upload', action='store', dest='UPLOAD', help='Upload Yeast ORF fasta file to use', metavar='')
    cmdparser.add_argument('-w', '--width', action='store', dest='WIDTH', help='width, default = 13', metavar='', type=int)      
    cmdparser.add_argument('-i', '--info', action='store_true', dest='INFO', help='Detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
             
    if cmdResults['INFO']:
        print("\n  Motifx.py ")
        print("\n  Automate submitting jobs to Motif-x website.")
        print("\n  Input  : Plain text file listing excel files to process, one excel file name per line.")
        print("\n\tExcel file format:")
        print("\tPpep	Group	Localized_Sequence	Motif_X_Input_Peptide")
        print("\tYGL076C_T8_S11	Induced	AAEKILtPEsQLKK	AAEKILT*PES*QLKK")
        print()
        print("\tColumn order is unimportant, column names must match above.")       
        print()
        print(" usage: Motifx.py -f inputfiles ")
        print()
        print(" Optional Arguments:")
        print("\t-o Minimum number of times each of your extracted motifs to occur in the data set (10)")
        print("\t-s P-value threshold for the binomial probability (.000001)")
        print("\t-u upload a newer version of SGD proteome (orf_trans.fasta) than used by motifx.")
        print("\t   If you choose to use the default from Motifx, then the follow on motifx matching")
        print("\t   may miss a few motifx, as Motifx uses an older version of SGD proteome.")
        print("\t-w Number of total characters in motif, (13)")
        print("\n  Output : Table with amino acid sequence, motif and png logo images plus log files for each input Excel file.") 
        print("\n\tTo see Python Docs and get a better explaination of the program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport Motifx")
        print("\thelp(Motifx)")
        print("\n\tSee Mike Place for help or suggestions.\n")
        sys.exit(1)
    
    # Get the inputfile, which contains a list of file names to process
    if cmdResults['FILE']:
        userFile = cmdResults['FILE']
    else:
        print('')
        cmdparser.print_help()
        sys.exit(1)
        
    # check for the occurance parameter 
    if cmdResults['OCC']:
        occurance = cmdResults['OCC']
        if not occurance.isdigit():
            print('\n\t-o occurance is a number, default = 10 \n')
            print('\tThe occurrence threshold refers to the minimum number of times'
                  '\n\tyou wish each of your extracted motifs to occur in the data set.\n')
            cmdparser.print_help()
            sys.exit(1)
    else:
        occurance = 10
            
    # check for significance parameter
    if cmdResults['SIG']:
        sig = float(cmdResults['SIG'])
        if sig > .1:
            print('\n\tYou are using a significance threshold greater than .1\n')
        elif sig == 0.0:
            print('\n\tYou are using a significance threshold of zero.\n')
    else:
        sig = '{:f}'.format(0.000001)
        
    # check for an alternate proteome file
    if cmdResults['UPLOAD']:
        proteome = cmdResults['UPLOAD']
        if not os.path.exists(proteome):
            print('\n\tAlternate proteome file does not exist.\n')
            cmdparser.print_help()
            sys.exit(1)
    else:
        proteome = 'default'
        
    # check of an alternate window width
    if cmdResults['WIDTH']:
        width = cmdResults['WIDTH']
    else:
        width = 13
        
        
    # open input file & process each file listed
    with open(userFile, 'r') as f:
        for line in f:
            line = line.rstrip()
            print("Processing file: %s" %(line))     # prints the name of the file currently being processed
            d = Motifx(line, occurance, sig, proteome, width )                         # create a Mofifx object 
            d.pepFile()                              # write peptide to temp file for use in job submission.            
            # submit job to motifx website
            for char in tqdm(d.centralRes):                # loop through all potential central characters
                resultsPage = d.submitMotifX(char)
                d.parseResults(resultsPage, char)

            d.writeResults()                         # write out results to file
            
            matchName = re.sub(r'.xlsx','-matched.txt', line)    # create matched motif and gene name output file name
            data = efm.exactFastaMatch(proteome, d.fileName + '-Motifx-results.txt', 'GeneList.txt')   # create efm object
            data.matchSeq()
            data.writeResult(matchName)
            
    
if __name__ == "__main__":
    main()
