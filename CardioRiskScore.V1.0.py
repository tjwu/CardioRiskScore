# Calculate risk score and make risk judgement by reading through sample's SNP vcf
# python CardioRiskScore.py 1000021992_HJGCTBCX2-1-IDMB4_snp.vcf
# Output 1000021992_HJGCTBCX2-1-IDMB4_RiskScore.csv which contains risk score and judgement
# Push to LIMS? 
import csv, os, sys

Infile = sys.argv[1]
VCFfile = csv.reader(open(Infile, "r"), delimiter= "\t")
#ScoreTable = csv.reader(open("RiskPosScore.csv", "r"), delimiter= "\t")
#ScoreTable = csv.reader(open("RiskPosScore.50.updateAug202018.csv", "r"), delimiter= "\t") # Look up first 50 rsIDs.
#ProxyScoreTable = csv.reader(open("RiskPosProxy.csv", "r"), delimiter= "\t") # Look up first 50 rsIDs.
OutputFN =  Infile.replace("vcf", "RiskScore.csv")
OutputFN =  os.path.basename(OutputFN)

W = csv.writer(open("CardioRisk/%s" % OutputFN, "w"), delimiter=",", lineterminator="\n")

# Compose output header line
def WriteHeader(OutFile):
    header = ["RSid", "Chr:Pos", "Ref", "Alt", "Risk Allele", "Genotype", "Ln(Published Odds Ratio)", "# of Risk Alleles", "# of Risk Alleles * Ln(OR)"]
    OutFile.writerow(header)

# Get variant position and risk score value
def RiskScore():
    ScoreTableFile = csv.reader(open("/hgsccl_software/devel/TJ/CardioRiskScore/RiskPosScore.50.updateAug202018.csv", "r"), delimiter= "\t") # Look up first 50 rsIDs.
    RSdict={}
    for j in ScoreTableFile:
        jRSid, jChr, jPos, jRef, jAlt, jRisk, jValue, jLnV = j[0:]    
        jkey = jChr + ":" + jPos
        RSdict[jkey] = [jRisk, jLnV, jRSid]
    return RSdict

def RStoScore():
    ScoreTableFile = csv.reader(open("/hgsccl_software/devel/TJ/CardioRiskScore/RiskPosScore.50.updateAug202018.csv", "r"), delimiter= "\t") # Look up first 50 rsIDs.
    RStoSdict={}
    for j in ScoreTableFile:
        jRSid, jChr, jPos, jRef, jAlt, jRisk, jValue, jLnV = j[0:]   
        RStoSdict[jRSid] = j
    return RStoSdict

def ProxyScore(ProxyScoreTableFile):
    PSGenedict = {}
    PSPosdict = {}
    for j in ProxyScoreTableFile:
        jRSid, jChr, jPos, jRef, jAlt, jRisk, jValue, jLnV, jGene = j[0:]    
        jkey = jChr + ":" + jPos
        #RSdict[jkey] = [jRisk, jLnV, jRSid]
        if PSPosdict.has_key(jkey):
            PSGenedict[jGene].append(jkey)
        else: 
            PSPosdict[jkey] = [jRisk, jLnV, jRSid]
        if PSGenedict.has_key(jGene):
            PSGenedict[jGene].append(jkey)
        else:
            PSGenedict[jGene] = [jkey]
    return PSGenedict, PSPosdict

# Use risk score and count of allele to calculate final score
def CalulateRiskScore(Count,Score):
    Final = float(Score) * Count
    return Final

# Risk category judgement Hight, Intermediate, Low. (More categories?)
def RiskJudgement(x):
    if x >=4.5824:
        RiskComment = "Top 5%"
    else:
        RiskComment = "Normal"
    #if x >= 4.326:
    #    RiskComment = "Top 20%"
    #if x >=3.96:
    #    RiskComment = "High Risk"
    #elif 3.96 > x >= 3.07:
    #    RiskComment = "Intermediate Risk"
    #else:
    #    RiskComment = "Low Risk"
    return RiskComment

# Create RSid-Gene dict to avoid count twice on same risk score row.
def RSidGene():
    RSidGenedict={}
    F = csv.reader(open("/hgsccl_software/devel/TJ/CardioRiskScore/Gene_RSid.csv","r"), delimiter = "\t")
    for line in F:
        Gene = line[0]
        for v in line[1:]:
            if v.find("rs")>=0:
                RSidGenedict[v] = Gene
    return RSidGenedict

# Create Gene-RSid dict to avoid count twice on same risk score row.
def GeneRSid():
    GeneRSiddict={}
    F = csv.reader(open("/hgsccl_software/devel/TJ/CardioRiskScore/Gene_RSid.csv","r"), delimiter = "\t")
    for line in F:
        Gene, RSid = line
        GeneRSiddict[Gene] = RSid
    return GeneRSiddict

# Read through sample's SNP vcf file, write position outcome and return final total count and total score
def CheckVCF(OutFile):
    VCFfile = csv.reader(open(Infile, "r"), delimiter= "\t")
    TotalCount = 0
    TotalScore = 0
    RSdict = RiskScore()
    RSidGenedict = RSidGene()
    DoneGenelist = []
    #PSGenedict, PSPosdict = ProxyScore(ProxyScoreTable)
    for i in VCFfile:
        if i[0].find("#") >=0:
            continue
        iChr, iPos, iID, iRef, iAlt, iQual, iFilter = i[0:7]
        Geno = i[-1].split(":")[0]
        #if iFilter.find("BLOCKAVG") >= 0:
        #    continue
        ikey = iChr + ":" + iPos
        # Read through vcf, get risk score
        if RSdict.has_key(ikey):
            Risk, LnV, RSid = RSdict[ikey]
            TempGene = RSidGenedict[RSid]
            if TempGene in DoneGenelist:
                continue
            else:
                DoneGenelist.append(TempGene)
            #if iAlt == Risk: modify to handle tri-allele
            #if Risk.find(iAlt) >= 0 or Risk.find(iRef) >= 0 :
            if iRef.find(Risk) >= 0:
                #Geno = i[-1].split(":")[0]
                if Geno == "0/0":
                    RiskAlleleCount = 2
                elif Geno == "0/1":
                    RiskAlleleCount = 1
                else:
                    RiskAlleleCount = 0
            elif iAlt.find(Risk) >= 0:
                if Geno == "1/1":
                    RiskAlleleCount = 2
                elif Geno == "0/1":
                    RiskAlleleCount = 1
                else:
                    RiskAlleleCount = 0
            else:
                RiskAlleleCount = 0

            Final = CalulateRiskScore(RiskAlleleCount, LnV)

            if iAlt.find(Risk) < 0 and iRef.find(Risk)< 0 and iAlt != ".":
                print "Check Ref, Alt and Risk", iRef, iAlt, Risk
                print i
            TotalScore += Final
            TotalCount += RiskAlleleCount

            Output = [RSid, "Chr"+ikey, iRef, iAlt, Risk, Geno, LnV, RiskAlleleCount, Final]
            OutFile.writerow(Output)
        # Lead SNP position is missing. Store proxy pos info
        #elif PSPosdict.has_key(ikey):
        #    Risk, LnV, RSid = PSPosdict[ikey]
        #    if Risk.find(iAlt) >= 0 or Risk.find(iRef) >= 0 :
        #        Geno = i[-1].split(":")[0]
        #        if Geno == "1/1" or Geno == "0/0":
        #            RiskAlleleCount = 2
        #            PSfinal = CalulateRiskScore(RiskAlleleCount, LnV)
        #            PSPosdict[ikey] = PSPosdict[ikey] + [RiskAlleleCount, PSfinal]
        #        else:
        #            RiskAlleleCount = 1
        #            PSfinal = CalulateRiskScore(RiskAlleleCount, LnV)
        #            PSPosdict[ikey] = PSPosdict[ikey] + [RiskAlleleCount, PSfinal]
        #    else:
        #        RiskAlleleCount = 0
        #        PSfinal = CalulateRiskScore(RiskAlleleCount, LnV)
        #        PSPosdict[ikey] = PSPosdict[ikey] + [RiskAlleleCount, PSfinal]
    #return [TotalCount, TotalScore, DoneGenelist, PSPosdict, PSGenedict]
    return [TotalCount, TotalScore, DoneGenelist]

# Write final outcome line to output file
def CommentTotalLine(TotalScore, TotalCount, OutFile):
    RiskComment = RiskJudgement(TotalScore)
    #RiskCommentTotal = ["Risk Comment", RiskComment, "", "", "", "Risk Score", TotalCount, round(TotalScore,3)]
    RiskCommentOutput1 = ["", "", "", "", "", "", "Risk Score", TotalCount, round(TotalScore,3)]
    RiskCommentOutput2 = ["", "", "", "", "", "", "","Risk Comment", RiskComment]
    W.writerow(RiskCommentOutput1)
    W.writerow(RiskCommentOutput2)

def CountMissing(MissingList, TotalScore, TotalCount, OutFile):
    GtoR = GeneRSid()
    RtoS = RStoScore()
    #VCFfile = csv.reader(open(Infile, "r"), delimiter= "\t")
    #print MissingList
    for m in MissingList:
        MissingRSid = GtoR[m]
        mRSid, mChr, mPos, mRef, mAlt, mRisk, mValue, mLnV = RtoS[MissingRSid]

        VCFfile = csv.reader(open(Infile, "r"), delimiter= "\t")
        for v in VCFfile:
            if v[0].find("#") >=0:
                continue
            vChr = v[0]
            if v[7].find("BLOCKAVG") >=0:
                vStart = v[1]
                vStop = v[7].split(";")[0].replace("END=","")
                if mChr == vChr and vStart <= mPos <= vStop:
                    if vStart == "1":
                        continue
                    DP = v[9].split(":")[3]
                    if int(DP) >= int(15):
                        if mRef == mRisk:
                            mRiskCount = 2
                        else:
                            mRiskCount = 0 
                        TotalScore += CalulateRiskScore(mRiskCount, mLnV)
                        TotalCount += mRiskCount
                        #print mRiskCount, TotalScore, m, mRSid
    FinalLine = ["", "", "", "", "", "", "Risk Score", TotalCount, round(TotalScore,3)]
    OutFile.writerow(FinalLine)

# If any gene is missing, report here ( working on RSid)
#def ReportMissingGene(DoneGenelist, OutFile, PSGenedict, PSPosdict, TotalScore, TotalCount):
def ReportMissingGene(DoneGenelist, OutFile, TotalScore, TotalCount):
    Genelist = ['ABCG8','ABO','ADAMTS7_1','ADAMTS7_2','ANKS1A','APOA5','APOB','BCAP29','CDKN2A','CDKN2BAS','COL4A1','COL4A1/COL4A2','CXCL12_1','CXCL12_2','CYP17A1','EDNRA','FLT1','FURIN/FES','GGCX/VAMP8','GUCY1A3','HDAC9','HHIPL1','HNF1A','IL6R','KCNE2','KCNK5','KIAA1462','LDLR','LIPA','LPA_1','LPA_2','MIA3','MRAS','PCSK9','PDGFD','PHACTR1','PLG','PPAP2B','RASD1','SH2B3','SLC22A3/LPAL2/LPA','SLC22A4/SLC22A5','SMG6','SORT1','TCF21','TRIB1','UBE2Z','WDR12','ZC3HC1','ZEB2--AC074093.1']
    cc = 0
    Missing = list(set(Genelist) - set(DoneGenelist))
    Proxylist = []
    # Report Missing Lead SNP 
    if len(Missing) != 0:
        Newline = ["Missing Lead SNP in these Gene(s):"] + Missing
        OutFile.writerow(Newline)
        CountMissing(Missing, TotalScore, TotalCount, OutFile)

    # Use Missing list and previous stored Proxy SNP info to get final score and comment
    #for m in Missing:
    #    pScore = 0
    #    pCount = 0
    #    if m in Proxylist:
    #        continue
    #    for pos in PSGenedict[m]:
    #        if PSPosdict.has_key(pos):
    #            pLnV, pRSid, pCount, pScore = PSPosdict[pos][1:]
    #            TotalScore += pScore
    #            TotalCount += pCount
    #            newrow = [pRSid, "Chr"+pos, pLnV, pCount, pScore]
    #            OutFile.writerow(newrow)
    #            Proxylist.append(m)
    ## Reoprt using proxy here. 
    #if len(Proxylist) !=0:
    #    ProxyFinalList = ["Gene(s) proxy SNP:"] + Proxylist
    #    OutFile.writerow(ProxyFinalList)
    #    CommentScore  = ["Final Comment", RiskJudgement(TotalScore), "Final Risk Score", TotalScore]
    #    OutFile.writerow(CommentScore)

def Run():
    WriteHeader(W)
    #TotalCount, TotalScore, DoneGenelist, PSPosdict, PSGenedict = CheckVCF(W)
    TotalCount, TotalScore, DoneGenelist = CheckVCF(W)
    CommentTotalLine(TotalScore, TotalCount, W)
    ReportMissingGene(DoneGenelist, W, TotalScore, TotalCount)
Run()
